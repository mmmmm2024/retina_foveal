#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一次変性時におけるAIIAC膜電位波形からPSDを計算
5–15 HzのバンドパワーをRod(%)×Cone(%)の行列にしてCSVとヒートマップ(PDF)で出力

処理:
- 1000–6000 ms を抽出（CROP_MS）
- FFT/片側PSD（Hann窓あり/なし切替、線形デトレンド任意）
- 5–15 Hz を積分してバンドパワー算出
- 11×11(Rod×Cone)の行列に配置し、色バー下限0で描画（X軸反転も可）

出力:
- AIIAC_bandpower_matrix_5-15Hz.csv
- AIIAC_bandpower_matrix_5-15Hz.pdf
"""

from __future__ import annotations
from pathlib import Path
import argparse
import re
import math
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

TARGET = "AIIAC"
# ======================= USER CONFIG (edit once) ============================
ROOT_DIR: str = f"Dim_{TARGET}_fovea"                 # "" → use this script's folder; or set absolute path
PATTERN: str = f"{TARGET}_R*_C*.txt"                  # filename glob pattern
RECURSIVE: bool = False                              # also search subfolders
# Analysis window & band (ms and Hz)
CROP_MS: Tuple[float, float] = (1000.0, 6000.0)
BAND_HZ: Tuple[float, float] = (5.0, 15.0)
DETREND_LINEAR: bool = False                          # optional linear detrend (after mean removal)
USE_HANN: bool = True
# Rod/ Cone grids (fixed to 11×11 as per project)
ROD_COUNTS: List[int] = [1, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400]
CONE_LEVELS: List[int] = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
# Plot controls
INVERT_X: bool = True                                 # True: X axis 100→0 (to match past figures)
FIGSIZE: Tuple[float, float] = (6.8, 6.2)
DPI: int = 150
# Colorbar controls
CB_CMAP: str = "Blues"                                # blue shades
CB_VMIN: float = 0.0                                  # start colorbar at 0 (band power ≥ 0 by definition)
CB_VMAX: Optional[float] = 0.002                       # set to a number to fix the upper bound; None → auto
# Output names
OUT_CSV: str = "AIIAC_bandpower_matrix_5-15Hz.csv"
OUT_PNG: str = "AIIAC_bandpower_matrix_5-15Hz.pdf"
# ===========================================================================

FILENAME_RE = re.compile(rf"^{TARGET}_R(\d+)_C(\d+)\.txt$", re.IGNORECASE)

# ----------------------- helpers (NEW) -----------------------
def trapz_compat(y, x=None, dx=1.0, axis=-1) -> float:
    """NumPyのバージョン差を吸収する台形積分。"""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x=x, dx=dx, axis=axis))
    if hasattr(np, "trapz"):
        return float(np.trapz(y, x=x, dx=dx, axis=axis))
    # ここには通常来ません
    raise AttributeError("Neither numpy.trapezoid nor numpy.trapz is available.")


def format_sig_fixed(x: float, sig: int = 2) -> str:
    """
    x を「有効数字 sig 桁」で固定小数（科学表記なし）で返す。
    例(sig=2):
      2      -> "2.0"
      3      -> "3.0"
      0.1    -> "0.10"
      0.02   -> "0.020"
      25     -> "25"
    """
    if not np.isfinite(x):
        return ""
    if x == 0:
        return "0.00"  # 0 は見た目安定のため

    ax = abs(x)
    exp = math.floor(math.log10(ax))
    decimals = max(sig - int(exp) - 1, 0)
    return f"{x:.{decimals}f}"


def choose_colorbar_scale(vmax: float) -> Tuple[float, int, str]:
    """
    vmax < 1 の場合：exp=floor(log10(vmax)) を採用して
      表示値 = 元の値 / 10^exp
      ラベル = Band power 5–15 Hz (10^{exp} mV^2)
    vmax >= 1 の場合：スケールなし（mV^2）
    """
    if not np.isfinite(vmax) or vmax <= 0:
        return 1.0, 0, r"Band power 5–15 Hz (mV$^2$)"

    if vmax >= 1.0:
        return 1.0, 0, r"Band power 5–15 Hz (mV$^2$)"

    exp = int(math.floor(math.log10(vmax)))   # 0.02 -> -2, 0.002 -> -3, 0.0002 -> -4
    scale = 10.0 ** exp                       # 10^-2 = 0.01
    label = rf"Band power 5–15 Hz (10$^{{{exp}}}$ mV$^2$)"
    return scale, exp, label

# ----------------------- original functions -----------------------
def linear_detrend(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - (A @ coef)


def periodogram_psd_onesided(x: np.ndarray, fs: float, use_hann: bool = True, remove_mean: bool = True):
    N = len(x)
    if remove_mean:
        x = x - np.mean(x)
    w = np.hanning(N) if use_hann else np.ones(N)
    xw = x * w
    U = (w**2).sum() / N
    X = np.fft.rfft(xw)
    Pxx_two = (np.abs(X) ** 2) / (fs * N * U)
    Pxx = Pxx_two.copy()
    if N % 2 == 0:
        Pxx[1:-1] *= 2.0
    else:
        Pxx[1:] *= 2.0
    f = np.fft.rfftfreq(N, d=1 / fs)
    return f, Pxx


def integrate_band(f: np.ndarray, Pxx: np.ndarray, fmin: float, fmax: float) -> float:
    if fmax <= fmin:
        return 0.0
    idx = (f >= fmin) & (f <= fmax)
    if not np.any(idx):
        return 0.0
    return float(trapz_compat(Pxx[idx], f[idx]))


def load_aiiac(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError(f"File {path} does not have two numeric columns.")
    return data[:, 0], data[:, 1]


def pick_folder(cli_folder: Optional[Path]) -> Path:
    if cli_folder is not None:
        return cli_folder
    if ROOT_DIR.strip():
        return Path(ROOT_DIR)
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()


def iter_files(folder: Path, pattern: str, recursive: bool) -> List[Path]:
    it = folder.rglob(pattern) if recursive else folder.glob(pattern)
    return sorted([p for p in it if p.is_file()], key=lambda p: p.name)


def rc_from_name(path: Path) -> Optional[Tuple[int, int]]:
    m = FILENAME_RE.match(path.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def compute_bandpower_for_file(path: Path, crop_ms: Tuple[float, float], band_hz: Tuple[float, float], detrend_linear: bool, use_hann: bool) -> float:
    t_ms, v_mV = load_aiiac(path)
    tmin, tmax = crop_ms
    mask = (t_ms >= tmin) & (t_ms <= tmax)
    if not np.any(mask):
        raise RuntimeError(f"No samples in crop range {tmin}-{tmax} ms for {path.name}")
    t = t_ms[mask]
    v = v_mV[mask]
    dt_ms = np.median(np.diff(t))
    fs = 1000.0 / dt_ms
    if detrend_linear:
        v = linear_detrend(v, t * 1e-3)
    f, Pxx = periodogram_psd_onesided(v, fs, use_hann=use_hann, remove_mean=True)
    return integrate_band(f, Pxx, band_hz[0], band_hz[1])


def build_matrix(powers: Dict[Tuple[int, int], float], rod_counts: List[int], cone_levels: List[int]) -> np.ndarray:
    Z = np.full((len(rod_counts), len(cone_levels)), np.nan, dtype=float)
    for i, R in enumerate(rod_counts):
        for j, C in enumerate(cone_levels):
            val = powers.get((R, C))
            if val is not None:
                Z[i, j] = val
    return Z


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Compute 5–15 Hz band power and plot a blue heatmap with colorbar min=0.")
    ap.add_argument("folder", nargs="?", type=Path, default=None, help="Folder with AIIAC_R*_C*.txt (default: ROOT_DIR or script folder)")
    ap.add_argument("--pattern", default=PATTERN)
    ap.add_argument("--recursive", action="store_true", default=RECURSIVE)
    ap.add_argument("--crop", nargs=2, type=float, metavar=("MS_MIN", "MS_MAX"), default=None)
    ap.add_argument("--band", nargs=2, type=float, metavar=("HZ_MIN", "HZ_MAX"), default=None)
    ap.add_argument("--detrend-linear", action="store_true", default=DETREND_LINEAR)
    ap.add_argument("--no-hann", action="store_true")

    # ★元の挙動は維持しつつ、OFFも可能にする（既定値は INVERT_X）
    ap.add_argument("--invert-x", dest="invert_x", action="store_true", help="X axis: 100→0 (right-to-left)")
    ap.add_argument("--no-invert-x", dest="invert_x", action="store_false", help="X axis: 0→100 (left-to-right)")
    ap.set_defaults(invert_x=INVERT_X)

    ap.add_argument("--out-csv", default=OUT_CSV)
    ap.add_argument("--out-png", default=OUT_PNG)

    # ★カラーバー上限もCLIで変えられるように（未指定ならCB_VMAX）
    ap.add_argument("--cb-vmax", type=float, default=CB_VMAX)

    args = ap.parse_args(argv)

    folder = pick_folder(args.folder)
    if not folder.exists() or not folder.is_dir():
        print(f"Error: {folder} is not a directory")
        return 2

    pattern = args.pattern or PATTERN
    recursive = bool(args.recursive)
    crop_ms = tuple(args.crop) if args.crop else CROP_MS
    band_hz = tuple(args.band) if args.band else BAND_HZ
    detrend_lin = bool(args.detrend_linear)
    use_hann = not bool(args.no_hann)
    invert_x = bool(args.invert_x)
    cb_vmax = args.cb_vmax  # Noneなら auto 扱いにしたい場合は --cb-vmax を外す

    files = iter_files(folder, pattern, recursive)
    if not files:
        print(f"No files matched pattern '{pattern}' in {folder}")
        return 0

    # Compute band power per file
    powers: Dict[Tuple[int, int], float] = {}
    for p in files:
        rc = rc_from_name(p)
        if rc is None:
            continue
        try:
            bp = compute_bandpower_for_file(p, crop_ms, band_hz, detrend_lin, use_hann)
        except Exception as e:
            print(f"[SKIP] {p.name}: {e}")
            continue
        powers[rc] = bp

    # Build matrix (rows=rod ascending, cols=cone ascending)
    Z = build_matrix(powers, ROD_COUNTS, CONE_LEVELS)

    # Export CSV (columns as Cone %, rows as Rod %)
    rod_pct = [r / max(ROD_COUNTS) * 100 for r in ROD_COUNTS]
    cone_pct_asc = [c / max(CONE_LEVELS) * 100 if max(CONE_LEVELS) > 0 else 0 for c in CONE_LEVELS]
    df = pd.DataFrame(Z, index=[f"{rp:.2f}" for rp in rod_pct], columns=[f"{cp:.0f}" for cp in cone_pct_asc])
    df.index.name = "Rod_%"
    df.to_csv(folder / args.out_csv)

    # Plot heatmap
    Zm = np.ma.masked_invalid(Z)

    # ---- NEW: decide exponent scaling based on vmax_used ----
    if cb_vmax is not None:
        vmax_used = float(cb_vmax)
    else:
        m = np.nanmax(Z)
        vmax_used = float(m) if np.isfinite(m) else 0.0

    scale, exp, cbar_label = choose_colorbar_scale(vmax_used)

    def cbar_fmt(val, pos):
        # show scaled value with 2 significant figures, no scientific notation
        return format_sig_fixed(val / scale, sig=2)
    # --------------------------------------------------------

    # ★ heatmap(正方形) と colorbar の高さを揃える
    CB_WIDTH_RATIO = 0.05

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1.0, CB_WIDTH_RATIO])

    ax  = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(
        Zm,
        origin="lower",
        extent=(-5, 105, -5, 105),
        interpolation="nearest",
        cmap=CB_CMAP,
        vmin=CB_VMIN,
        vmax=vmax_used if cb_vmax is not None else None,
        aspect="equal",
    )

    try:
        ax.set_box_aspect(1)
    except Exception:
        pass

    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbar_label)

    # ★重要：formatter を cb.formatter に直接入れて update_ticks
    cb.formatter = mticker.FuncFormatter(cbar_fmt)
    cb.update_ticks()

    ax.set_xlabel("Cone survival rate (%)")
    ax.set_ylabel("Rod survival rate (%)")
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_yticks(np.arange(0, 101, 10))
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)

    if invert_x:
        ax.invert_xaxis()

    out_fig = folder / args.out_png
    # サイズが毎回変わるのを避けたいなら tight は外すのが安定
    fig.savefig(out_fig)

    print(f"Saved figure: {out_fig}")
    print(f"Saved CSV: {folder / args.out_csv}")

    expected = len(ROD_COUNTS) * len(CONE_LEVELS)
    got = np.isfinite(Z).sum()
    if got < expected:
        print(f"Note: computed {got}/{expected} cells. Missing files will be NaN in the CSV/heatmap.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
