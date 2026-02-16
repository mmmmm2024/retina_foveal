#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 mGC発火率ヒートマップ + knee検出（Cones×Rods）
 ON/OFF mGCの膜電位（{TARGET}_R*_C*.txt）からスパイクを検出
 発火率(Hz)の行列を作ってヒートマップ化、さらに急変点(knee)を推定して線で重ね描き

 入力:
 - {TARGET}_fovea フォルダ内の CSV 形式txt（先頭2列が time(ms), voltage(mV)）
   例: ON_GC_R400_C20.txt など（Rod=R, Cone=C）

 出力:
 - ヒートマップPDF（knee線つき）
 - 発火率行列CSV
 - knee推定結果（Rodごと / Coneごと）
 - 変化量（drop）指標CSV（絶対/相対）

処理:
 1) 閾値 + ヒステリシス閾値でスパイク数をカウント
 2) 解析窓(start_time–end_time)で発火率(Hz)を計算
 3) Rod×Cone の行列Zを作成して保存
 4) 行/列方向の最大の変化から knee を推定（任意でPAVAで単調化）
"""

from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable  # ★追加

# モードと機能フラグ
MODE = "ON"    # "ON" or "OFF"
ENABLE_STEEPEST = False
ENABLE_RAPID3D  = False

# ★PAVA を使うかどうか（True: PAVAあり / False: PAVAなし）
USE_PAVA = False

# 入力
TARGET = "ON_GC" if MODE == "ON" else "OFF_GC"
ROOT_DIR = Path(f"{TARGET}_fovea")
PATTERN = f"{TARGET}_R*_C*.txt"
RECURSIVE = False

# スパイク検出
thr_hi = 0.0
thr_lo = -20.0
start_time = 1000.0
end_time   = 6000.0

# 表示とアルゴリズム設定
MOVING_AVG_WIN = 1

# ★正方形PDF
FIGSIZE = (7.5, 7.5)

DPI = 300
VMIN = 0.0
VMAX = 7.0 if MODE == "ON" else 17.0

CMAP = LinearSegmentedColormap.from_list("white_to_blue", ["#FFFFFF", "#1f4e97"])

ROD_MAX  = 400.0
CONE_MAX = 20.0

cone_vals = [20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0]
rod_vals_display = [400, 360, 320, 280, 240, 200, 160, 120, 80, 40, 0]

# 出力ファイル
PDF_PATH        = ROOT_DIR / f"{MODE}mGC_ConesxRods_robust_knee.pdf"
CSV_MAT         = ROOT_DIR / f"{MODE}mGC_heatmap_matrix.csv"
CSV_KNEE_ROW    = ROOT_DIR / f"{MODE}mGC_knee_by_rod.csv"
CSV_KNEE_COL    = ROOT_DIR / f"{MODE}mGC_knee_by_cone.csv"
CSV_REL_ROD     = ROOT_DIR / f"{MODE}mGC_rel_drop_by_rod.csv"
CSV_REL_CONE    = ROOT_DIR / f"{MODE}mGC_rel_drop_by_cone.csv"

RAPID_SMOOTH_WIN = 1
RAPID_PERCENTILE = 90.0
CSV_GRADMAG     = ROOT_DIR / f"{MODE}mGC_gradient_mag.csv"
CSV_RAPID_MASK  = ROOT_DIR / f"{MODE}mGC_rapid_region_mask.csv"

FILENAME_RE = re.compile(rf"^{re.escape(TARGET)}_R(\d+)_C(\d+)\.txt$", re.IGNORECASE)

# 補助関数
def iter_files(folder: Path) -> list[Path]:
    it = folder.rglob(PATTERN) if RECURSIVE else folder.glob(PATTERN)
    return sorted([p for p in it if p.is_file()], key=lambda p: p.name)

def rc_from_name(p: Path) -> tuple[int, int] | None:
    m = FILENAME_RE.match(p.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))

def load_trace(path: Path) -> tuple[np.ndarray, np.ndarray]:
    try:
        data = np.loadtxt(path, delimiter=",", skiprows=1)
    except Exception:
        data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError(f"{path.name}: not 2-column numeric data")
    return data[:, 0], data[:, 1]

def firing_rate_hz(t_ms: np.ndarray, v_mV: np.ndarray) -> float:
    mask = (t_ms >= start_time) & (t_ms <= end_time)
    if not np.any(mask):
        return np.nan
    v = v_mV[mask]

    spike_count = 0
    armed = True
    for i in range(1, len(v)):
        if armed:
            if v[i-1] <= thr_hi and v[i] > thr_hi:
                spike_count += 1
                armed = False
        else:
            if v[i] < thr_lo:
                armed = True

    duration_s = (end_time - start_time) / 1000.0
    return spike_count / duration_s if duration_s > 0 else np.nan

def moving_average_1d(y, win=3):
    y = np.asarray(y, float)
    if win is None or win <= 1 or len(y) < 2:
        return y.copy()
    win = int(win)
    if win % 2 == 0:
        win += 1
    pad = win // 2
    ypad = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(win) / win
    return np.convolve(ypad, kernel, mode="valid")

def moving_average_2d(A, win=3):
    if win is None or win <= 1:
        return np.array(A, float, copy=True)
    B = np.apply_along_axis(moving_average_1d, 1, np.asarray(A, float), win)
    C = np.apply_along_axis(moving_average_1d, 0, B, win)
    return C

def pava_increasing(y, w=None):
    y = list(map(float, y))
    if w is None:
        w = [1.0] * len(y)
    else:
        w = list(map(float, w))
    blocks = [{"sw": w[i], "sy": w[i]*y[i], "len": 1} for i in range(len(y))]
    i = 0
    while i < len(blocks) - 1:
        a = blocks[i]["sy"] / blocks[i]["sw"]
        b = blocks[i+1]["sy"] / blocks[i+1]["sw"]
        if a > b:
            blocks[i]["sw"] += blocks[i+1]["sw"]
            blocks[i]["sy"] += blocks[i+1]["sy"]
            blocks[i]["len"] += blocks[i+1]["len"]
            del blocks[i+1]
            if i > 0:
                i -= 1
        else:
            i += 1
    out = []
    for b in blocks:
        out.extend([b["sy"]/b["sw"]] * b["len"])
    return np.array(out, float)

def pava_decreasing(y, w=None):
    return -pava_increasing(-np.asarray(y, float), w=w)

def knee_by_max_gradient(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2:
        return np.nan
    dx = np.diff(x); dy = np.diff(y)
    with np.errstate(divide="ignore", invalid="ignore"):
        slope = np.abs(dy / dx)
    j = int(np.nanargmax(slope))
    return 0.5 * (x[j] + x[j+1])

def max_drop_metrics(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2:
        return np.nan, np.nan, np.nan, np.nan
    drops = -np.diff(y)
    j = int(np.nanargmax(drops))
    drop_abs = float(drops[j])
    dyn = float(np.nanmax(y) - np.nanmin(y))
    rel_range = drop_abs / dyn if dyn > 0 else np.nan
    rel_level = drop_abs / y[j] if y[j] > 0 else np.nan
    knee_pos = 0.5 * (x[j] + x[j+1])
    return knee_pos, drop_abs, rel_range, rel_level

def gradient_magnitude_percent(S, rod_vals, cone_vals):
    S = np.asarray(S, float)
    nr, nc = S.shape
    rod_pct  = np.array(rod_vals,  float) / ROD_MAX  * 100.0
    cone_pct = np.array(cone_vals, float) / CONE_MAX * 100.0

    Gx = np.zeros_like(S, float)
    Gy = np.zeros_like(S, float)

    for j in range(nc):
        if j == 0:
            dx = cone_pct[j] - cone_pct[j+1]
            Gx[:, j] = (S[:, j] - S[:, j+1]) / dx
        elif j == nc-1:
            dx = cone_pct[j-1] - cone_pct[j]
            Gx[:, j] = (S[:, j-1] - S[:, j]) / dx
        else:
            dx = cone_pct[j-1] - cone_pct[j+1]
            Gx[:, j] = (S[:, j-1] - S[:, j+1]) / dx

    for i in range(nr):
        if i == 0:
            dy = rod_pct[i] - rod_pct[i+1]
            Gy[i, :] = (S[i, :] - S[i+1, :]) / dy
        elif i == nr-1:
            dy = rod_pct[i-1] - rod_pct[i]
            Gy[i, :] = (S[i-1, :] - S[i, :]) / dy
        else:
            dy = rod_pct[i-1] - rod_pct[i+1]
            Gy[i, :] = (S[i-1, :] - S[i+1, :]) / dy

    return np.sqrt(Gx**2 + Gy**2)

# メイン処理
def main():
    if not ROOT_DIR.exists():
        raise SystemExit(f"Folder not found: {ROOT_DIR}")

    files = iter_files(ROOT_DIR)
    if not files:
        raise SystemExit(f"No files matched: {PATTERN} in {ROOT_DIR}")

    rates: dict[tuple[int, int], float] = {}
    for p in files:
        rc = rc_from_name(p)
        if rc is None:
            continue
        t, v = load_trace(p)
        rates[rc] = firing_rate_hz(t, v)

    Z = np.full((len(rod_vals_display), len(cone_vals)), np.nan, float)
    for i, r in enumerate(rod_vals_display):
        rr = 1 if r == 0 else r
        for j, c in enumerate(cone_vals):
            Z[i, j] = rates.get((rr, c), np.nan)

    pd.DataFrame(Z, index=rod_vals_display, columns=cone_vals).to_csv(CSV_MAT, index_label="Rod")

    # knees: row-wise
    cols = np.array(cone_vals, float)
    x_line = np.full(len(rod_vals_display), np.nan, float)
    cone_knees = np.full(len(rod_vals_display), np.nan, float)

    rod_drop_abs  = np.full(len(rod_vals_display), np.nan, float)
    rod_rel_range = np.full(len(rod_vals_display), np.nan, float)
    rod_rel_level = np.full(len(rod_vals_display), np.nan, float)

    for i in range(Z.shape[0]):
        y_row = Z[i, :]
        mask = ~np.isnan(y_row)
        if mask.sum() < 2:
            continue
        x_cones = cols[mask]
        y_vals  = y_row[mask]

        y_ma  = moving_average_1d(y_vals, win=MOVING_AVG_WIN)
        y_use = pava_decreasing(y_ma) if USE_PAVA else y_ma

        x_star = knee_by_max_gradient(x_cones, y_use)
        cone_knees[i] = x_star

        _, dabs, rrange, rlevel = max_drop_metrics(x_cones, y_use)
        rod_drop_abs[i]  = dabs
        rod_rel_range[i] = rrange
        rod_rel_level[i] = rlevel

        if np.isnan(x_star):
            continue
        if x_star >= cols[0]:
            x_line[i] = 0.0
        elif x_star <= cols[-1]:
            x_line[i] = len(cols) - 1.0
        else:
            k = np.where((cols[:-1] >= x_star) & (x_star >= cols[1:]))[0]
            x_line[i] = (k[0] + 0.5) if len(k) else float(np.argmin(np.abs(cols - x_star)))

    pd.DataFrame({"Rod": rod_vals_display, "knee_cone_value": cone_knees}).to_csv(CSV_KNEE_ROW, index=False)

    # knees: col-wise
    rows = np.array(rod_vals_display, float)
    rod_knees_per_col = np.full(len(cone_vals), np.nan, float)

    col_drop_abs  = np.full(len(cone_vals), np.nan, float)
    col_rel_range = np.full(len(cone_vals), np.nan, float)
    col_rel_level = np.full(len(cone_vals), np.nan, float)

    for j in range(Z.shape[1]):
        y_col = Z[:, j]
        mask = ~np.isnan(y_col)
        if mask.sum() < 2:
            continue
        x_rods = rows[mask]
        y_vals = y_col[mask]

        y_ma  = moving_average_1d(y_vals, win=MOVING_AVG_WIN)
        y_use = pava_decreasing(y_ma) if USE_PAVA else y_ma

        kpos, dabs, rrange, rlevel = max_drop_metrics(x_rods, y_use)
        rod_knees_per_col[j] = kpos
        col_drop_abs[j]  = dabs
        col_rel_range[j] = rrange
        col_rel_level[j] = rlevel

    pd.DataFrame({"Cone": cone_vals, "knee_rod_value": rod_knees_per_col}).to_csv(CSV_KNEE_COL, index=False)

    pd.DataFrame({
        "Rod": rod_vals_display,
        "knee_cone_value": cone_knees,
        "abs_drop_Hz": rod_drop_abs,
        "rel_drop_range": rod_rel_range,
        "rel_drop_level": rod_rel_level,
    }).to_csv(CSV_REL_ROD, index=False)

    pd.DataFrame({
        "Cone": cone_vals,
        "knee_rod_value": rod_knees_per_col,
        "abs_drop_Hz": col_drop_abs,
        "rel_drop_range": col_rel_range,
        "rel_drop_level": col_rel_level,
    }).to_csv(CSV_REL_CONE, index=False)

    # Optional rapid 
    if ENABLE_RAPID3D:
        Z_smooth = moving_average_2d(Z, win=RAPID_SMOOTH_WIN)
        Gmag = gradient_magnitude_percent(Z_smooth, rod_vals_display, cone_vals)
        thr = np.nanpercentile(Gmag, RAPID_PERCENTILE)
        rapid_mask = (Gmag >= thr)
        pd.DataFrame(Gmag, index=rod_vals_display, columns=cone_vals).to_csv(CSV_GRADMAG, index_label="Rod")
        pd.DataFrame(rapid_mask.astype(int), index=rod_vals_display, columns=cone_vals).to_csv(CSV_RAPID_MASK, index_label="Rod")

    # Plot (square, no overlap, aligned colorbar)
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)

    gs = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[1.0, 0.18],
        hspace=0.06,
        left=0.10, right=0.96,
        top=0.93, bottom=0.07
    )

    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[1, 0])
    ax_leg.axis("off")

    Z_masked = np.ma.masked_invalid(Z)
    norm = Normalize(vmin=VMIN, vmax=VMAX)

    im = ax.imshow(
        Z_masked,
        origin="upper",
        interpolation="nearest",
        cmap=CMAP,
        norm=norm
    )

    ax.set_aspect("equal", adjustable="box")

    xticks = np.arange(len(cone_vals))
    xlabels_pct = [f"{int(round(100*c/CONE_MAX))}" for c in cone_vals]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels_pct)

    yticks = np.arange(len(rod_vals_display))
    ylabels_pct = [f"{int(round(100*r/ROD_MAX))}" for r in rod_vals_display]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels_pct)

    ax.set_xlabel("Cone Survival rate(%)", labelpad=1)
    ax.set_ylabel("Rod Survival rate(%)")
    ax.set_title(f"{MODE} mGC firing rate [Hz]", pad=6)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4.2%", pad=0.06)
    cbar = fig.colorbar(im, cax=cax, label="Firing rate [Hz]")
    cbar.set_ticks(np.arange(int(np.floor(VMIN)), int(np.ceil(VMAX)) + 1, 1))
    cbar.ax.tick_params(pad=1)

    # knee lines
    y_idx = np.arange(len(rod_vals_display), dtype=float)

    valid_row = ~np.isnan(x_line)
    h_row, = ax.plot(
        x_line[valid_row], y_idx[valid_row],
        color="red", linewidth=2.0, label="knee along Cone"
    )

    idxs = np.arange(len(rod_vals_display), dtype=float)
    rod_asc, idx_asc = rows[::-1], idxs[::-1]
    y_line_col = np.array([
        np.interp(rv, rod_asc, idx_asc) if np.isfinite(rv) else np.nan
        for rv in rod_knees_per_col
    ], dtype=float)

    x_idx = np.arange(len(cone_vals), dtype=float)
    valid_col = ~np.isnan(y_line_col)
    h_col, = ax.plot(
        x_idx[valid_col], y_line_col[valid_col],
        linestyle="-", color="#f5b000", linewidth=2.0, label="knee along Rod"
    )

    handles = [h_row, h_col]
    labels  = [h_row.get_label(), h_col.get_label()]
    ax_leg.legend(
        handles, labels,
        loc="center",
        frameon=True,
        ncol=1,
        fontsize=9,
        borderpad=0.5,
        labelspacing=0.35,
        handlelength=2.8
    )

    fig.savefig(PDF_PATH)
    plt.close(fig)

    print("Saved:")
    for p in [PDF_PATH, CSV_MAT, CSV_KNEE_ROW, CSV_KNEE_COL, CSV_REL_ROD, CSV_REL_CONE]:
        print(" ", p)

if __name__ == "__main__":
    main()
