#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
設定だけで cell / synapse モード切り替え
刺激区間の強調表示
---------------------------------------------------------------------------
実行方法:
  python plot_pdf_stim_modes.py

このプログラムがすること:
  - フォルダ内から、指定パターンに一致するtxtファイルを検索する
  - 先頭2列を (time(ms), value) として読み込む
  - トレース全体をプロットする
  - 刺激区間だけ別の色で重ね描きして強調する
  - 必要に応じて、刺激バー（黒いバー）を描画する
  - 各ファイルと同じ場所に、同じベース名の PDF を保存する

入力ファイル形式:
  少なくとも2列（time(ms), value）があること
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# ============================ CONFIG (EDIT HERE) ============================
# ---- Mode ----
# "cell"   : membrane potential traces (ylabel + fixed ylim)
# "syn"    : synapse/other variable traces (empty ylabel + auto ylim by default)
MODE: str = "cell"  # "cell" or "syn"

# ---- Target (used to build default folder/pattern and filename sorting) ----
# Examples: "OFF_GC", "ON_GC", "AIIAC", "ONCB2ONGC_u"
TARGET: str = "OFF_GC"

# ---- Folder / pattern ----
# If ROOT_DIR is empty, defaults are used:
#   cell:  "{TARGET}_fovea"
#   syn :  "{TARGET}"
ROOT_DIR: str = ""

# If PATTERN is empty, default is: "{TARGET}_R*_C*.txt"
PATTERN: str = ""

# Search subfolders?
RECURSIVE: bool = False

# ---- Plot window ----
XLIM_MS: Tuple[float, float] = (1000.0, 6000.0)

# ---- Labels / Y-limits by mode ----
CELL_YLABEL: str = "Membrane potential (mV)"
CELL_YLIM: Tuple[float, float] = (-100.0, 20.0)

SYN_YLABEL: str = ""                   # often blank for synapse variable plots
SYN_YLIM: Optional[Tuple[float, float]] = None  # None = autoscale

# ---- Stimulus ----
STIM_START_MS: float = 2000.0
STIM_END_MS: float = 3000.0
DRAW_STIM_BAR: bool = True

# ---- Style ----
BASE_COLOR: str = "#274A78"   # blue-ish
STIM_COLOR: str = "#FF6F40"   # orange
LINEWIDTH: float = 1.5

STIM_BAR_COLOR: str = "black"
STIM_BAR_HEIGHT_FRAC: float = 0.015
STIM_BAR_Y_FRAC: float = 0.02
# ==========================================================================


def build_filename_re(target: str) -> re.Pattern:
    # Accept: <TARGET>_R###_C##.txt (case-insensitive R/C)
    return re.compile(rf"^{re.escape(target)}_[Rr](\d+)_[Cc](\d+)\.txt$", re.IGNORECASE)


def parse_rc_from_name(path: Path, filename_re: re.Pattern) -> Tuple[int, int]:
    m = filename_re.match(path.name)
    if not m:
        return (10**9, 10**9)
    return int(m.group(1)), int(m.group(2))


def iter_files(folder: Path, pattern: str, recursive: bool, filename_re: re.Pattern) -> List[Path]:
    it = folder.rglob(pattern) if recursive else folder.glob(pattern)
    files = [p for p in it if p.is_file()]
    files.sort(key=lambda p: parse_rc_from_name(p, filename_re))
    return files


def read_trace(txt_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read CSV-like text file, use first two columns as (t, y)."""
    try:
        df = pd.read_csv(txt_path, comment="#")
    except Exception as e:
        raise RuntimeError(f"Failed to read {txt_path}: {e}")

    if df.shape[1] < 2:
        raise ValueError(f"{txt_path} does not have at least two columns")

    t = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float, copy=False)
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float, copy=False)

    if np.isnan(t).any() or np.isnan(y).any():
        raise ValueError(f"{txt_path} contains non-numeric entries in the first two columns")

    return t, y


def plot_and_save_pdf(
    txt_path: Path,
    xlim: Tuple[float, float],
    ylim: Optional[Tuple[float, float]],
    stim: Tuple[float, float],
    xlabel: str,
    ylabel: str,
    draw_stim_bar: bool,
) -> Path:
    t, y = read_trace(txt_path)
    stim_start, stim_end = float(stim[0]), float(stim[1])

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)

    # Base trace
    ax.plot(t, y, color=BASE_COLOR, lw=LINEWIDTH, zorder=1)

    # Stimulus-highlight overlay
    if np.isfinite(stim_start) and np.isfinite(stim_end) and stim_end > stim_start:
        mask = (t >= stim_start) & (t <= stim_end)
        if np.any(mask):
            ax.plot(t[mask], y[mask], color=STIM_COLOR, lw=LINEWIDTH, zorder=2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    # Stimulus bar
    if draw_stim_bar and (stim_end > stim_start):
        ymin, ymax = ax.get_ylim()
        yr = (ymax - ymin) if (ymax > ymin) else 1.0
        rect_h = yr * STIM_BAR_HEIGHT_FRAC
        y_pos = ymin + yr * STIM_BAR_Y_FRAC
        ax.add_patch(
            patches.Rectangle(
                (stim_start, y_pos),
                width=(stim_end - stim_start),
                height=rect_h,
                color=STIM_BAR_COLOR,
                zorder=10,
            )
        )

    out_pdf = txt_path.with_suffix(".pdf")
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def main() -> int:
    mode = MODE.strip().lower()
    if mode not in ("cell", "syn"):
        print(f"Error: MODE must be 'cell' or 'syn' (got {MODE!r})", file=sys.stderr)
        return 2

    target = TARGET.strip()
    if not target:
        print("Error: TARGET is empty", file=sys.stderr)
        return 2
       
    default_root = f"{target}_fovea" if mode == "cell" else f"{target}"
    default_pattern = f"{target}_R*_C*.txt"

    folder = Path(ROOT_DIR) if ROOT_DIR.strip() else Path(default_root)
    pattern = PATTERN if PATTERN.strip() else default_pattern

    if not folder.exists() or not folder.is_dir():
        print(f"Error: folder is not a directory: {folder}", file=sys.stderr)
        return 2

    ylabel = CELL_YLABEL if mode == "cell" else SYN_YLABEL
    ylim = CELL_YLIM if mode == "cell" else SYN_YLIM
    stim = (STIM_START_MS, STIM_END_MS)

    filename_re = build_filename_re(target)
    files = iter_files(folder, pattern, RECURSIVE, filename_re)

    print(f"MODE={mode}  TARGET={target}")
    print(f"FOLDER={folder}  PATTERN={pattern!r}  RECURSIVE={RECURSIVE}")
    print(f"XLIM_MS={XLIM_MS}  STIM={stim}  YLIM={ylim}  DRAW_STIM_BAR={DRAW_STIM_BAR}")
    print(f"Files found: {len(files)}")

    if not files:
        print("No files matched. Check FOLDER/PATTERN/TARGET in CONFIG.")
        return 0

    n_ok = 0
    for p in files:
        try:
            out_pdf = plot_and_save_pdf(
                txt_path=p,
                xlim=XLIM_MS,
                ylim=ylim,
                stim=stim,
                xlabel="Time (ms)",
                ylabel=ylabel,
                draw_stim_bar=DRAW_STIM_BAR,
            )
            print(f"[OK] {p.name} -> {out_pdf.name}")
            n_ok += 1
        except Exception as e:
            print(f"[ERR] {p.name}: {e}", file=sys.stderr)

    print(f"Done. {n_ok}/{len(files)} PDFs created.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
