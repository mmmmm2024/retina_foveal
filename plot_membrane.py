#!/usr/bin/env python3
"""
フォルダ内の波形txt(CSV形式・2列)を一括で読み込み、
時間(ms)–膜電位(mV)の折れ線グラフを作って、各ファイルと同名のPDFを保存
（xlim/ylimで表示範囲を固定、必要ならサブフォルダも探索）
"""


from __future__ import annotations
import argparse
from pathlib import Path
import re
import sys
from typing import Iterable, Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt

TARGET = "AIIAC"
# TARGET = "ON_GC"
# === USER CONFIG (edit here once) ============================================
# If ROOT_DIR is empty (""), the script uses the folder where THIS FILE sits.
# If you prefer a fixed location, set ROOT_DIR to that folder path.
# ROOT_DIR: str = f"Dim_{TARGET}_RB-AII75%_AII-ONCB0%_AII-AII100%_ONCB-ONCB100%"          # e.g., "/Users/you/data/AIIAC" or "C:/data/AIIAC"; empty → script folder
# ROOT_DIR: str = f"{TARGET}_RB-AII50%_AII-ONCB0%_AII-AII100%_ONCB-ONCB100%"
# ROOT_DIR: str = "Dim_AIIAC_fovea"
ROOT_DIR: str = "Dim_AIIAC_8-9M_gRBC2ACx_gjAC2CBy"
RECURSIVE: bool = False      # True to also search subfolders
# PATTERN: str = f"{TARGET}_R*_C*.txt"
PATTERN: str = f"{TARGET}_x*_y*.txt"
# PATTERN: str = "RBC_R*_C*.txt"
# PATTERN: str = "ONCB_R*_C*.txt"
# PATTERN: str = f"{TARGET}_R*_C*.txt"
XLIM: Optional[Tuple[float, float]] = (1000, 1300)   # e.g., (1000, 6000) or None
YLIM: Optional[Tuple[float, float]] =(-65, -40) # e.g., (-70, -40) or None
# ============================================================================

FILENAME_RE = re.compile(rf"^{TARGET}_R(\d+)_C(\d+)\.txt$", re.IGNORECASE)
# FILENAME_RE = re.compile(r"^RBC_R(\d+)_C(\d+)\.txt$", re.IGNORECASE)
# FILENAME_RE = re.compile(r"^ONCB_R(\d+)_C(\d+)\.txt$", re.IGNORECASE)

def parse_rc_from_name(path: Path) -> Tuple[int, int]:
    """Return (R, C) numbers from a filename like AIIAC_R001_C00.txt.
    If it doesn't match, return a large tuple so sorting pushes it to the end.
    """
    m = FILENAME_RE.match(path.name)
    if not m:
        return (10**9, 10**9)
    r = int(m.group(1))
    c = int(m.group(2))
    return (r, c)


def iter_files(folder: Path, pattern: str, recursive: bool) -> Iterable[Path]:
    if recursive:
        files = sorted(folder.rglob(pattern), key=parse_rc_from_name)
    else:
        files = sorted(folder.glob(pattern), key=parse_rc_from_name)
    return [p for p in files if p.is_file()]


def read_trace(txt_path: Path) -> Tuple[pd.Series, pd.Series]:
    """Read a two-column CSV with a header. Return (time_ms, voltage_mV)."""
    try:
        df = pd.read_csv(txt_path, comment="#")
    except Exception as e:
        raise RuntimeError(f"Failed to read {txt_path}: {e}")

    if df.shape[1] < 2:
        raise ValueError(f"{txt_path} does not have at least two columns")

    # Be robust to exact column names; assume first two columns are time and voltage
    time_col = df.columns[0]
    volt_col = df.columns[1]
    return df[time_col], df[volt_col]


def plot_and_save_pdf(txt_path: Path, out_pdf: Optional[Path] = None,
                      xlim: Optional[Tuple[float, float]] = None,
                      ylim: Optional[Tuple[float, float]] = None) -> Path:
    """Create a line plot and save to PDF next to the source file.
    Returns the PDF path.
    """
    t_ms, v_mV = read_trace(txt_path)

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.plot(t_ms.values, v_mV.values)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Membrane potential (mV)")

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    if out_pdf is None:
        out_pdf = txt_path.with_suffix(".pdf")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def resolve_default_folder(cli_folder: Optional[Path]) -> Path:
    """Pick the folder to process from CLI or CONFIG, with sensible fallbacks."""
    if cli_folder is not None:
        return cli_folder
    if ROOT_DIR.strip():
        return Path(ROOT_DIR)
    # Fallback: script's directory (or CWD if __file__ not defined)
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()


def main(argv: Optional[Iterable[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Batch-plot AIIAC traces to PDFs.")
    p.add_argument("folder", nargs="?", type=Path, default=None,
                   help="Folder with AIIAC_R*_C*.txt (default: CONFIG ROOT_DIR or script folder)")
    p.add_argument("--pattern", default=PATTERN,
                   help=f"Glob pattern to match files (default: {PATTERN!r})")
    # Note: nargs=2 cannot have a tuple default. We'll merge CLI with CONFIG below.
    p.add_argument("--xlim", nargs=2, type=float, metavar=("XMIN", "XMAX"), default=None,
                   help="Optional x-axis limits in ms, e.g., --xlim 1000 6000")
    p.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"), default=None,
                   help="Optional y-axis limits in mV, e.g., --ylim -70 -40")
    p.add_argument("--recursive", action="store_true", default=RECURSIVE,
                   help=f"Search subfolders recursively (default: {RECURSIVE})")

    args = p.parse_args(argv)

    folder: Path = resolve_default_folder(args.folder)
    if not folder.exists() or not folder.is_dir():
        print(f"Error: {folder} is not a directory", file=sys.stderr)
        return 2

    # Merge CLI with CONFIG
    pattern = args.pattern or PATTERN
    xlim = tuple(args.xlim) if args.xlim else (tuple(XLIM) if XLIM else None)
    ylim = tuple(args.ylim) if args.ylim else (tuple(YLIM) if YLIM else None)
    recursive = bool(args.recursive)

    files = iter_files(folder, pattern, recursive)
    if not files:
        print(f"No files matched pattern '{pattern}' in {folder}")
        return 0

    print(f"Found {len(files)} files in {folder}. Writing PDFs next to their sources...")
    n_ok = 0
    for f in files:
        try:
            out_pdf = plot_and_save_pdf(f, xlim=xlim, ylim=ylim)
            print(f"[OK] {f.name} -> {out_pdf.name}")
            n_ok += 1
        except Exception as e:
            print(f"[ERR] {f.name}: {e}", file=sys.stderr)

    print(f"Done. {n_ok}/{len(files)} PDFs created.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
