#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
二次変性時のAIIAC膜電位波形(AIIAC_x*_y*.txt)から 5–15 Hz のバンドパワーを計算
gRBC2AC(%) × gjAC2CB(%) のスイープ結果を2DヒートマップとCSVで出力

入力: ROOT_DIR内の AIIAC_x{gRBC2AC}_y{gjAC2CB}.txt
処理: 1000–6000 ms を抽出 → FFT/PSD → 5–15 Hz を積分
出力: AIIAC_bandpower_gRBC2AC_vs_gjAC2CB.csv / .pdf
"""

from __future__ import annotations
from pathlib import Path
import re
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======================= CONFIG ============================
ROOT_DIR = "Dim_AIIAC_8-9M_gRBC2ACx_gjAC2CBy"

CROP_MS = (1000.0, 6000.0)
BAND_HZ = (5.0, 15.0)
USE_HANN = True
DETREND_LINEAR = False

GRBC2AC_LEVELS = list(range(0, 101, 5))   # Y-axis
GJAC2CB_LEVELS = list(range(0, 101, 5))   # X-axis

CB_CMAP = "Blues"
CB_VMIN = 0.0
CB_VMAX = 25

OUT_CSV = "AIIAC_bandpower_gRBC2AC_vs_gjAC2CB.csv"
OUT_PNG = "AIIAC_bandpower_gRBC2AC_vs_gjAC2CB.pdf"
# ===========================================================

FILENAME_RE = re.compile(r"^AIIAC_x(\d+)_y(\d+)\.txt$", re.IGNORECASE)


def load_aiiac(path: Path):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1]


def periodogram_psd(x, fs, use_hann=True):
    N = len(x)
    x = x - np.mean(x)
    w = np.hanning(N) if use_hann else np.ones(N)
    xw = x * w
    U = (w**2).sum() / N

    X = np.fft.rfft(xw)
    Pxx = (np.abs(X)**2) / (fs * N * U)
    Pxx[1:-1] *= 2
    f = np.fft.rfftfreq(N, d=1/fs)
    return f, Pxx


def bandpower(path: Path) -> float:
    t_ms, v = load_aiiac(path)
    mask = (t_ms >= CROP_MS[0]) & (t_ms <= CROP_MS[1])
    t = t_ms[mask]
    v = v[mask]

    fs = 1000.0 / np.median(np.diff(t))
    f, Pxx = periodogram_psd(v, fs, USE_HANN)

    idx = (f >= BAND_HZ[0]) & (f <= BAND_HZ[1])
    return float(np.trapz(Pxx[idx], f[idx]))


def main():
    folder = Path(ROOT_DIR)
    Z = np.full(
        (len(GRBC2AC_LEVELS), len(GJAC2CB_LEVELS)),
        np.nan
    )

    for p in folder.glob("AIIAC_x*_y*.txt"):
        m = FILENAME_RE.match(p.name)
        if not m:
            continue

        x = int(m.group(1))
        y = int(m.group(2))

        if x not in GRBC2AC_LEVELS or y not in GJAC2CB_LEVELS:
            continue

        i = GRBC2AC_LEVELS.index(x)
        j = GJAC2CB_LEVELS.index(y)

        try:
            Z[i, j] = bandpower(p)
        except Exception as e:
            print(f"[SKIP] {p.name}: {e}")

    # CSV
    df = pd.DataFrame(
        Z,
        index=[f"{v}" for v in GRBC2AC_LEVELS],
        columns=[f"{v}" for v in GJAC2CB_LEVELS]
    )
    df.index.name = "gRBC2AC_%"
    df.columns.name = "gjAC2CB_%"
    df.to_csv(folder / OUT_CSV)

    # Plot
    fig, ax = plt.subplots(figsize=(6.8, 6.2), dpi=150)
    im = ax.imshow(
        np.ma.masked_invalid(Z),
        origin="lower",
        extent=(-2.5, 102.5, -2.5, 102.5),
        cmap=CB_CMAP,
        vmin=CB_VMIN,
        vmax=CB_VMAX
    )

    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Band power 5–15 Hz (mV$^2$)")

    ax.set_xticks(GJAC2CB_LEVELS)
    ax.set_yticks(GRBC2AC_LEVELS)

    ax.set_xlabel("gap conductance between AIIAC and CB(%)")
    ax.set_ylabel("conductance from RBC to AIIAC(%)")
    # ax.set_title("AIIAC oscillation power (5–15 Hz)")
    ax.set_aspect("equal")

    fig.savefig(folder / OUT_PNG, bbox_inches="tight")
    print(f"Saved PNG: {OUT_PNG}")
    print(f"Saved CSV: {OUT_CSV}")


if __name__ == "__main__":
    main()
