"""
密度モデルから細胞を2D平面上にランダム配置
KDTreeで各Ganglion(Ganglion_Lee)の受容野内に入る細胞数をカウント
偏心度(deg)ビンごとに平均/SDとGC受容野半径・面積の統計をExcelに出力

手順:
- diff_exp の密度モデル → 10–15mmの同心円リングで点群生成
- 各細胞タイプの点群にKDTreeを作成
- 各Ganglion点について、半径モデル(midget/parasol)のRF内の近傍点数をカウント
- Ecc_degを丸めてビン化し、平均・SD・n_GCを集計してシート保存
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from collections import OrderedDict

np.random.seed(0)
ARC_W   = 0.1
R_MIN   = 10.0
R_MAX   = 15.0
RADII   = np.arange(R_MIN, R_MAX, ARC_W)
MM_PER_DEG = 0.29
mm2deg = lambda r_mm: r_mm / MM_PER_DEG


# GC 樹状突起半径モデル（μm → mm）
def r_midget(r_mm):  return (8.64 * r_mm**1.04) / 2 / 1000
def r_parasol(r_mm): return (70.2 * r_mm**0.60) / 2 / 1000

# 密度モデル
def diff_exp(x, c1, k1, c2, k2, c3=0, k3=0):
    return c1*np.exp(k1*x) + c2*np.exp(k2*x) + c3*np.exp(k3*x)

PARAMS = OrderedDict([
    # ("Cone",      (3.673e5, -7.828e0,   -2.000e5, -2.000e3,   2.034e4,  -2.164e-1)),
    # ("S-Cone",    (5.048e3, -3.014e0,   -1.100e4, -7.869e0,   1.455e3,  -1.370e-1)),
    # ("Rod",       (5.855e5, -1.388e-1,  -5.989e5, -2.998e-1)),
    # ("Ganglion",  (5.717e5, -1.031e0,   -6.000e5, -1.261e0,   -5.527e1, -9.658e1)),
    ("Ganglion_Lee",  (8.475e5, -1.258e-1,   -8.691e5, -1.556e0)),
    # ("OFFMBC",    (4.208e5, -1.225e0,   -4.551e5, -1.387e0,   1.230e4,  -6.986e-2)),
    # ("ONBC",      (1.788e5, -4.755e-1,  -1.836e5, -5.618e-1,   8.894e3,  -1.832e-2)),
    # ("DB3a",      (4.275e3, -1.949e-1,  -9.614e3, -5.831e0,   3.336e2,   5.811e-3)),
    # ("DB3b",      (7.061e3, -1.511e-1,  -7.869e3, -3.120e0,  -3.653e3,  -3.347e-1)),
    # ("Horizontal",(2.051e4, -3.850e-1,  -4.246e4, -4.196e0,  2.672e3,  -1.097e-2)),
    # ("H1",        (1.659e4, -3.344e-1,  -4.227e6, -1.633e1,   1.280e3,   4.952e-2)),
    # ("H2",        (5.131e3, -2.197e-1,  -5.440e3, -6.706e-1,  -8.006e-1, -5.171e1)),
    # ("GlyAC",     (1.614e4, -1.260e-1,  -3.944e5, -5.923e0,  -6.498e2, -5.675e1)),
    ("AIIAC",    (1.297e5, -1.114e0,  -1.392e5, -1.217e0,    2.553e3,  -7.214e-2)),
    # ("Müller",    (2.011e5, -1.845e0,  -2.369e5, -2.362e0,    1.379e4,  -2.815e-2)),
])

density_fn = {n: (lambda p: (lambda r: diff_exp(r, *p)))(v) for n, v in PARAMS.items()}

# ランダム点生成
def sample_points(radii, dens_fn):
    pts, rs = [], []
    for r in radii:
        n = int(dens_fn(r) * (2 * np.pi * r * ARC_W))
        if n == 0:
            continue
        ang   = np.random.rand(n) * 2 * np.pi
        r_jit = np.random.uniform(r - ARC_W/2, r + ARC_W/2, n)
        pts.append(np.c_[r_jit*np.cos(ang), r_jit*np.sin(ang)])
        rs.extend(r_jit)
    if not pts:
        return np.empty((0, 2)), np.array([])
    return np.vstack(pts), np.asarray(rs)

# 全細胞サンプリング & KD-Tree
pts_dict  = {}
radius_gc = None
for cell, fn in density_fn.items():
    pts, r_arr = sample_points(RADII, fn)
    pts_dict[cell] = {
        "pts": pts,
        "tree": cKDTree(pts) if pts.size else None
    }
    if cell == "Ganglion_Lee":
        radius_gc = r_arr

# 統計計算 → Excel
def build_stats(r_fn, sheet, writer):
    cols = ["Ecc_mm", "Ecc_deg"] + list(density_fn.keys())
    rows = []
    radius_list = []  # ★追加
    area_list = []

    for (x, y), r in zip(pts_dict["Ganglion_Lee"]["pts"], radius_gc):
        rf = r_fn(r)
        radius_list.append(rf * 1000)  # mm → μm
        area_um2 = np.pi * (rf * 1000) ** 2
        area_list.append(area_um2)

        counts = []
        for cell in density_fn.keys():
            tree = pts_dict[cell]["tree"]
            counts.append(0 if tree is None
                          else len(tree.query_ball_point([x, y], rf)))
        rows.append([r, mm2deg(r), *counts])

    df = pd.DataFrame(rows, columns=cols)
    df["Radius_um"] = radius_list      # ★追加
    df["Area_um2"] = area_list
    df["Ecc_deg_bin"] = df["Ecc_deg"].round(2)
    gc_counts = df.groupby("Ecc_deg_bin").size()

    agg = OrderedDict()
    for cell in density_fn.keys():
        agg[f"Mean_{cell}"] = (cell, "mean")
        agg[f"SD_{cell}"]   = (cell, "std")
    agg["Mean_Radius_um"] = ("Radius_um", "mean")  # ★追加
    agg["SD_Radius_um"]   = ("Radius_um", "std")   # ★追加
    agg["Mean_Area_um2"]  = ("Area_um2", "mean")
    agg["SD_Area_um2"]    = ("Area_um2", "std")

    stat = (
        df.groupby("Ecc_deg_bin")
          .agg(**agg)
          .reset_index()
          .rename(columns={"Ecc_deg_bin": "Ecc_deg"})
    )
    stat["n_GC"] = stat["Ecc_deg"].map(gc_counts)
    stat.insert(0, "Ecc_mm", stat["Ecc_deg"] * MM_PER_DEG)
    stat.sort_values("Ecc_deg").to_excel(writer, sheet_name=sheet, index=False)


# with pd.ExcelWriter("GC_RF_CellStats_1-15mm_2.xlsx", engine="openpyxl") as writer:
#     build_stats(r_midget,  "Midget",  writer)
#     build_stats(r_parasol, "Parasol", writer)

# print("GC_RF_CellStats_1-15mm.xlsx 保存完了")

# # ------------------------------------------------------------
# # 6. Ganglion_Lee の数のグラフ作成
# # ------------------------------------------------------------
# import matplotlib.pyplot as plt

# # 半径ごとのGanglion_Lee密度を計算
# r_vals = np.linspace(0, 15, 300)  # 10–15 mmを100点で
# dens_fn = density_fn["Ganglion_Lee"]

# # 1mm幅のリング面積 × 密度 = 細胞数（近似）
# counts = []
# for r in r_vals:
#     ring_area = 2 * np.pi * r * 1.0  # 1 mm幅のリング
#     counts.append(dens_fn(r) * ring_area)

# # ---- 通常の 0–15 mm グラフ ----
# plt.figure(figsize=(7,5))
# plt.plot(r_vals, counts, color="blue", label="Ganglion_Lee")
# plt.xlabel("Eccentricity (mm)")
# plt.ylabel("Cell count per 1mm annulus")
# plt.title("Ganglion_Lee Cell Count (0–15 mm)")
# plt.legend()
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.savefig("Ganglion_Lee_CellCount_0-15mm.pdf")
# plt.close()

# # ---- 左右対称 -15〜15 mm グラフ ----
# x_vals = np.concatenate([-r_vals[::-1], r_vals])
# y_vals = np.concatenate([counts[::-1], counts])

# plt.figure(figsize=(9,5))
# plt.plot(x_vals, y_vals, color="red", label="Ganglion_Lee (symmetric)")
# plt.xlabel("Eccentricity (mm)")
# plt.ylabel("Cell count per 1mm annulus")
# plt.title("Ganglion_Lee Cell Count (-15 to 15 mm)")
# plt.legend()
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.savefig("Ganglion_Lee_CellCount_-15to15mm.pdf")
# plt.close()

# print("Ganglion_Lee の細胞数グラフを 0–15 mm と -15〜15 mm の pdf で保存しました。")
