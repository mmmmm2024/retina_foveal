"""
密度モデル（Rod / AIIAC / Ganglion_Lee）を偏心度範囲で評価
Rod密度の値や比（例: dens(1mm) と dens(12mm) の比、密度から面積換算した比）を出力

処理:
- Leeらの細胞の密度関数
- 偏心度 r(mm) を入力して Rod 密度 dens_r(r) を計算
- dens_r(1), dens_r(12) やその比率・換算値を出力

"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from collections import OrderedDict
import matplotlib.pyplot as plt

np.random.seed(0)

ARC_W   = 0.1
R_MIN   = 11.0
R_MAX   = 13.0
RADII   = np.arange(R_MIN, R_MAX, ARC_W)

MM_PER_DEG = 0.29
mm2deg = lambda r_mm: r_mm / MM_PER_DEG

# 樹状突起半径の倍率
RADIUS_SCALE = 1.2   # 例: 0.8, 1.2, 2.0 など


# Midget GC 樹状突起半径モデル（μm → mm）
def r_midget(r_mm):
    return (8.64 * r_mm**1.04) / 2 / 1000  # mm


# 密度モデル
def diff_exp(x, c1, k1, c2, k2, c3=0, k3=0):
    return c1*np.exp(k1*x) + c2*np.exp(k2*x) + c3*np.exp(k3*x)

PARAMS = OrderedDict([
    ("Ganglion_Lee", (8.475e5, -1.258e-1, -8.691e5, -1.556e0)),
    ("AIIAC",        (1.297e5, -1.114e0,  -1.392e5, -1.217e0, 2.553e3, -7.214e-2)),
    ("Rod",          (5.855e5, -1.388e-1, -5.989e5, -2.998e-1)),
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

# ------------------------------------------------------------
# 5. 統計計算 → Excel（Midget GC 専用 + 半径倍率適用）
# ------------------------------------------------------------
# def build_stats_midget(sheet, writer, radius_scale=1.0):
#     cols = ["Ecc_mm", "Ecc_deg"] + list(density_fn.keys())
#     rows = []
#     radius_list = []
#     area_list = []

#     gc_pts = pts_dict["Ganglion_Lee"]["pts"]
#     if gc_pts.size == 0:
#         # GCが1つも生成されなかった場合でも落ちないように
#         empty = pd.DataFrame(columns=cols + ["Radius_um", "Area_um2"])
#         empty.to_excel(writer, sheet_name=sheet, index=False)
#         return

#     for (x, y), r in zip(gc_pts, radius_gc):
#         rf = radius_scale * r_midget(r)   # ★ここで倍率を反映（rf: mm）
#         radius_um = rf * 1000             # mm → μm
#         area_um2  = np.pi * (radius_um ** 2)

#         radius_list.append(radius_um)
#         area_list.append(area_um2)

#         counts = []
#         for cell in density_fn.keys():
#             tree = pts_dict[cell]["tree"]
#             counts.append(0 if tree is None else len(tree.query_ball_point([x, y], rf)))
#         rows.append([r, mm2deg(r), *counts])

#     df = pd.DataFrame(rows, columns=cols)
#     df["Radius_um"] = radius_list
#     df["Area_um2"]  = area_list

#     df["Ecc_deg_bin"] = df["Ecc_deg"].round(2)
#     gc_counts = df.groupby("Ecc_deg_bin").size()

#     agg = OrderedDict()
#     for cell in density_fn.keys():
#         agg[f"Mean_{cell}"] = (cell, "mean")
#         agg[f"SD_{cell}"]   = (cell, "std")
#     agg["Mean_Radius_um"] = ("Radius_um", "mean")
#     agg["SD_Radius_um"]   = ("Radius_um", "std")
#     agg["Mean_Area_um2"]  = ("Area_um2", "mean")
#     agg["SD_Area_um2"]    = ("Area_um2", "std")

#     stat = (
#         df.groupby("Ecc_deg_bin")
#           .agg(**agg)
#           .reset_index()
#           .rename(columns={"Ecc_deg_bin": "Ecc_deg"})
#     )
#     stat["n_GC"] = stat["Ecc_deg"].map(gc_counts)
#     stat.insert(0, "Ecc_mm", stat["Ecc_deg"] * MM_PER_DEG)

#     # 仕上げ
#     stat = stat.sort_values("Ecc_deg")
#     stat.to_excel(writer, sheet_name=sheet, index=False)

#     # （任意）生データも欲しいならコメント解除
#     # df.to_excel(writer, sheet_name=f"{sheet}_raw", index=False)

# out_xlsx = f"GC_RF_CellStats_Midget_10-15mm_scale{RADIUS_SCALE:.2f}.xlsx"
# with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
#     build_stats_midget("Midget", writer, radius_scale=RADIUS_SCALE)

# print(f"保存完了: {out_xlsx}")

# ------------------------------------------------------------
# 6. Ganglion_Lee の数のグラフ作成（任意で残す）
# ------------------------------------------------------------
# r_vals = np.linspace(0, 15, 300)
# dens_fn = density_fn["Ganglion_Lee"]

# counts = []
# for r in r_vals:
#     ring_area = 2 * np.pi * r * 1.0  # 1 mm 幅リング
#     counts.append(dens_fn(r) * ring_area)

# plt.figure(figsize=(7,5))
# plt.plot(r_vals, counts, label="Ganglion_Lee")
# plt.xlabel("Eccentricity (mm)")
# plt.ylabel("Cell count per 1mm annulus")
# plt.title("Ganglion_Lee Cell Count (0–15 mm)")
# plt.grid(True, linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.savefig("Ganglion_Lee_CellCount_0-15mm.pdf")
# plt.close()

# print("Ganglion_Lee の細胞数グラフ(pdf)も保存しました。")

# def dens_gc(r_mm: float) -> float:
#     return float(density_fn["Ganglion_Lee"](r_mm))
# print(dens_gc(1))
# print(dens_gc(12))
# print(40/dens_gc(1))
# print(40/dens_gc(12))
# print(40 * dens_gc(12)/dens_gc(1))

def dens_r(r_mm: float) -> float:
    return float(density_fn["Rod"](r_mm))
print(dens_r(1))
print(dens_r(12))
print(400/dens_r(1))
print(2520/dens_r(12))
print(120960/dens_r(12))
print(40 * dens_r(12)/dens_r(1))