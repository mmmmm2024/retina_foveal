"""
4子午線のCone/Rod密度データ(4meridians.xlsx)から2D密度マップを補間
指定GC位置(θx,θy)の受容野内Cone/Rod数と、投影受容野面積(mm²)から
「1細胞あたり刺激面積(area_per_cone/rod)」を計算
さらに視野格子上の多数GCについて同じ量を計算し、ExcelとヒートマップPDFを出力

入力: 4meridians.xlsx（Cones per sq mm / Rods per sq mm）
出力: Cone_Rod_Heatmaps.pdf, GC_10000_summary.xlsx, Area_per_Cell_Heatmaps_10000.pdf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def mm2deg(r):
    return 3.556*r + 0.05993*r**2 - 7.358e-3*r**3 + 3.027e-4*r**4

def deg2mm(t):
    return 0.2762*t - 1.937e-4*t**2 + 1.28e-5*t**3

offset = {"Temp":1.5, "Sup":0.5, "Nasal":-1.5, "Inf":-0.5}
def four_meridians_to_points(df):
    P, V = [], []
    for _, row in df.iterrows():
        ecc = row["Ecc_mm"]
        for d, off in offset.items():
            val = row[d]
            if np.isnan(val): continue
            ecc_deg = mm2deg(ecc + off) - mm2deg(off)
            x, y = {"Temp":(-ecc_deg,0), "Nasal":(ecc_deg,0),
                    "Sup":(0,ecc_deg),   "Inf": (0,-ecc_deg)}[d]
            P.append((x, y)); V.append(val)
    return np.array(P), np.array(V)

def rf_radius_deg(ecc_mm, t):
    diam_um = 70.2*ecc_mm**0.6 if t=="parasol" else 8.64*ecc_mm**1.04
    return mm2deg(diam_um / 2000)

# Rod / Cone ヒートマップ
xls = "4meridians.xlsx"
cones = pd.read_excel(xls, sheet_name="Cones per sq mm", header=7).iloc[:, [1,2,3,4,5]]
rods  = pd.read_excel(xls, sheet_name="Rods per sq mm",  header=7).iloc[:, [1,2,3,4,5]]
cones.columns = rods.columns = ["Ecc_mm","Sup","Inf","Temp","Nasal"]

cone_P, cone_V = four_meridians_to_points(cones)
rod_P , rod_V  = four_meridians_to_points(rods)

gx, gy = np.mgrid[-100:100:500j, -100:100:500j]
cone_map = griddata(cone_P, cone_V, (gx, gy), method="cubic", fill_value=0)
rod_map  = griddata(rod_P , rod_V , (gx, gy), method="cubic", fill_value=0)

# ----- 可視化 & PDF 保存 -----
fig = plt.figure(figsize=(12, 5))
for i, (data, title, cmap) in enumerate([(cone_map, "Cone density", "viridis"),
                                         (rod_map , "Rod  density", "viridis")]):
    ax = fig.add_subplot(1, 2, i + 1)
    ax.contourf(gx, gy, data, levels=100, cmap=cmap)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("Horizontal (deg)")
    ax.set_ylabel("Vertical (deg)")
    fig.colorbar(ax.collections[0], ax=ax, label="cells / mm²")
fig.tight_layout()
fig.savefig("Cone_Rod_Heatmaps.pdf")

# -------------------------------------------------
# 2–5. 個別 GC 計算 (例: 40°,40°)
# -------------------------------------------------
theta_x, theta_y = 40, 40
gc_type = "midget"
view_mm = 500

ecc_mm = np.hypot(deg2mm(theta_x), deg2mm(theta_y))
rf_deg = rf_radius_deg(ecc_mm, gc_type)
mask   = (gx - theta_x)**2 + (gy - theta_y)**2 <= rf_deg**2

deg_cell = 120 / 500
mm_per_deg = view_mm * np.pi / 180
cell_area_mm2 = (deg_cell * mm_per_deg)**2
cone_n = np.sum(cone_map[mask] * cell_area_mm2)
rod_n  = np.sum( rod_map[mask] * cell_area_mm2)

phi = np.linspace(0, 2*np.pi, 360)
edge_tx = theta_x + rf_deg*np.cos(phi)
edge_ty = theta_y + rf_deg*np.sin(phi)
edge_xmm = view_mm * np.tan(np.radians(edge_tx))
edge_ymm = view_mm * np.tan(np.radians(edge_ty))
stim_area = 0.5 * np.abs(np.dot(edge_xmm, np.roll(edge_ymm, -1)) -
                         np.dot(edge_ymm, np.roll(edge_xmm, -1)))
area_per_cone = stim_area / cone_n if cone_n else np.nan
area_per_rod  = stim_area / rod_n  if rod_n  else np.nan

print(f"GC @ ({theta_x}°, {theta_y}°)  [{gc_type}]")
print(f" RF angle radius  : {rf_deg:.3f}°")
print(f" RF projected area: {stim_area:.3f} mm²")
print(f" Cones in RF      : {cone_n:.0f}")
print(f" Rods  in RF      : {rod_n :.0f}")
print(f" Area / cone (mm²): {area_per_cone:.4e}")
print(f" Area / rod  (mm²): {area_per_rod :.4e}")


# (-50°,-50°) 〜 (+49°,+49°) 100×100＝10 000 GC を Excel 出力
rows_10k = []
x_vals = np.arange(-70, 70)          # 100 等間隔 (-50 … +49)
y_vals = np.arange(-70, 70)

for tx in x_vals:
    for ty in y_vals:
        ecc_mm  = np.hypot(deg2mm(tx), deg2mm(ty))
        rf_d    = rf_radius_deg(ecc_mm, gc_type)           # gc_type = "parasol" or "midget"
        mask    = (gx - tx)**2 + (gy - ty)**2 <= rf_d**2   # ギザギザ近似

        cone_cnt = np.sum(cone_map[mask] * cell_area_mm2)
        rod_cnt  = np.sum( rod_map[mask] * cell_area_mm2)

        phi = np.linspace(0, 2*np.pi, 360)
        ex  = tx + rf_d*np.cos(phi)
        ey  = ty + rf_d*np.sin(phi)
        ex_mm = view_mm * np.tan(np.radians(ex))
        ey_mm = view_mm * np.tan(np.radians(ey))
        area  = 0.5 * np.abs(np.dot(ex_mm, np.roll(ey_mm, -1)) -
                             np.dot(ey_mm, np.roll(ex_mm, -1)))

        rows_10k.append({
            "theta_x_deg":   tx,
            "theta_y_deg":   ty,
            "ecc_mm":        ecc_mm,
            "rf_angle_deg":  rf_d,
            "rf_area_mm2":   area,
            "cones":         cone_cnt,
            "rods":          rod_cnt,
            "area_per_cone": area / cone_cnt if cone_cnt else np.nan,
            "area_per_rod":  area / rod_cnt  if rod_cnt  else np.nan
        })

pd.DataFrame(rows_10k).to_excel("GC_10000_summary.xlsx", index=False)

# ---- 2500 GC の「1細胞あたり刺激面積」ヒートマップ ----------------
# 既に rows → DataFrame にまとめてある場合はそれを再利用
# 例: gc_df = pd.read_excel("GC_2500_summary.xlsx")

gc_df = pd.DataFrame(rows_10k)          # rows がメモリにある場合

# 50×50 グリッド（x=0..49, y=0..49）
X = np.arange(70)
Y = np.arange(70)
Z_cone = np.full((70, 70), np.nan)
Z_rod  = np.full((70, 70), np.nan)

for _, r in gc_df.iterrows():
    ix, iy = int(r["theta_x_deg"]), int(r["theta_y_deg"])
    Z_cone[iy, ix] = r["area_per_cone"]
    Z_rod [iy, ix] = r["area_per_rod"]

# ヒートマップ (–50°..+49°, 100×100) 可視化
gc10k_df = pd.DataFrame(rows_10k)           # 先ほど作成した 10 000 GC データ

# 2-D 配列へ整形
pivot_cone = gc10k_df.pivot(index="theta_y_deg",
                            columns="theta_x_deg",
                            values="area_per_cone").sort_index(ascending=True)
pivot_rod  = gc10k_df.pivot(index="theta_y_deg",
                            columns="theta_x_deg",
                            values="area_per_rod").sort_index(ascending=True)

Z_cone = pivot_cone.to_numpy()              # shape = (100,100)
Z_rod  = pivot_rod .to_numpy()

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

pcm1 = ax[0].imshow(Z_cone, origin="lower", cmap="plasma",
                    extent=[-70, 69, -70, 69], aspect="equal")
ax[0].set_title("Area per cone (mm²)")
ax[0].set_xlabel("theta_x (deg)")
ax[0].set_ylabel("theta_y (deg)")
fig.colorbar(pcm1, ax=ax[0])

pcm2 = ax[1].imshow(Z_rod, origin="lower", cmap="viridis",
                    extent=[-70, 69, -70, 69], aspect="equal")
ax[1].set_title("Area per rod (mm²)")
ax[1].set_xlabel("theta_x (deg)")
ax[1].set_ylabel("theta_y (deg)")
fig.colorbar(pcm2, ax=ax[1])

plt.tight_layout()
plt.savefig("Area_per_Cell_Heatmaps_10000.pdf")   # PDF に保存
plt.show()
