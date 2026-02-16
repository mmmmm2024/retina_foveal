"""
RGC（神経節細胞）密度の「全体マップ」を作るプログラム
----------------------------------------------------
実行:
python rgc_density_polarmap_curcio_watson.py

Excelに入っているCurcio(1990)のRGC密度データ（4方向の子午線データ）を
Watsonの式で補正したうえで、網膜全体の密度分布を推定
極座標のヒートマップとしてPDFに保存

- Excelから、Temporal / Superior / Nasal / Inferior の4方向データを読む

- 偏心度を mm → deg（視野角）に変換する（Watson Appendix A6 の式）

- 「光軸と視軸のズレ（offset）」を方向ごとに補正して、
    視野上での位置（deg）になるように調整する
    ※Temporal/Superior/Nasal/Inferiorで補正量が違う

- 面積補正係数（Watson Appendix A7）を使って密度を補正する
（mm^2あたりの密度を、deg^2に対応する形へ変換するイメージ）

- 各方向のデータ点を、極座標(r,θ)から(x,y)に変換して平面上に配置する

- 4方向にしか点がないので、griddata(method="cubic")で2次元補間を行い、
網膜全体の密度マップ（連続的な面のデータ）を作る

- 極座標のヒートマップとして表示し、PDFに保存する

入力:
- Curcio_JCompNeurol1990_GCtopo_F6.xlsx
    必要な列:
    mm
    Temp.mean_GC/sq mm, Supe.mean_GC/sq mm, Nasa.mean_GC/sq mm, Infe.mean_GC/sq mm
    ※列名は "Temp/Supe/Nasa/Infe" の4文字に対応している必要があります

出力:
- rgc_density_cubic.pdf（極座標の密度ヒートマップ）
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Watson Appendix A6: mm → deg
def mm_to_deg(r_mm):
    return (
        3.556 * r_mm +
        0.05993 * r_mm**2 -
        0.007358 * r_mm**3 +
        0.0003027 * r_mm**4
    )

# Watson Appendix A7: 面積補正係数
def area_conversion_factor(r_deg):
    return (
        0.0752 +
        5.846e-5 * r_deg -
        1.064e-5 * r_deg**2 +
        4.116e-8 * r_deg**3
    )

# 光軸と視軸のズレ D (mm)
offsets = {
    'Temporal': 1.5,
    'Superior': 0.5,
    'Nasal': -1.5,
    'Inferior': -0.5
}

# 方位に対応する角度（度）
angle_map = {
    'Temporal': 0,
    'Superior': 90,
    'Nasal': 180,
    'Inferior': 270
}

# Excelファイル読み込み
file_path = "Curcio_JCompNeurol1990_GCtopo_F6.xlsx"
data = pd.read_excel(file_path)

# データ格納用
all_points = []
all_values = []
all_theta = []
all_r = []

# 各子午線処理
for key in offsets.keys():
    ecc_mm = data['mm']
    raw_density = pd.to_numeric(data[f'{key[:4]}.mean_GC/sq mm'], errors='coerce')
    
    r_optic = ecc_mm - offsets[key]
    deg_optic = mm_to_deg(r_optic)
    deg_corr = mm_to_deg(offsets[key])
    ecc_deg_corrected = deg_optic + deg_corr

    area_factor = area_conversion_factor(ecc_deg_corrected)
    corrected_density = raw_density * area_factor

    # 有効データのみ
    valid = ~np.isnan(ecc_deg_corrected) & ~np.isnan(corrected_density)
    r_vals = ecc_deg_corrected[valid].values
    d_vals = corrected_density[valid].values
    theta_rad = np.radians(angle_map[key])

    # 極→直交変換
    x_vals = r_vals * np.cos(theta_rad)
    y_vals = r_vals * np.sin(theta_rad)

    all_points.extend(np.stack([x_vals, y_vals], axis=1))
    all_values.extend(d_vals)
    all_theta.extend([theta_rad] * len(r_vals))
    all_r.extend(r_vals)

# numpyに変換
points = np.array(all_points)
values = np.array(all_values)
r_all = np.array(all_r)
theta_all = np.array(all_theta)

# --- 極座標補間 ---
r_grid = np.linspace(0.001, 30, 300)
theta_grid = np.linspace(-np.pi, np.pi, 300)
theta_mesh, r_mesh = np.meshgrid(theta_grid, r_grid)
x_polar = r_mesh * np.cos(theta_mesh)
y_polar = r_mesh * np.sin(theta_mesh)
polar_z = griddata(points, values, (x_polar, y_polar), method='cubic')

# --- 直交座標補間 ---
# grid_x, grid_y = np.meshgrid(np.linspace(-60, 60, 300), np.linspace(-60, 60, 300))
# grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')

# --- カラーマップ設定（NaNは白） ---
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='white')

# --- 描画 ---
fig = plt.figure(figsize=(12, 12))

# 極座標
ax1 = fig.add_subplot(polar=True) #1, 2, 1, polar=True)
pc = ax1.pcolormesh(theta_mesh, r_mesh, polar_z, shading='auto', cmap=cmap) #,
                    #vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z))
# 方向ラベルを描画（r=最大の位置に表示）
label_r = r_grid.max() + 5 # 半径の外側に文字を置く

# 各方向にテキスト追加
for label, angle_deg in angle_map.items():
    angle_rad = np.radians(angle_deg)
    ax1.text(angle_rad, label_r, label, ha='center', va='center',
             fontsize=10, fontweight='bold', color='black')
# ax1.set_title("RGC Density (Polar)")
plt.colorbar(pc, ax=ax1, pad=0.15, shrink=0.9)

# 直交座標
# ax2 = fig.add_subplot(1, 2, 2)
# im = ax2.imshow(grid_z, extent=(-30, 30, -30, 30), origin='lower', #cmap=cmap,
#                 vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z))
# ax2.set_title("RGC Density (Cartesian)")
# ax2.set_xlabel("X (deg)")
# ax2.set_ylabel("Y (deg)")
# plt.colorbar(im, ax=ax2)

plt.tight_layout()
plt.savefig("rgc_density_cubic.pdf")