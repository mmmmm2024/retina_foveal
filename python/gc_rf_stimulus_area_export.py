"""
Curcio_RGC_Xvalues.xlsx の Ecc_mm から midget/parasol の受容野半径を計算
mm→deg（Watson式）→刺激平面(mm, 視距離500mm)へ変換して受容野面積/半径を求めてCSV出力
各GCの受容野を円で可視化（表示のみ）。

出力: GC_Midget_StimulusAreas.csv / GC_Parasol_StimulusAreas.csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 視度変換関数（mm → deg）
def mm_to_deg(r_mm):
    return (
        3.556 * r_mm +
        0.05993 * r_mm**2 -
        0.007358 * r_mm**3 +
        0.0003027 * r_mm**4
    )

# 受容野のdeg → 平面刺激距離（mm）に変換 
def deg_to_stimulus_mm(radius_deg, viewing_distance_mm=500):
    radius_rad = np.deg2rad(radius_deg)
    return abs(viewing_distance_mm * np.tan(radius_rad))

# 平面刺激上の受容野面積を計算
def calculate_stimulus_area(r_gc_mm, r_dend_mm, viewing_distance_mm=500):
    r_min = max(r_gc_mm - r_dend_mm, 0)
    r_max = r_gc_mm + r_dend_mm
    deg_min = mm_to_deg(r_min)
    deg_max = mm_to_deg(r_max)
    x_min = deg_to_stimulus_mm(deg_min, viewing_distance_mm)
    x_max = deg_to_stimulus_mm(deg_max, viewing_distance_mm)
    area = np.pi * abs(x_max**2 - x_min**2)
    center_x = deg_to_stimulus_mm(mm_to_deg(r_gc_mm), viewing_distance_mm)
    return center_x, area

# 樹状突起直径式（μm）
def dendritic_radius_midget(r_mm):
    return (8.64 * (r_mm ** 1.04)) / 2 / 1000  # → mm

def dendritic_radius_parasol(r_mm):
    return (70.2 * (r_mm ** 0.6)) / 2 / 1000  # → mm

df = pd.read_excel("Curcio_RGC_Xvalues.xlsx")
df = df.dropna(subset=["Ecc_mm", "X_mm"])  # 基本列があるものに限定


def compute_gc_stimulus_areas(df, radius_function, output_csv):
    rows = []
    for _, row in df.iterrows():
        r_gc = row["Ecc_mm"]
        x_mm = row["X_mm"]
        r_dend = radius_function(r_gc)

        # 偏心度（deg）
        ecc_deg = mm_to_deg(r_gc) - mm_to_deg(1.5)

        # 受容野の範囲（mm）
        r_min = max(r_gc - r_dend, 0)
        r_max = r_gc + r_dend

        # 受容野の範囲（deg）
        deg_min = mm_to_deg(r_min)
        deg_max = mm_to_deg(r_max)

        # 平面刺激での中心・面積
        x_center, area = calculate_stimulus_area(r_gc, r_dend)
        
        # 面積 → 半径
        radius = np.sqrt(area / np.pi)
        
        rows.append([
            r_gc,         # 偏心度(mm)
            ecc_deg,      # 偏心度(deg)
            deg_min,      # 受容野最小角度
            deg_max,      # 受容野最大角度
            x_mm,         # 元のX位置
            x_center,     # 平面刺激上の中心
            area,          # 平面刺激上の受容野面積
            radius        # 半径（mm）
        ])
    
    df_out = pd.DataFrame(rows, columns=[
        "Eccentricity (mm)",
        "Eccentricity (deg)",
        "Receptive Field Min (deg)",
        "Receptive Field Max (deg)",
        "X_mm (from data)",
        "X_center (computed)",
        "Stimulus Area (mm²)",
        "radius"        # 半径

    ])
    
    df_out.to_csv(output_csv, index=False)
    return df_out


# Midget型
df_midget = compute_gc_stimulus_areas(df, dendritic_radius_midget, "GC_Midget_StimulusAreas.csv")

#  Parasol型 
df_parasol = compute_gc_stimulus_areas(df, dendritic_radius_parasol, "GC_Parasol_StimulusAreas.csv")


import matplotlib.pyplot as plt
import numpy as np

def plot_all_gc_receptive_fields(df_out, margin=100):
    """
    各GCの平面刺激受容範囲を、図形的に描画（はみ出さないよう自動スケーリング）。
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # 描画中心
    center_y = 0  # GCは横方向に並ぶ前提
    xs = []
    radii = []

    for _, row in df_out.iterrows():
        gc_x = row["X_center (computed)"]
        radius = np.sqrt(row["Stimulus Area (mm²)"] / np.pi)
        xs.append(gc_x)
        radii.append(radius)

        circle = plt.Circle((gc_x, center_y), radius, color='skyblue', alpha=0.6, edgecolor='k')
        ax.add_patch(circle)
        ax.plot(gc_x, center_y, 'o', color='blue')

    # 中心位置（原点）に × を描画
    ax.plot(0, center_y, 'x', color='black', markersize=10, label="Stimulus Center")

    # 自動スケール調整（左右端のGC＋半径に基づく）
    x_min = min(x - r for x, r in zip(xs, radii)) - margin
    x_max = max(x + r for x, r in zip(xs, radii)) + margin
    y_min = center_y - max(radii) - margin
    y_max = center_y + max(radii) + margin

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')
    ax.set_yticks([])
    # ax.set_xlabel("Stimulus plane X position (mm)")
    ax.set_title("Stimulus")
    # ax.legend()
    plt.tight_layout()
    plt.show()

plot_all_gc_receptive_fields(df_midget)    # Midget型GC
plot_all_gc_receptive_fields(df_parasol) # Parasol型GC（必要なら）