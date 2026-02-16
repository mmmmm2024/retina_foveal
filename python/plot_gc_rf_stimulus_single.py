"""
指定したGC位置(θx, θy [deg])について、偏心度(mm)を計算
midget/parasolの受容野半径を求めて、刺激平面(視距離500mm)へ投影した受容野円を描画
"""

import numpy as np
import matplotlib.pyplot as plt

# 偏心角 (deg) → 偏心度 (mm)（Watson 2014式）
def deg_to_mm(theta_deg):
    return 0.2762 * theta_deg - 0.0001937 * theta_deg**2 + 1.280e-5 * theta_deg**3

# mm → deg 変換関数
def mm_to_deg(r_mm):
    return (
        3.556 * r_mm +
        0.05993 * r_mm**2 -
        0.007358 * r_mm**3 +
        0.0003027 * r_mm**4
    )

# 偏心角から偏心度（mm）を計算
def compute_eccentricity_mm(theta_x, theta_y):
    x_mm = deg_to_mm(theta_x)
    y_mm = deg_to_mm(theta_y)
    return np.sqrt(x_mm**2 + y_mm**2)

# Parasol型GCの受容野半径（mm）
def parasol_rf_radius_deg(ecc_mm):
    return 70.2 * ecc_mm**0.6 / 2 / 1000

# Midget型GCの受容野半径（mm）
def midget_rf_radius_deg(ecc_mm):
    return 8.64 * ecc_mm**1.04 / 2 / 1000

# 視覚刺激空間への変換（500mm離れた位置に投影）
def project_to_visual_field(theta_x, theta_y, viewing_distance=500):
    x = viewing_distance * np.tan(np.radians(theta_x))
    y = viewing_distance * np.tan(np.radians(theta_y))
    return x, y

# 任意のGC位置（網膜偏心角）を指定（例：30°, 20°）
theta_center_x = 40
theta_center_y = 40

# 偏心度（mm）を計算し、Parasol型とMidget型の半径を取得
ecc_mm = compute_eccentricity_mm(theta_center_x, theta_center_y)
r_parasol_deg = parasol_rf_radius_deg(ecc_mm)
r_midget_deg = midget_rf_radius_deg(ecc_mm)

# 受容野の境界（角度）を生成
angles = np.linspace(0, 2 * np.pi, 360)
# Parasol
theta_x_edge_p = theta_center_x + mm_to_deg(r_parasol_deg) * np.cos(angles)
theta_y_edge_p = theta_center_y + mm_to_deg(r_parasol_deg) * np.sin(angles)
# Midget
theta_x_edge_m = theta_center_x + mm_to_deg(r_midget_deg) * np.cos(angles)
theta_y_edge_m = theta_center_y + mm_to_deg(r_midget_deg) * np.sin(angles)

# 500mm離れた視空間に投影
x_visual_p, y_visual_p = project_to_visual_field(theta_x_edge_p, theta_y_edge_p)
x_visual_m, y_visual_m = project_to_visual_field(theta_x_edge_m, theta_y_edge_m)
cx, cy = project_to_visual_field(theta_center_x, theta_center_y)

# 可視化（視軸を原点とした平面）
plt.figure(figsize=(6, 6))
plt.plot(x_visual_p, y_visual_p, label='Parasol', color='navy')
plt.plot(x_visual_m, y_visual_m, label='Midget', color='darkorange')
plt.scatter(cx, cy, color='red', label='GC')
plt.axhline(0, color='gray')
plt.axvline(0, color='gray')
plt.axis('equal')

plt.xlabel('transverse')
plt.ylabel('sagittal')
plt.title('Visual Stimulus Areas of GC Receptive Fields')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
