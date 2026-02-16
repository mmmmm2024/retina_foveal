"""
視野座標(deg)に10°刻みでGC中心を格子状に配置し、
各位置の偏心度(mm)から midget/parasol の受容野半径を計算して、
受容野の円を刺激平面(mm, 視距離500mm)に投影して可視化する。
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

# 描画開始
plt.figure(figsize=(10, 10))

# θの配置範囲（a, b = -60～60、10ずつ）
for a in range(-6, 7):
    for b in range(-6, 7):
        theta_center_x = 10 * a
        theta_center_y = 10 * b

        # 偏心度の計算
        ecc_mm = compute_eccentricity_mm(theta_center_x, theta_center_y)

        # Parasolの外周
        r_parasol_deg = parasol_rf_radius_deg(ecc_mm)
        angles = np.linspace(0, 2 * np.pi, 100)
        theta_x_edge = theta_center_x + mm_to_deg(r_parasol_deg) * np.cos(angles)
        theta_y_edge = theta_center_y + mm_to_deg(r_parasol_deg) * np.sin(angles)
        x_parasol, y_parasol = project_to_visual_field(theta_x_edge, theta_y_edge)
        plt.plot(x_parasol, y_parasol, color='navy', alpha=0.3, label='Parasol' if (a == -6 and b == -6) else "")

        # Midgetの外周
        r_midget_deg = midget_rf_radius_deg(ecc_mm)
        theta_x_edge_m = theta_center_x + mm_to_deg(r_midget_deg) * np.cos(angles)
        theta_y_edge_m = theta_center_y + mm_to_deg(r_midget_deg) * np.sin(angles)
        x_midget, y_midget = project_to_visual_field(theta_x_edge_m, theta_y_edge_m)
        plt.plot(x_midget, y_midget, color='darkorange', alpha=0.3, label='Midget' if (a == -6 and b == -6) else "")

        # 中心点（GCの位置）
        cx, cy = project_to_visual_field(theta_center_x, theta_center_y)
        plt.scatter(cx, cy, color='red', s=5)

# 軸や凡例などの整備
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
