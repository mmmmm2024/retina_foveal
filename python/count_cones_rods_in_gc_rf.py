"""
4子午線の視細胞密度データ（Cone/Rod）から、2D密度マップを補間し作成
指定したGC位置（deg）と受容野半径（midget/parasol近似）に対して、
受容野内のCone数・Rod数を積算して表示

入力: 
4meridians.xlsx（"Cones per sq mm", "Rods per sq mm"）
出力: 
ヒートマップ表示
受容野内のCone/Rod推定数（print）
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# データ読み込み
file_path = "4meridians.xlsx"         # Cones / Rods の 4子午線データ
cones = pd.read_excel(file_path, sheet_name="Cones per sq mm", header=7)
rods  = pd.read_excel(file_path, sheet_name="Rods per sq mm",  header=7)

# 使う列だけ抽出し名前を統一
cones = cones.iloc[:, [1,2,3,4,5]] ; rods = rods.iloc[:, [1,2,3,4,5]]
cones.columns = rods.columns = ["Ecc_mm","Sup","Inf","Temp","Nasal"]

# mm↔deg 変換
def mm_to_deg(r_mm):
    return 3.556*r_mm + 0.05993*r_mm**2 - 0.007358*r_mm**3 + 0.0003027*r_mm**4
def deg_to_mm(theta_deg):                   # 近似逆変換
    return 0.2762*theta_deg - 0.0001937*theta_deg**2 + 1.28e-5*theta_deg**3

# 4子午線 → 2D ヒートマップ補間
offset = {"Temp":1.5,"Sup":0.5,"Nasal":-1.5,"Inf":-0.5}
def interp_points(df):
    P,V = [],[]
    for _,row in df.iterrows():
        ecc = row["Ecc_mm"]
        for d,off in offset.items():
            val = row[d]
            if np.isnan(val): continue
            ecc_deg = mm_to_deg(ecc+off)-mm_to_deg(off)
            x,y = {"Temp":(-ecc_deg,0),"Nasal":(ecc_deg,0),
                   "Sup":(0,ecc_deg),"Inf":(0,-ecc_deg)}[d]
            P.append((x,y)); V.append(val)
    return np.array(P),np.array(V)

cone_P, cone_V = interp_points(cones)
rod_P,  rod_V  = interp_points(rods)

#　グリッド補間
g_x,g_y = np.mgrid[-60:60:500j, -60:60:500j]      # 500×500, 1セル≈0.24°
cone_grid = griddata(cone_P, cone_V, (g_x,g_y), method="cubic", fill_value=0)
rod_grid  = griddata(rod_P,  rod_V,  (g_x,g_y), method="cubic", fill_value=0)

# 任意 GC の設定
theta_x, theta_y = 10, 0              # GC中心 (deg)
view_dist_mm     = 500                # スクリーンまでの距離

def rf_radius_deg(ecc_mm, gc_type="parasol"):
    if gc_type=="parasol":
        return 70.2*ecc_mm**0.6/2/1000
    return 8.64*ecc_mm**1.04/2/1000   # midget

ecc_mm = np.hypot(deg_to_mm(theta_x), deg_to_mm(theta_y))
rf_deg = rf_radius_deg(ecc_mm,"parasol")   # ← midget にする場合は引数変更

# マスク & 細胞数積算
mask = (g_x-theta_x)**2 + (g_y-theta_y)**2 <= rf_deg**2

deg_step = 120/500                       # 1セル辺長 (deg)
mm_per_deg = view_dist_mm*np.pi/180
cell_area_mm2 = (deg_step*mm_per_deg)**2

cone_num = np.sum(cone_grid[mask]*cell_area_mm2)
rod_num  = np.sum( rod_grid[mask]*cell_area_mm2)

rf_radius_mm = view_dist_mm*np.tan(np.deg2rad(rf_deg))
rf_area_mm2  = np.pi*rf_radius_mm**2
cone_area_per = rf_area_mm2/cone_num if cone_num>0 else None
rod_area_per  = rf_area_mm2/rod_num  if rod_num>0 else None

# 可視化
plt.figure(figsize=(7,7))
plt.contourf(g_x,g_y,cone_grid,levels=100,cmap="viridis")
plt.colorbar(label="Cone density (cells/mm²)")
circ = plt.Circle((theta_x,theta_y),rf_deg,color="white",fill=False,lw=2)
plt.gca().add_patch(circ)
plt.scatter(theta_x,theta_y,color="red",s=30)
plt.title(f"GC RF @({theta_x}°, {theta_y}°)\nCones:{int(cone_num)}, Rods:{int(rod_num)}")
plt.axis("equal"); plt.xlabel("Horizontal (deg)"); plt.ylabel("Vertical (deg)")
plt.tight_layout(); plt.show()

# print-out
print(f"Cones in RF : {cone_num:.0f}")
print(f"Rods  in RF : {rod_num:.0f}")
print(f"Area per cone (mm²) : {cone_area_per:.3f}")
print(f"Area per rod  (mm²) : {rod_area_per:.3f}")
