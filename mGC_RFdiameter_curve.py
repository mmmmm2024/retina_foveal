"""
midget GC（mGC）の受容野直径モデル（8.64 * x^1.04）を
偏心度 x(mm) = 0〜15 の範囲で計算してプロットし、PDFに保存。

- 青線: 受容野直径 y = 8.64 * x^1.04（μm）
- 赤破線: 直線近似の比較用 y = 8.64 * x（μm）
出力: ReceptiveField_Midget_0-15mm.pdf
"""

import numpy as np
import matplotlib.pyplot as plt

# mGCの受容野直径モデル（偏心度x[mm] → 直径[μm]）
def receptive_field_diameter(x):
    return 8.64 * x**1.04  # μm

# x の範囲 (0〜15mm, 0.01刻み)
x = np.linspace(0, 15, 500)
y = receptive_field_diameter(x)

# プロット
plt.figure(figsize=(7,5))
plt.plot(x, y, color="blue", linewidth=1, label="mGC")
# y = x の赤い補助線（単位はそのまま：mm 対 μm）
plt.plot(x, 8.64 * x, color="red", linestyle="--", linewidth=1, label=r"$y=8.64x$")

plt.xlabel("Eccentricity (mm)", fontsize=12)
plt.ylabel("Receptive field diameter (μm)", fontsize=12)
# plt.title("Midget GC Receptive Field Diameter", fontsize=14)
plt.legend()
plt.xlim(0, 15)
plt.ylim(0, 140)
plt.grid(True, linestyle="--", alpha=0.6)

# 保存
plt.tight_layout()
plt.savefig("ReceptiveField_Midget_0-15mm.pdf")
plt.show()
