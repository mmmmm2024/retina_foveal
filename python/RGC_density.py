"""
Watson(2014)のRGC密度モデルとCurcio(1990)実測データを読込
mm↔deg変換・面積補正・offset補正を適用し、密度/累積数/比較図を作成しPDF保存

出力PDF: 
    Watson密度(log-log)
    Curcio密度(log-log)
    Watson累積RGC数
    Watson vs Curcio比較
    deg↔mm変換
    offset補正の影響
    面積補正係数(A7)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Watson (2014) の数式モデルに基づく解析

# 子午線ごとの色指定（論文図と対応）
colors = {
    'Temporal': 'red',
    'Superior': 'blue',
    'Nasal': 'green',
    'Inferior': 'gray'
}

offsets = {
    'Temporal': 1.5,
    'Superior': 0.5,
    'Nasal': -1.5,
    'Inferior': -0.5
}

# 視野中心（fovea）での受容野密度（Watson論文より）
dgf0 = 33163.2  # cells / deg^2

# 子午線ごとのモデルパラメータ（Equation 4 の a, r2, re）
params = {
    'Temporal': {'a': 0.9851, 'r2': 1.058, 're': 22.14},
    'Superior': {'a': 0.9935, 'r2': 1.035, 're': 16.35},
    'Nasal':    {'a': 0.9729, 'r2': 1.084, 're': 7.633},
    'Inferior': {'a': 0.9960, 'r2': 0.9932, 're': 12.13}
}

# Watsonモデルによる密度プロット
r_log = np.logspace(-1, 2, 500)  # 偏心 0.1〜100 deg
plt.figure(figsize=(7, 5))
for name, p in params.items():
    a, r2, re = p['a'], p['r2'], p['re']
    dgf = dgf0 * (a * (1 + r_log / r2)**-2 + (1 - a) * np.exp(-r_log / re))
    plt.plot(r_log, dgf, label=name, color=colors[name])

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Eccentricity (deg)')
plt.ylabel('Receptive Field Density (cells / deg$^2$)')
plt.title('Watson Model: RGCf Density vs Eccentricity')
plt.legend()
plt.grid(True, which="both", ls='--', lw=0.4)
plt.tight_layout()
plt.savefig("rgc_density_loglog.pdf")

# Curcio & Allen (1990) の実測データ処理（Drasdoモデルによる変換）
file_path = "Curcio_JCompNeurol1990_GCtopo_F6.xlsx"
df_raw = pd.read_excel(file_path, sheet_name=0, header=7)
df_curcio = df_raw.iloc[:, [0, 1, 3, 5, 7]].copy()
df_curcio.columns = ['Ecc_mm', 'Temp', 'Sup', 'Nasal', 'Inf']

# 偏心距離の補正式（Watson Appendix A6）
def ecc_mm_to_deg(r_mm):
    return (
        3.556 * r_mm +
        0.0599302 * r_mm**2 -
        0.00735803 * r_mm**3 +
        0.000302704 * r_mm**4
    )

# 面積補正係数 a(r₀)：mm²/deg² → Appendix A7 より（Watson 2014）
def area_conversion_factor(r_deg):
    return (
        0.0752 +
        5.846e-5 * r_deg -
        1.064e-5 * r_deg**2 +
        4.116e-8 * r_deg**3
    )

# mm → deg に変換（補正式に従って）
df_curcio['Ecc_mm'] = pd.to_numeric(df_curcio['Ecc_mm'], errors='coerce')
df_curcio = df_curcio.dropna(subset=['Ecc_mm'])
df_curcio['Ecc_deg'] = ecc_mm_to_deg(df_curcio['Ecc_mm']) 

# RGC密度（cells/mm²）→ cells/deg² に変換（面積補正式で）
for col in ['Temp', 'Sup', 'Nasal', 'Inf']:
    df_curcio[col] = pd.to_numeric(df_curcio[col], errors='coerce') / area_conversion_factor(df_curcio['Ecc_deg'])

df_plot = df_curcio[['Ecc_deg', 'Temp', 'Sup', 'Nasal', 'Inf']].dropna()


# Curcio 実測データ
plt.figure(figsize=(7, 5))
plt.plot(df_plot['Ecc_deg'], df_plot['Temp'], '--', label='Temporal (Curcio)', color='red')
plt.plot(df_plot['Ecc_deg'], df_plot['Sup'], '--', label='Superior (Curcio)', color='blue')
plt.plot(df_plot['Ecc_deg'], df_plot['Nasal'], '--', label='Nasal (Curcio)', color='green')
plt.plot(df_plot['Ecc_deg'], df_plot['Inf'], '--', label='Inferior (Curcio)', color='gray')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Eccentricity (deg)')
plt.ylabel('Measured RGC Density (cells / deg$^2$)')
plt.title('Curcio & Allen (1990): RGC Density vs Eccentricity')
plt.legend()
plt.grid(True, which="both", ls='--', lw=0.4)
plt.tight_layout()
plt.savefig("rgc_curcio.pdf")

# Watsonモデルによる累積RGC数
r_lin = np.linspace(0, 20, 500)
dr = r_lin[1] - r_lin[0]
plt.figure(figsize=(7, 5))
for name, p in params.items():
    a, r2, re = p['a'], p['r2'], p['re']
    dgf = dgf0 * (a * (1 + r_lin / r2)**-2 + (1 - a) * np.exp(-r_lin / re))
    cumulative = np.cumsum(dgf * 2 * np.pi * r_lin) * dr  # Watson式：累積は 2πr·密度
    plt.plot(r_lin, cumulative / 1000, label=name, color=colors[name])  # ×1000単位で表示

plt.xlabel('Eccentricity (deg)')
plt.ylabel('Cumulative RGC (×1000 cells)')
plt.title('Watson Model: Cumulative RGC Count')
plt.xlim(0, 20)
plt.ylim(0, 700)
plt.legend()
plt.grid(True, ls='--', lw=0.4)
plt.tight_layout()
plt.savefig("rgc_cumulative_linear.pdf")

# Watson と Curcio の密度重ね描き
plt.figure(figsize=(7, 5))

# Watsonモデル（実線）
for name, p in params.items():
    a, r2, re = p['a'], p['r2'], p['re']
    dgf = dgf0 * (a * (1 + r_log / r2)**-2 + (1 - a) * np.exp(-r_log / re))
    plt.plot(r_log, dgf, label=f'{name} (Watson)', color=colors[name], linewidth=2)

# Curcio実測データ（赤い点で表示）
plt.plot(df_plot['Ecc_deg'], df_plot['Temp'], 'o', label='Temporal (Curcio)', color='red')
plt.plot(df_plot['Ecc_deg'], df_plot['Sup'], 'o', label='Superior (Curcio)', color='blue')
plt.plot(df_plot['Ecc_deg'], df_plot['Nasal'], 'o', label='Nasal (Curcio)', color='green')
plt.plot(df_plot['Ecc_deg'], df_plot['Inf'], 'o', label='Inferior (Curcio)', color='gray')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Eccentricity (deg)')
plt.ylabel('RGC Density (cells / deg$^2$)')
plt.title('Watson Model vs Curcio Data: RGC Density Comparison')
plt.legend()
plt.xlim(0.01, 100)
plt.ylim(bottom=1)
plt.grid(True, which="both", ls='--', lw=0.4)
plt.tight_layout()
plt.savefig("rgc_comparison.pdf")

# 視野角度（deg）→ 網膜距離（mm）の変換とプロット

# Appendix式 (A5) に基づく変換式：r0mm(r0)
def deg_to_mm(r_deg):
    return (
        0.268 * r_deg +
        0.0003427 * r_deg**2 -
        8.3309e-6 * r_deg**3
    )

# グラフ作成範囲（0〜60度程度が妥当）
r_deg_range = np.linspace(0, 100, 500)
r_mm = deg_to_mm(r_deg_range)

# （参考）直線近似：1 deg ≈ 0.268 mm → 逆数
linear_approx = r_deg_range * 0.268

# プロット
plt.figure(figsize=(7, 5))
plt.plot(r_deg_range, r_mm, color='black')
plt.plot(r_deg_range, linear_approx, '--', label='Linear Approx (1 deg ≈ 0.268 mm)', color='red')
plt.xlabel('Eccentricity (deg)')
plt.ylabel('Retinal distance from visual axis (mm)')
plt.title('Appendix 6: deg to mm')
plt.grid(True, ls='--', lw=0.4)

plt.xlim(0, 100)
plt.ylim(0, 26.8)
plt.tight_layout()
plt.savefig("deg_to_mm.pdf")

# Appendix A6 の式：r0(mm) → r0(deg)
def mm_to_deg(r_mm):
    return (
        3.556 * r_mm +
        0.0599302 * r_mm**2 -
        0.00735803 * r_mm**3 +
        0.000302704 * r_mm**4
    )

# mmの範囲（通常は0〜6 mmで十分）
r_mm_range = np.linspace(0, 23, 500)
rdeg = mm_to_deg(r_mm_range)

# （参考）直線近似：1 deg ≈ 0.268 mm → 逆数
linear_approx = r_mm_range / 0.268

plt.figure(figsize=(7, 5))
plt.plot(r_mm_range, rdeg, color='black', lw=2)
plt.plot(r_mm_range, linear_approx, '--', label='Linear Approx (1 deg ≈ 0.268 mm)', color='red')

plt.xlabel('Retinal Distance (mm)')
plt.ylabel('Eccentricity (deg)')
plt.title('Appendix 6: mm to deg')
plt.grid(True, ls='--', lw=0.4)
plt.xlim(0, 23)
plt.ylim(0, 100)

plt.tight_layout()
plt.savefig("mm_to_deg.pdf")

# Figure A3 相当（視軸補正後の子午線ごとの変換)
plt.figure(figsize=(7, 5))
for name, D in offsets.items():
    r_mm_vis = r_mm - D
    rdeg_corrected = mm_to_deg(r_mm_vis) + mm_to_deg(D)
    plt.plot(r_mm, rdeg_corrected, label=name, color=colors[name])

plt.xlabel('Distance from Visual Axis (mm)')
plt.ylabel('Eccentricity (deg)')
plt.grid(True, ls='--', lw=0.4)
plt.xlim(0, 23)
plt.ylim(0, 100)
plt.legend()
plt.tight_layout()
plt.savefig("offset_mm_deg.pdf")

# Figure A4 相当（補正あり / なし の比較）
plt.figure(figsize=(7, 5))
for name, D in offsets.items():
    r_mm_vis = r_mm - D
    rdeg_no_corr = mm_to_deg(r_mm_vis)
    rdeg_corr = mm_to_deg(r_mm_vis) + mm_to_deg(D)
    plt.plot(rdeg_no_corr, rdeg_corr, label=name, color=colors[name])

plt.plot([0, 100], [0, 100], 'k--', label='No Correction = Correction')
plt.xlabel('Eccentricity without Offset Correction (deg)')
plt.ylabel('Eccentricity with Offset Correction (deg)')
# plt.title('Figure A4: Offset Correction Comparison')
plt.grid(True, ls='--', lw=0.4)
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.legend()
plt.tight_layout()
plt.savefig("offset_deg_mm.pdf")

# Appendix A7: 面積補正係数 a(r₀)（単位: mm²/deg²）
# 偏心度（視野角度）の範囲：0〜60 deg
r_deg_range = np.linspace(0, 100, 500)
area_factor = area_conversion_factor(r_deg_range)

plt.figure(figsize=(7, 5))
plt.plot(r_deg_range, area_factor, color='darkorange', lw=2)

plt.xlabel('Eccentricity (deg)')
plt.ylabel('Area Conversion Factor (mm² / deg²)')
plt.title('Appendix A7: Area ratio')
plt.grid(True, ls='--', lw=0.4)
plt.tight_layout()
plt.savefig("area_conversion_factor.pdf")