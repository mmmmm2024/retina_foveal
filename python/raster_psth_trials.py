"""
複数試行の膜電位データからスパイクを検出
ラスタープロットとPSTH（ビン幅bin_size）を作成してPDF保存

処理:
- data/ON_GC/ON_GC[0]_v_{amp}_{trial}.txt (trial=1..n_trials) を読み込む
- start_time〜end_time の範囲で、電位が閾値を跨いだ時刻をスパイクとして検出
- Raster: (trial番号, spike_time) を点で描画 → data/raster_ONGC_{amp}.pdf
- PSTH : 全試行のスパイク時刻をヒストグラム化 → data/PSTH_ONGC_{amp}.pdf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 解析パラメータの設定
start_time = 1000    # 解析開始時刻 [ms]
end_time = 6000      # 解析終了時刻 [ms]
bin_size = 50        # PSTHのビン幅 [ms]
bins = np.arange(start_time, end_time + bin_size, bin_size)  # ヒストグラム用のビン

n_trials = 10        # 試行数（seed=1～10）

# 各試行ごとのスパイク時刻を格納するリスト
# raster_data: 各要素は (trial_number, spike_time) のタプル
raster_data = []
# 全試行のスパイク時刻をひとまとめに（PSTH用）
all_spike_times = []

# 閾値
threshold = 0.0
amp = 10

# 各試行のファイルを読み込み、スパイク検出
for trial in range(1, n_trials + 1):
    filename = f'data/ON_GC/ON_GC[0]_v_{int(amp)}_{trial}.txt'
    # CSVファイルを読み込み（ヘッダー行がある前提）
    df = pd.read_csv(filename)
    
    # 時間列と電位列の名称はファイル中のヘッダーに合わせる
    # ここでは "time(ms)" と "voltage(mV)" となっているものとする
    # 解析対象の時間領域に絞る
    df = df[(df['time(ms)'] >= start_time) & (df['time(ms)'] <= end_time)]
    
    # numpy配列に変換
    time_array = df['time(ms)'].values
    voltage_array = df['voltage(mV)'].values

    # スパイク検出（下から上への閾値通過を検出）
    # ここでは、あるサンプルで閾値以下で、次のサンプルで閾値以上になった場合をスパイクとする
    spike_times = []
    for i in range(1, len(voltage_array)):
        if voltage_array[i-1] < threshold and voltage_array[i] >= threshold:
            spike_time = time_array[i]
            spike_times.append(spike_time)
            all_spike_times.append(spike_time)
            raster_data.append((trial, spike_time))
    
    print(f"Trial {trial}: {len(spike_times)} spikes detected.")

# ----- ラスタープロットの作成 -----
plt.figure(figsize=(10, 6))
# 各スパイクについて、x軸を時刻、y軸を試行番号としてプロット
for trial, spike_time in raster_data:
    plt.plot(spike_time, trial, 'k.', markersize=5)

output_filename = f"data/raster_ONGC_{int(amp)}.pdf"

plt.xlabel("Time (ms)")
plt.ylabel("Trial")
plt.title("Raster Plot of ONGC")
plt.xlim(start_time, end_time)
plt.ylim(0.5, n_trials + 0.5)
plt.tight_layout()
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# ----- PSTH (Peristimulus Time Histogram) の作成 -----
# 全試行のスパイク時刻からヒストグラムを作成
counts, bin_edges = np.histogram(all_spike_times, bins=bins)
# 各ビンの中心位置を計算
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

output = f"data/PSTH_ONGC_{int(amp)}.pdf"

plt.figure(figsize=(10, 4))
plt.bar(bin_centers, counts, width=bin_size, align='center')
plt.xlabel("Time (ms)")
plt.ylabel("Spike Count")
plt.title("PSTH of ONGC")
plt.xlim(start_time, end_time)
plt.tight_layout()
plt.savefig(output, dpi=300, bbox_inches='tight')
