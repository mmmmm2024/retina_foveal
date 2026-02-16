import pandas as pd
import numpy as np

# ファイル名
filename = "ON_GC_200.txt"

# ヒステリシスしきい値 (mV)
thr_hi = 0.0     # 上抜けでスパイク開始
thr_lo = -20.0   # ここまで下降したら次の検出を許可

# 時間範囲 (ms)
start_time = 1000.0
end_time   = 6000.0

# ファイルを読み込み（元コードに合わせて1行スキップ）
# ※ header=None と names を指定して、先頭データ行がヘッダー扱いされないようにする
data = pd.read_csv(filename, skiprows=1, header=None, names=["time", "voltage"])

# 数値化（非数は NaN → 除去）で型比較エラーを防止
data["time"] = pd.to_numeric(data["time"], errors="coerce")
data["voltage"] = pd.to_numeric(data["voltage"], errors="coerce")
data = data.dropna(subset=["time", "voltage"])

# 指定時間範囲で切り出し（ここは元コードの意図を踏襲）
seg = data[(data["time"] >= start_time) & (data["time"] <= end_time)].reset_index(drop=True)

# ===== ヒステリシスによるスパイク検出 =====
t = seg["time"].to_numpy()
v = seg["voltage"].to_numpy()

spike_times = []
armed = True  # True=次のスパイク検出が可能
for i in range(1, len(v)):
    if armed:
        # 上向き閾値を超えたらスパイク
        if v[i-1] <= thr_hi and v[i] > thr_hi:
            spike_times.append(t[i])  # 必要なら線形補間に変更可
            armed = False
    else:
        # 下向き閾値を下回ったら次の検出を許可
        if v[i] < thr_lo:
            armed = True

# 以降は元の集計フローに合わせる
spike_times = pd.Series(spike_times)
valid_spikes = spike_times  # 時間範囲内で既に抽出済み

# スパイク数と発火率を計算
spike_count = len(valid_spikes)
duration_s  = (end_time - start_time) / 1000.0  # 観測時間 (秒)
firing_rate = spike_count / duration_s if duration_s > 0 else float("nan")  # 発火率 (Hz)

# 結果を表示
print(filename)
print(f"Spike count (1000ms-6000ms): {spike_count}")
print(f"Firing rate (1000ms-6000ms): {firing_rate:.2f} Hz")
