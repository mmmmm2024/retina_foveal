"""
NEURON網膜回路モデルを実行し、指定したシナプスの変数をCSV形式で保存するスクリプト。

処理:
- hocファイルを読み込み、CoreNEURONを有効化してシミュレーションを実行
- 刺激、ノイズ、リボンシナプス、抑制性シナプス、各種ギャップ結合を設定して初期化
- シナプス配列名（例: ONCB2ONGC）と pre/postを指定して計測対象を選択
- 指定変数（例: p1 / u / w / g / i など）を、t>=1000 ms から time(ms) と一緒に保存

入力:
- hoc: src/parameters_new.hoc, createcells.hoc, src/netconnection_fovea.hoc など

出力:
- {SYN_ARRAY_NAME}_{syn_name}_R{Num_R}_C{Num_C_RP}.txt
  例: ONCB2ONGC_p1_R400_C20.txt（ヘッダ: time(ms), p1）
"""

import neuron
from neuron import h
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
import matplotlib.patches as patches
from neuron import coreneuron

# --- Load NEURON hoc files ---
print("=========== init.py =============")
h.load_file("stdrun.hoc")
coreneuron.enable = True 
coreneuron.nthread = 8
coreneuron.cell_permute = 1 
h.cvode.cache_efficient(1)
h.load_file("src/parameters_new.hoc")
h.load_file("createcells.hoc")
h.load_file("src/netconnection_fovea.hoc")

# --- Simulation Control ---
def step():
    h.fadvance()

def init():
    h.finitialize()
    h.fcurrent()
    h.dt = h.step_dt
    h.tstop = h.tstop

def start(amp, R, C):
    # Call functions defined in hoc files
    h.iclamps(amp)       # in createcells.hoc
    h.noise()             # in createcells.hoc
    h.noise_set()         # in createcells.hoc
    h.Ribbon_syn(R, C)  # in netconnection.hoc
    h.Ribbon_syn_set()  # in netconnection.hoc
    h.Gly_syn(C)           # in netconnection.hoc
    h.Gly_syn_set()       # in netconnection.hoc
    h.Cone_GJ()
    h.Cone_GJ_set()
    h.R_C_GJ()             # in netconnection.hoc
    h.R_C_GJ_set()         # in netconnection.hoc
    h.AC_ONBC_GJ()             # in netconnection.hoc
    h.AC_ONBC_GJ_set()         # in netconnection.hoc
    h.AC_GJ()             # in netconnection.hoc
    h.AC_GJ_set()         # in netconnection.hoc
    h.OFFGC_GJ()          # in netconnection.hoc
    h.OFFGC_GJ_set()      # in netconnection.hoc
    h.ONCB_GJ()           # in netconnection.hoc
    h.ONCB_GJ_set()       # in netconnection.hoc
    h.OFFCB_GJ()         # in netconnection.hoc
    h.OFFCB_GJ_set()      # in netconnection.hoc
    init()

def get_object_name(obj):
    """
    NEURONオブジェクトの名前を動的に取得する関数。
    """
    obj_repr = str(obj)
    name = obj_repr.split(".")
    return name

# setting objects
object_names = {}
object_names[h.Cones[0]] = "Cones"
object_names[h.Rods[0]] = "Rods"
object_names[h.R_BC[0]] = "R_BC"
object_names[h.ON_CBC[0]] = "ON_CBC"
object_names[h.OFF_CBC[0]] = "OFF_CBC"
object_names[h.AIIAC[0]] = "AIIAC"
object_names[h.ON_GC[0]] = "ON_GC"
object_names[h.OFF_GC[0]] = "OFF_GC"


# === REFACTORED: Static plot code ===
def plot_static():
    times, values = record_times_and_voltages()  # 関数名はそのまま使う（中身は“シナプス値”になる）
    x1, y1, x2, y2, x3, y3 = [], [], [], [], [], []
    for t, v in zip(times, values):
        if 1.0 <= t < 2.0:
            x1.append(t); y1.append(v)
        elif 2.0 <= t <= 3.0:
            x2.append(t); y2.append(v)
        elif 3.0 < t <= 6.0:
            x3.append(t); y3.append(v)

    fig, ax = plt.subplots(figsize=(10, 5))
    line1, = ax.plot(x1, y1, color='#274A78', lw=1.5, zorder=1)
    line2, = ax.plot(x2, y2, color='#FF6F40', lw=1.5, zorder=1)
    line3, = ax.plot(x3, y3, color='#274A78', lw=1.5, zorder=1)

    ax.set_xlim(1, 1.3)
    ax.set_ylabel(f"{syn_name} (synapse variable)")  # ←ラベルだけ変更（方法は同じ）

    ymin, ymax = ax.get_ylim()
    rect_h = (ymax - ymin) * 0.015
    y_pos = ymin + (ymax - ymin) * 0.02
    stim_rect = patches.Rectangle((2.0, y_pos), width=1.0, height=rect_h,
                                  color='black', zorder=10)
    ax.add_patch(stim_rect)

    output_filename = f"{syn_label}_R{h.Num_R:.0f}_C{h.Num_C_RP:.0f}.pdf"
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')


# === NEW FUNCTION: Export synapse data to text file ===
def export_data_txt(filename):
    """
    Run a fresh simulation and write time (ms) and synapse variable
    for t >= 1000 ms into a CSV text file with headers.
    （計測ループや書き方は元のまま。参照先だけシナプスに変更）
    """
    with open(filename, 'w') as f:
        f.write(f"time(ms),{syn_name}\n")
        while h.t < h.tstop:
            step()
            if h.t >= 1000:
                value = getattr(target_synapse, syn_name)  # ★ここだけが本質的な変更
                f.write(f"{h.t:.4f},{value:.6f}\n")

    print(f"Data exported to {filename}")


def record_times_and_voltages():
    """
    関数名はそのまま（呼び出し側を崩さないため）。
    中身は “シナプス値” を返すように変更。
    """
    times, values = [], []
    while h.t < h.tstop:
        step()
        if h.t >= 1000:
            times.append(h.t / 1000.0)
            values.append(getattr(target_synapse, syn_name))  # ★ここだけ変更
    return times, values


# -------------------------------
# Synapse selector (ここを編集するだけで計測対象シナプスを変えられる)
# -------------------------------
def pick_synapse(syn_array_name: str, post_idx: int, pre_idx: int):
    """
    まず [post][pre] を試して、だめなら [pre][post] を試す
    """
    syn_arr = getattr(h, syn_array_name)
    try:
        return syn_arr[post_idx][pre_idx], "post_pre"
    except Exception:
        return syn_arr[pre_idx][post_idx], "pre_post"


# run_timer_simple.py
from time import perf_counter
from datetime import datetime, timezone, timedelta

start_wall = datetime.now(timezone.utc)     # 開始時刻（UTC）
t0 = perf_counter()

init()
start(h.AMP, h.Num_R, h.Num_C)

# もともとの target_object は残してOK（使わなくても動く）
target_object = h.OFF_GC[0]
name = 'v'
object_name = object_names.get(target_object)

# ★ここからが追加：計測対象シナプスの指定
SYN_ARRAY_NAME = "ONCB2ONGC"   # hoc 側の配列名に合わせる
POST_IDX = 0
PRE_IDX  = 0
target_synapse, syn_order = pick_synapse(SYN_ARRAY_NAME, POST_IDX, PRE_IDX)

syn_name  = "p1"              # 例: "isyn", "i", "g", "u", "P1", "w" など
syn_label = f"{SYN_ARRAY_NAME}_{syn_name}"

print(f"[INFO] recording synapse: {syn_label}")
# （任意）利用可能な変数候補を見たいとき：print([x for x in dir(target_synapse) if x.startswith('_ref_')])


if __name__ == "__main__":
    print(f"{syn_label}_R{h.Num_R:.0f}_C{h.Num_C_RP:.0f}.txt")
    export_data_txt(f"{syn_label}_R{h.Num_R:.0f}_C{h.Num_C_RP:.0f}.txt")
    print("complete")

elapsed = perf_counter() - t0
end_wall = datetime.now(timezone.utc)       # 終了時刻（UTC）

# 日本時間に変換（UTC+9）
JST = timezone(timedelta(hours=9))
start_jst = start_wall.astimezone(JST)
end_jst = end_wall.astimezone(JST)

minutes, seconds = divmod(int(elapsed), 60)

print(f"Start : {start_jst.isoformat()}")
print(f"End   : {end_jst.isoformat()}")
print(f"Elapsed: {minutes}分 {seconds}秒")
