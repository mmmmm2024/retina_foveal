"""
NEURON網膜回路モデルを実行し、指定した細胞の膜電位トレースを保存するスクリプト。

処理:
- hocファイルを読み込み（parameters / createcells / netconnection）
- CoreNEURONを有効化してシミュレーションを高速実行
- 指定セルの膜電位を、t>=1000 ms から保存

出力:
- {cellname}_R{Num_R}_C{Num_C}.txt（ヘッダ: time(ms), voltage(mV)）

補足:
- plot_static() は刺激区間を色分けしてPDF保存
- animate_recording() は膜電位の動画(mp4)保存
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
    # if h.ENABLE_GRAPHICAL_INTERFACE:
    #     h.Plot()
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
    times, voltages = record_times_and_voltages()
    # segment into 1-2s, 2-3s, 3-4s
    x1, y1, x2, y2, x3, y3 = [], [], [], [], [], []
    for t, v in zip(times, voltages):
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

    # ax.set_xlim(1, 6)
    ax.set_xlim(1, 1.3)
    # ax.set_ylim(-50, -40)
    # ax.set_ylim(-80, -10)
    # plt.yticks([-100, -80, -60, -40, -20, 0, 20])
    ax.set_ylabel("Membrane Potential (mV)")
    # ax.set_title(f"{object_name}")

    # add stimulus rectangle
    ymin, ymax = ax.get_ylim()
    rect_h = (ymax - ymin) * 0.015
    y_pos = ymin + (ymax - ymin) * 0.02
    stim_rect = patches.Rectangle((2.0, y_pos), width=1.0, height=rect_h,
                                  color='black', zorder=10)
    ax.add_patch(stim_rect)

    # output_filename = f"{object_name}_RC{h.R_cov:.0f}_v_fovea.pdf"
    output_filename = f"{object_name}_R{h.Num_R:.0f}_C{h.Num_C_RP:.0f}.pdf"
    # output_filename = f"{object_name}_{h.Num_R:.0f}_vol_fovea.pdf"
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    # plt.show()

# === REFACTORED: Animation code ===
def animate_recording():
    from matplotlib.animation import FuncAnimation
    times, voltages = [], []

    fig, ax = plt.subplots(figsize=(10, 5))
    line, = ax.plot([], [], lw=1.5)
    ax.set_xlim(1, 4)
    ax.set_ylim(-75, -40)
    ax.set_ylabel("Membrane Potential (mV)")
    ax.set_title(f"{object_name}")

    def update(frame):
        step()
        if h.t >= 1000:
            t_s = h.t / 1000.0
            times.append(t_s)
            voltages.append(getattr(target_object.soma, name))
            line.set_data(times, voltages)
        return [line]
    # Add stimulus rectang

    frame_count = int(h.tstop / h.dt)
    anim = FuncAnimation(fig, update, frames=frame_count, interval=0, blit=True)
    anim.save(f"membrane_{object_name}.mp4", writer="ffmpeg", fps=100)
    plt.close(fig)

# === NEW FUNCTION: Export data to text file ===
def export_data_txt(filename):
    """
    Run a fresh simulation and write time (ms) and membrane potential (mV)
    for t >= 1000 ms into a CSV text file with headers.
    """
    with open(filename, 'w') as f:
        f.write("time(ms),voltage(mV)\n")
        while h.t < h.tstop:
            step()
            if h.t >= 1000:
                voltage = getattr(target_object.soma, name)
                f.write(f"{h.t:.4f},{voltage:.4f}\n")

    print(f"Data exported to {filename}")


def record_times_and_voltages():
    times, voltages = [], []
    while h.t < h.tstop:
        step()
        if h.t >= 1000:
            times.append(h.t / 1000.0)
            # voltages.append(h.AC2OFFCB[0][0].i)
            voltages.append(getattr(target_object.soma, name))
    return times, voltages


# run_timer_simple.py
from time import perf_counter
from datetime import datetime, timezone, timedelta

start_wall = datetime.now(timezone.utc)     # 開始時刻（UTC）
t0 = perf_counter()

init()
start(h.AMP, h.Num_R, h.Num_C)
# h.block_ih(h.gihbar)

#target_object = h.OFF_GC[0]
target_object = h.OFF_GC[0]
name = 'v'
object_name = object_names.get(target_object)

if __name__ == "__main__":
    # print(f"{object_name}_{h.Num_R:.0f}")
    # plot_static()
    export_data_txt(f"{object_name}_{h.Num_R:.0f}.txt")
    # animate_recording()
    print(f"{object_name}_{h.Num_R:.0f}")
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