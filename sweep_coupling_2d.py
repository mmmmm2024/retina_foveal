#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import product
from pathlib import Path
import subprocess, sys, re, shutil

BASE = Path(__file__).resolve().parent
PY   = sys.executable

PARAM_PATH = BASE / "src" / "parameters_new.hoc"
INIT = str((BASE / "init.py").resolve())
TARGET = "OFF_GC"

# ===== sweep settings (ONLY these two) =====
x_list = [i / 20 for i in range(0, 21)]  # 0.0 ... 1.0
y_list = [i / 20 for i in range(0, 21)]  # 0.0 ... 1.0

# x_list = [4 / 20]  # 0.0 ... 1.0
# y_list = [11 / 20]  # 0.0 ... 1.0

# base values (fixed)
BASE_g_RBC2AC = 0.0012
BASE_gj_AC2CB = 0.0005

RESULTS_DIR = BASE / f"Dim_{TARGET}_8-9M_gRBC2ACx_gjAC2CBy"
# RESULTS_DIR = BASE / f"Dim_{TARGET}_8-9M_gRBC2ACx_gjAC2CBy"
RESULTS_DIR.mkdir(exist_ok=True)

ORIGINAL_TEXT = PARAM_PATH.read_text(encoding="utf-8")


def replace_var(src: str, name: str, value_str: str) -> str:
    """
    hocの「name = value ...」の value 部分だけ差し替える（コメント等は維持）。
    """
    pattern = re.compile(
        rf'^(\s*{re.escape(name)}\s*=\s*)([^/\n]*)(.*)$',
        re.MULTILINE
    )
    if not pattern.search(src):
        raise RuntimeError(f"[ERROR] {name} が {PARAM_PATH.name} 内に見つかりません。")

    return pattern.sub(lambda m: f"{m.group(1)}{value_str}{m.group(3)}", src, count=1)


def make_param_text(x: float, y: float) -> str:
    """
    g_RBC2AC と gj_AC2CB だけを書き換えたhocテキストを返す。
    それ以外は ORIGINAL_TEXT のまま固定。
    """
    text = ORIGINAL_TEXT
    text = replace_var(text, "g_RBC2AC", f"{BASE_g_RBC2AC:.6g} * {x:.2f}")
    text = replace_var(text, "gj_AC2CB", f"{BASE_gj_AC2CB:.6g} * {y:.2f}")
    return text


def write_parameters(x: float, y: float) -> None:
    new_text = make_param_text(x, y)

    # 念のためバックアップを1回だけ作る
    bak = PARAM_PATH.with_suffix(PARAM_PATH.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(PARAM_PATH, bak)

    PARAM_PATH.write_text(new_text, encoding="utf-8")


def restore_parameters() -> None:
    bak = PARAM_PATH.with_suffix(PARAM_PATH.suffix + ".bak")
    if bak.exists():
        shutil.copy2(bak, PARAM_PATH)


try:
    for x, y in product(x_list, y_list):
        ix, iy = int(round(x * 100)), int(round(y * 100))
        print(f"[RUN] x={x:.2f}, y={y:.2f}  "
              f"(g_RBC2AC={BASE_g_RBC2AC}*{x:.2f}, gj_AC2CB={BASE_gj_AC2CB}*{y:.2f})")

        # 1) hoc書き換え（2変数のみ）
        write_parameters(x, y)

        # 2) 実行前に「既にあるtxt」を記録して、実行後に増えた分だけ拾う
        before = {p.resolve() for p in BASE.glob("*.txt")}

        res = subprocess.run(
            [PY, INIT],
            cwd=BASE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        if res.returncode != 0:
            log_path = RESULTS_DIR / f"FAIL_x{ix:02d}_y{iy:02d}.log"
            log_path.write_text(res.stdout, encoding="utf-8")
            tail = "\n".join(res.stdout.splitlines()[-40:])
            print(f"[WARN] FAILED x{ix:02d}_y{iy:02d} (code={res.returncode}) "
                  f"log -> {log_path.name}\n--- log tail ---\n{tail}\n--- end ---")
            continue

        after = [p for p in BASE.glob("*.txt") if p.resolve() not in before]
        if not after:
            # 念のためフォールバック（最新txt）
            out_candidates = list(BASE.glob("*.txt"))
            if not out_candidates:
                print(f"[WARN] x{ix:02d}_y{iy:02d} で .txt 出力が見つかりませんでした。")
                continue
            out_txt = max(out_candidates, key=lambda p: p.stat().st_mtime)
        else:
            # 新規に増えたtxtのうち最新
            out_txt = max(after, key=lambda p: p.stat().st_mtime)

        # 3) 条件が分かる名前で保存（x,yは0..10の整数化）
        new_name = RESULTS_DIR / f"{TARGET}_x{ix:03d}_y{iy:03d}.txt"
        if new_name.exists():
            new_name.unlink()
        shutil.move(str(out_txt), str(new_name))
        print(f"[OK] -> {new_name.name}")

    print("[DONE] all sweeps finished.")

finally:
    # 実験後にparameters_new.hocを元に戻したい場合は有効
    restore_parameters()
