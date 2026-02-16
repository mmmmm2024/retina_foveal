#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import product
from pathlib import Path
import subprocess, sys, re, shutil

BASE = Path(__file__).resolve().parent
PY   = sys.executable

PARAM_PATH = BASE / "src" / "parameters_new.hoc"

INIT = str((BASE / "init.py").resolve())
TARGET = "AIIAC"

Num_R_list    = [1, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400]
Num_C_RP_list = list(range(20, -1, -2))   # 20,18,...,0

# Num_R_list    = [1]
# Num_C_RP_list = [0]  

# 結果を出力したテキストをまとめるフォルダ
RESULTS_DIR = BASE / f"Dim_{TARGET}_fovea"
RESULTS_DIR.mkdir(exist_ok=True)

# --- parameters.hoc をいじるための準備 ---
# 元の parameters.hoc をテンプレートとして読み込んでおく

ORIGINAL_TEXT = PARAM_PATH.read_text(encoding="utf-8")

def make_param_text(num_r: int, num_c_rp: int, g_r2rb: float) -> str:
    """
    Num_R / Num_C_RP / g_R2RB だけ値を差し替えた文字列を返す。
    他の部分（コメントなど）は一切変えない。
    """
    text = ORIGINAL_TEXT  # 作業用

    def replace_var(src: str, name: str, value_str: str) -> str:
        # 例: Num_R    = 400     // コメント
        pattern = re.compile(
            rf'^(\s*{name}\s*=\s*)([^/\n]*)(.*)$',
            re.MULTILINE
        )

        if not pattern.search(src):
            print(f"[WARN] {name} が {PARAM_PATH.name} 内に見つかりませんでした。")
            return src  # 見つからなければ何も変えない

        def repl(m: re.Match) -> str:
            # m.group(1): "Num_R =" など左側
            # m.group(2): もとの数値部
            # m.group(3): コメントなど行末
            return f"{m.group(1)}{value_str}{m.group(3)}"

        return pattern.sub(repl, src, count=1)

    # 3つの変数だけ順番に書き換え
    text = replace_var(text, "Num_R",    str(num_r))
    text = replace_var(text, "Num_C_RP", str(num_c_rp))
    text = replace_var(text, "g_R2RB",   f"{g_r2rb:.8g}")

    return text  # ★ここが超重要：必ず文字列を返す


def write_parameters(num_r: int, num_c_rp: int, g_r2rb: float) -> None:
    new_text = make_param_text(num_r, num_c_rp, g_r2rb)
    if new_text is None:
        raise RuntimeError("make_param_text が None を返しました。")

    # 念のためバックアップを1回だけ作っておく
    bak = PARAM_PATH.with_suffix(PARAM_PATH.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(PARAM_PATH, bak)

    PARAM_PATH.write_text(new_text, encoding="utf-8")

# --- 実行ループ ---

for R, C in product(Num_R_list, Num_C_RP_list):
    # 条件ごとの g_R2RB
    g = 0.0 if R == 1 else 1e-5

    print(f"[RUN] Num_R={R:3d}, Num_C_RP={C:2d}, g_R2RB={g}")

    # 1) parameters.hoc の変数を上書き
    write_parameters(R, C, g)

    # 2) init.py を実行（毎回メモリ解放）
    res = subprocess.run([PY, INIT],
                         cwd=BASE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         text=True)

    if res.returncode != 0:
        # エラー時はログ末尾だけ表示して次へ
        tail = "\n".join(res.stdout.splitlines()[-20:])
        print(f"[WARN] FAILED R{R:03d}_C{C:02d} (code={res.returncode})\n--- log tail ---\n{tail}\n--- end ---")
        continue

    # init.py が出力した .txtをフォルダに移す（名前に R/C を付ける）
    out_candidates = list(BASE.glob("*.txt"))
    if not out_candidates:
        print(f"[WARN] R{R}_C{C} で .txt 出力が見つかりませんでした。")
        continue

    out_txt = max(out_candidates, key=lambda p: p.stat().st_mtime)

    # 出力ファイル名を AIIAC_Rxxx_Cyy.txt に固定
    # new_name = RESULTS_DIR / f"AIIAC_R{R:03d}_C{C:02d}.txt"
    new_name = RESULTS_DIR / f"{TARGET}_R{R:03d}_C{C:02d}.txt"

    if new_name.exists():
        new_name.unlink()
    shutil.move(str(out_txt), str(new_name))
    print(f"[OK] -> {new_name.name}")

print("[DONE] all sweeps finished.")
