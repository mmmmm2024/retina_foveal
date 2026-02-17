# ヒト網膜数理モデル（偏心度 1 mm）

## プロジェクトの使命（Mission Statement）
中心窩近傍（偏心度 1 mm）の**ヒト網膜回路**を,マウス網膜モデルをベースに数理モデル化し,**網膜色素変性症（RP）**に伴う回路変化を,**神経節細胞（GC）**および **AIIアマクリン細胞（AIIAC）**の膜電位応答から解析できるようにする.

---

## システムアーキテクチャ（System Architecture）

### データフロー（入力 → 出力）
- 網膜は **視細胞（Rod/Cone）→（水平細胞）→ 双極細胞 → AIIアマクリン細胞 → 神経節細胞** からなる層状回路で構成される.  
- **視細胞に photocurrent を与える**ことで,光刺激に対する膜電位変化を再現する.  
- 細胞間結合は以下で表現する.  
  - **化学シナプス**：神経伝達物質放出により後細胞へ影響  
  - **興奮性/抑制性シナプス**：リボンシナプスモデル（4-state）および 3-state モデル  
  - **ギャップ結合**：ギャップ結合モデル  

---

## ビルドと展開（Deployment & Build）
依存関係を解決し,NEURON 実行環境を整えた上で **MOD（NMODL）ファイルをコンパイル**してから実行する.

```bash
python3 -m venv .venv && source .venv/bin/activate && pip
install --upgrade pip && pip install numpy matplotlib pandas
neuron==8.2.4 && git clone
https://github.com/mmmmm2024/retina02.git && cd retina
&& nrnivmodl mod
```

---

## ディレクトリ構造（Directory Structure）
- `cell/` : モデル本体（細胞形状・接続定義などの中核）
- `mod/` : 計算エンジン（NMODLで書かれたチャネル・シナプスなどのCコード群）
- `src/` : 実験パラメータ（条件設定・パラメータ管理）
- `python/` : 解析スクリプト（ネットワーク構造や応答の可視化・定量）

---

## 実験シナリオの実行（Experimental Scenarios）
- `init.py`：NEURON網膜回路モデルを実行し,指定した細胞の膜電位トレースをCSVとして保存する.  
- `init_syn.py`：NEURON網膜回路モデルを実行し,指定したシナプス変数（例：p1/u/w/g/i など）を時系列CSVとして保存する.  
- `sweep_rodcone.py`：Rod数×Cone数の2次元グリッドでシミュレーションを一括実行し,結果を条件ごとに整理して保存する.  
- `sweep_syn_rodcone.py`：Rod×Cone条件をスイープしつつ,指定シナプスの変数（例：isyn など）を保存する実験を一括実行する.  
- `sweep_coupling_2d.py`：回路の結合パラメータを2次元でスイープしてシミュレーションを実行し,条件ごとの出力を収集する.  
