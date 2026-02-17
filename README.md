# retina_foveal
# 👁️ Human Retina Model (Eccentricity: 1mm & 12mm)
### ヒト網膜数理モデル：生理学的シミュレーション・プラットフォーム

![Retina Layer Diagram](https://img.shields.io/badge/Platform-NEURON%20%7C%20Python-blue)
![Language](https://img.shields.io/badge/Language-C%20%7C%20Python-orange)

> **"Code is read much more often than it is written."** > 本プロジェクトは、ヒト網膜という「精密な並列信号処理回路」を数理的に再現したものです。C言語を主戦場とするエンジニアが、Pythonという制御レイヤーを通じて生物物理学的現象を解析できるように設計されています。

---

## 1. 概要 (Overview)

本プロジェクトは、ヒト網膜の **中心窩（偏心度 1mm）** および **周辺部（偏心度 12mm）** における神経回路を再現した数理モデルです。

* **エンジン**: **NEURON (C/C++ベース)**。物理的なイオンチャネルや膜電位の計算を担当。
* **コントローラー**: **Python**。実験条件の制御、並列実行、データ解析を担当。
* **ターゲット**: 網膜生理学やPythonの経験がないエンジニア（C言語のポインタや構造体を理解していることを前提とします）。

### コア機能
1.  **連続実験（Parameter Sweep）**: 条件を動的に変更し、数百〜数千パターンの実験を自動実行。
2.  **変性シミュレーション**: 「一次変性（光受容体の消失）」および「二次変性（ネットワークの崩壊）」のプロセスを追跡。
3.  **計算機的検証**: 生理学データ（細胞密度・接続数）が、計算モデル上で正しく再現されているかを自動検証。

---

## 2. クイックスタート (Quickstart)

Cエンジニアにとっての「ビルド（Build）」と「実行（Execute）」のプロセスです。

### 2.1 依存ライブラリのインストール
```bash
# 仮想環境の作成と有効化
python3 -m venv .venv
source .venv/bin/activate
```
# 必要なパッケージのインストール
pip install --upgrade pip
pip install numpy matplotlib pandas neuron==8.2.4

2.2 モデルのコンパイル (Compilation)

NEURONの低レイヤーコード（.mod ファイル）をコンパイルします。
C言語の make に相当する重要な作業です。

nrnivmodl mod


[!IMPORTANT]
実行後、arm64/ または x86_64/ フォルダが生成されていれば成功です。
mod/ 内を書き換えた際は、必ず再実行（再コンパイル）してください。

2.3 動作確認 (Smoke Run)

最小構成でシミュレーションを実行し、環境が正常か確認します。

```python init.py

3. プロジェクト構造 (Project Structure)

コードの役割を「回路の階層」として理解するためのマップです。

ディレクトリ/ファイル	役割 (Role)
cell/	モデル本体。細胞の形状や接続の定義。解析の核心部。
mod/	計算エンジン。NMODL形式で書かれたCコード群（チャネル・シナプス）。
src/	ユーティリティ。Python側のヘルパー関数やラッパー。
combine/	データ統合。複数条件の集計や、偏心度を跨ぐデータ処理用。
results/	出力先。シミュレーション結果（CSV, PNG等）が保存される。
*.hoc	NEURON独自の言語で書かれた細胞構築定義ファイル。
4. 実験の実行方法 (How to Run Experiments)

用途に合わせて以下のスクリプトを使い分けてください。

生存率スイープ

python sweep_survival.py


変性進行（細胞死）に伴う回路の変化をシミュレートします。

結合強度解析

python sweep_coupling.py


Gap Junctionなどの接続条件を振り、応答特性の変化を解析します。

可視化（ヒートマップ）

python aiiac_heatmap_rc_peripheral.py


解析結果を視覚的なマップとして出力し、傾向を把握します。

5. 設計哲学と10年後の後輩へ (Philosophy)

このシステムは、単一の「巨大な網膜」をシミュレートするものではありません。

なぜ一体化しないのか:
網膜の物理的な接続特性は偏心度（中心からの距離）によって劇的に変化します。
これらを一つの関数にまとめると、パラメータの混濁を招き、バグの温床となります。

拡張のルール:
新しい偏心度（例：5mm）を追加する場合は、既存の 1mm や 12mm の構成を「テンプレート」として複製・修正してください。
これにより、過去の実験との比較可能性（再現性）が担保されます。

6. トラブルシューティング (Troubleshooting)

ImportError: No module named 'neuron'
仮想環境が有効化されていますか？ source .venv/bin/activate を確認してください。

mechanism not found エラー
nrnivmodl mod を実行しましたか？ 実行ディレクトリがリポジトリルートであることを確認してください。

シミュレーション結果が更新されない
mod/ 内のCコードを修正しましたか？ 再コンパイルしない限り、変更は反映されません。
