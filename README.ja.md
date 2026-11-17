# FEAL-8_Re-Call-for-Attack

**Re: Call for Attack — 1989 年の FEAL-8 懸賞問題（予選問題）の鍵を求めた**

[English](README.md)

```
KEY = A5EA6EC2FB0CC741
```

1989 年 8 月、NTT は "The FEAL-8 Cryptosystem and a Call for Attack"（宮口庄司、CRYPTO '89）で FEAL-8 の懸賞問題を発表した。予選問題は、平文 8192 バイトと、それを秘密の 64 ビット鍵 1 本・ECB モードの FEAL-8 で暗号化した暗号文 8192 バイトから、鍵を求めるものである。応募期間は 1991 年 8 月 31 日に終わっている。2012 年に松井充氏はこの問題を未解決のままと述べており、解いたという公表は見当たらない。

このリポジトリでは、鍵、問題データ、鍵を求めた手法、ソースコード、探索のログを公開する。鍵が求まったのは 2026 年 11 月 17 日。

## 鍵の確認

```sh
python3 src/verify_key.py A5EA6EC2FB0CC741 data/plaintext.bin data/ciphertext.bin
# key a5ea6ec2fb0cc741: 1024/1024 blocks match
# block0 0909814081408265 -> 0b2ca544e0c77e30 file 0b2ca544e0c77e30
```

`verify_key.py` は単体で動く FEAL-8 の実装で、最初に標準のテストベクタ（鍵 0123456789ABCDEF、平文 0 → CEEF2C86F2490752）を確かめる。隣の鍵 A5EA6EC2FB0CC740 では 1 ブロックも一致しない。

## 求め方

分割（partitioning）を組み合わせた既知平文の線形解読である。Biham–Carmeli（SAC 2014）と先小山ら（IEICE 2016）の手法（どちらも FEAL-8X が対象）を組み合わせた。詳しくは [docs/attack.ja.md](docs/attack.ja.md)。

1. **データ**：1024 ブロック、異なる組は 981（既知平文 2^10 組。選択平文攻撃は使えない）
2. **段階 1**：EK0（15 ビット）、DK7（22 ビット）、DK6 の 8 ビットの仮説 2^45 個すべてについて、第 6 段の 8 ビット分割統計量を、2 つの族の 5 段線形近似で同時に計算し、z0² + z1² ≥ r1² のものを残す
3. **段階 2**：残り 15 ビットを補い、5 段近似（相関 1/4、正解鍵で各族約 8σ）の統計量を XOR と popcount で計算し、r ≥ 8.5 を候補とする
4. **段階 3**：近似を両方向に 1 段ずらして EK0、EK1、DK6、DK5 を求め、3 段近似で組を判定し（正解は約 16σ）、鍵スケジュールを逆算する
5. **両方向**：復号は暗号化と同じ構造なので、平文と暗号文を入れ替えても同じ攻撃が使える

鍵は、逆方向・r1 = 3.44 の探索で、RTX 5090 ×8 を使って開始から約 40 分で見つかった（片方向の全空間は約 1 時間）。同じ機材で総当たりすると平均約 35 日かかる。しきい値を厳しくしたそれまでの探索では、わずかの差で見逃していた。その理由は正解鍵の統計量から分かる（[results/experiments/actual_key_statistics.txt](results/experiments/actual_key_statistics.txt)）。

## 再現

```sh
make -C src && make -C src check        # ビルド（OpenMP 対応の gcc）と、鍵・データの確認
bash scripts/reproduce_quick.sh         # CPU で 5〜10 分。鍵を含んでいた DK7 で段階 1〜3 をやり直し、鍵を復元する

make -C src gpu ARCH=sm_120             # CUDA（RTX 30xx は sm_86、RTX 40xx は sm_89）
bash scripts/run_search.sh rev 3.44     # 全 GPU で片方向の全探索（RTX 5090 ×8 で約 1 時間）
```

## 構成

```
data/      懸賞問題の平文・暗号文（.bin、16 進、文章部の文字）、印刷された CBC-MAC、鍵
docs/      攻撃の説明（英語・日本語）、参考文献
src/       攻撃（段階 1〜3、CPU 版・CUDA 版）、鍵の検証、試験データ作成、補助実験
scripts/   全探索の実行スクリプトと短時間の再現スクリプト
results/   すべての探索のログ（鍵を出したものを含む）と実験結果
```

## データの出典

問題データは、松井甲子雄『コンピュータによる暗号解読法入門』（森北出版、1990 年、ISBN 4-627-83350-4）の付図6〜9 から手で書き起こし、同書に載っている CBC-MAC 64 個との照合で確かめた（すべて一致）。詳しくは [data/README.md](data/README.md)。

## ライセンス

このリポジトリのために書いたコードと文書は [Unlicense](LICENSE) でパブリックドメインとする。`data/` の問題データは、検証と研究のために上記の出版物から転載したもので、平文に含まれる文章は宮口・白石・清水（NTT）による FEAL-8 の論文の冒頭である。これらには Unlicense は適用されない。
