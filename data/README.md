# Challenge data / 懸賞問題のデータ

| File | Contents |
|---|---|
| `plaintext.bin` | Plaintext, 8192 bytes = 1024 blocks (付図6) |
| `ciphertext.bin` | Ciphertext, 8192 bytes = 1024 blocks (付図7) |
| `key.txt` | Recovered key `A5EA6EC2FB0CC741` |
| `plaintext.hex`, `ciphertext.hex` | One block per line: block number, byte offset, 64-bit block in hex |
| `plaintext_text.txt` | Bytes 0x0000–0x17FE of the plaintext decoded from Shift-JIS (CP932) to UTF-8, CR+LF → LF |
| `mac_fig8.txt`, `mac_fig9.txt` | CBC-MAC values printed with the problem for the plaintext (付図8) and the ciphertext (付図9) |
| `sha256sum.txt` | SHA-256 of the files above (`sha256sum -c sha256sum.txt`) |

The ciphertext is FEAL-8 in ECB mode: `C_i = E_K(P_i)` for every block, with K = `A5EA6EC2FB0CC741`.
The first block is `0909814081408265` → `0B2CA544E0C77E30`.

暗号文は鍵 `A5EA6EC2FB0CC741` による FEAL-8・ECB で、全ブロックについて `C_i = E_K(P_i)` が成り立つ。第 1 ブロックは `0909814081408265` → `0B2CA544E0C77E30`。

```sh
cd src
python3 verify_key.py A5EA6EC2FB0CC741 ../data/plaintext.bin ../data/ciphertext.bin   # 1024/1024 blocks match
python3 check_mac.py ../data                                                          # 32/32 and 32/32
```

## Source / 出典

The data were transcribed by hand from photographs of the appendix figures 6–9 (around pp. 197–206) of
松井甲子雄『コンピュータによる暗号解読法入門』森北出版, 1990, ISBN 4-627-83350-4 (2nd printing, 1997).
The book prints the plaintext and ciphertext as hex dumps together with CBC-MAC values for checking.

松井甲子雄『コンピュータによる暗号解読法入門』（森北出版、1990 年、ISBN 4-627-83350-4。第 1 版第 2 刷 1997 年）の付図6〜9（pp. 197–206 付近）を撮影し、手で書き起こした。同書には平文・暗号文の 16 進ダンプと、照合用の CBC-MAC が載っている。

## Checking the transcription / 書き起こしの確認

The printed MAC values follow this rule (all 64 values match, so the rule and the transcription are
both confirmed):

- K = 0, IV = 0, x_i = FEAL-8_K(D_i ⊕ x_{i−1}), chained over all 1024 blocks without restarting.
- After every 32 blocks (x_32, x_64, …, x_1024), the left 32 bits of x are printed: 32 values for the
  plaintext (付図8) and 32 for the ciphertext (付図9).

印刷された MAC は次の規則で計算されている（64 個すべてが一致したので、規則と書き起こしの両方が確認できた）。

- K = 0、IV = 0、x_i = FEAL-8_K(D_i ⊕ x_{i−1}) を 1024 ブロック通しで連鎖させる（32 ブロックごとに IV に戻さない）。
- 32 ブロックごと（x_32, x_64, …, x_1024）に x の左 32 ビットを書き出す。平文（付図8）と暗号文（付図9）で 32 個ずつ。

Four reading errors in the plaintext were found with the MAC and corrected (MAC groups #15, #22, #25
and #31). The ciphertext matched on the first attempt. Typographical
errors of the original text are kept as printed (e.g. "mudulus", "streucuter", "microproccessors",
"Shanon", 「混ぜ合わる」).

平文では MAC の 4 グループ（#15、#22、#25、#31）で読み違いが見つかり、MAC が合うように直した。暗号文は最初から全部一致した。原文の誤植（mudulus、streucuter、microproccessors、Shanon、「混ぜ合わる」など）は印刷どおりに残している。

## Structure / データの構造

- Plaintext 0x0000–0x17FE: the beginning of the FEAL-8 paper by Miyaguchi, Shiraishi and Shimizu
  (Japanese, then English), Shift-JIS with CR+LF line ends. 0x17FF is 0x1A (MS-DOS end of file).
- Plaintext 0x1800–0x1FFF: 2048 random-looking bytes.
- 981 distinct blocks (725 in the text part, 256 in the random part). 31 plaintext blocks occur more
  than once (74 positions in total); at all of these positions the ciphertext blocks repeat in the same way, and no two different plaintext
  blocks share a ciphertext block, as expected for ECB. The most frequent block, `838b8353838a8359`
  (Shift-JIS katakana), occurs 7 times.

- 平文 0x0000–0x17FE：宮口・白石・清水の FEAL-8 論文の冒頭（日本語、続いて英語）。Shift-JIS、改行は CR+LF。0x17FF は 0x1A（MS-DOS の EOF）。
- 平文 0x1800–0x1FFF：乱数状の 2048 バイト。
- 異なるブロックは 981 個（文章部 725、乱数部 256）。2 回以上現れる平文ブロックは 31 種類（計 74 か所）あり、そのすべてで暗号文も同じように重複している。異なる平文ブロックが同じ暗号文になることはなく、ECB と整合する。最も多いのは `838b8353838a8359`（Shift-JIS の片仮名）の 7 回。
