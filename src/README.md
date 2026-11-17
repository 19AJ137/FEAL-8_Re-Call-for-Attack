# Source code / ソースコード

The attack is described in [`../docs/attack.md`](../docs/attack.md) ([日本語](../docs/attack.ja.md)).
Comments in the C/CUDA sources are in Japanese.
C／CUDA のソース中のコメントは日本語。

## Build / ビルド

```sh
make             # CPU programs (gcc with OpenMP)
make gpu         # CUDA programs; ARCH=sm_86 (RTX 30xx, default), sm_89 (RTX 40xx), sm_120 (RTX 50xx)
make check       # verify the key and the CBC-MAC of the data
```

Compiled with gcc 13 and CUDA 13.5. The GPU programs were run on RTX 3070 and RTX 5090 with CUDA 13.0.

## The attack / 攻撃本体

| File | Purpose |
|---|---|
| `feal8.h` | FEAL-8, equivalent representation, approximation masks, partition tables, data loading |
| `approx8.h` | The eight 6-round approximations of Sakikoyama et al. (two families); stages 1–2 use their 5-round parts |
| `attackC2.c` | Stages 1–2 on the CPU (OpenMP) |
| `attackC2_gpu.cu` | Stages 1–2 on the GPU; output identical to `attackC2.c` |
| `stage3.c` | Stage 3: remaining subkeys, 3-round test and key-schedule inversion → 64-bit key |

```sh
# stages 1-2: DK7 indices [start, start+count); the GPU version appends to out.txt, the CPU version overwrites it
./attackC2_gpu real P.bin C.bin start count out.txt [r1 rpre r2]     # defaults 3.9175 2.486 8.5
./attackC2     real P.bin C.bin start count out.txt [r1 rpre r2]     # defaults 3.9175 2.486 6.0
# stage 3: reads the first six columns of the candidate file
RECOVER=1 [REV=1] [T1=1.0 T2=4.5] ./stage3 real P.bin C.bin out.txt
```

- Forward direction: `P.bin` = `../data/plaintext.bin`, `C.bin` = `../data/ciphertext.bin`.
- Reverse direction: swap the two files and set `REV=1` for `stage3`. `stage3` treats any value of
  `REV` as reverse, so leave it unset for the forward direction.
- The candidate file has the columns `t e kappa e8 r  r-stat z0 z1`: DK7 index (22 bits), EK0 (15 bits),
  κ, the remaining DK7 and DK6 bits, and the stage-2 statistics.
- `stage3` options (environment variables): `T1`, `T2` (thresholds of steps (a) and (c); default 2.0 and
  5.0), `T3` (3-round test; default 8.0), `RECOVER=1` (run step (e) and print `KEY ...`), `REV=1`,
  `ONEMAP=1` (do not also try the inverse mapping in step (e)), `LIMIT=n` (first n candidates only).
- `attackC2_gpu` keeps at most `CAP` stage-1 survivors per DK7 (default 65,536, enough for r1 ≥ 3.44;
  build with `-DCAP=262144` for r1 = 3.0).
- Self-test mode with a random key derived from a seed:
  `./attackC2 self ../data/plaintext.bin 11 64` (true DK7 + 63 random DK7) and
  `./stage3 self ../data/plaintext.bin 11` (true candidate only). `mkct_seed` writes the ciphertext for
  the same key so that the `real` path can be tested end to end.

The scripts in [`../scripts`](../scripts) run these programs for a full search and for a quick
reproduction of the successful step.

## Utilities / 補助

| File | Purpose |
|---|---|
| `verify_key.py` | Independent Python implementation; checks a key against all 1024 blocks |
| `feal8.py` | Python FEAL-8 and CBC-MAC (used to check the transcription) |
| `check_mac.py` | Checks `data/*.bin` against the printed CBC-MAC values |
| `mkct.c` | `./mkct KEY P.bin out.bin` — ECB-encrypt a plaintext file under a given key |
| `mkct_seed.c` | `./mkct_seed P.bin seed out.bin` — same with the key used by the `self` modes |
| `mkdt.c` | `./mkdt P.bin seed out.bin` — "decrypt" instead (test data for C = D_K(P)) |
| `bs_gpu.cu` | Bitsliced brute force on the GPU (`./bs_gpu test`, `./bs_gpu bench`) for comparison |

## Supporting experiments / 補助実験（`analysis/`）

These programs produced the measurements quoted in the documentation. Their outputs are in
[`../results/experiments`](../results/experiments).

| File | Purpose |
|---|---|
| `attackC.c`, `attackC_gpu.cu` | Earlier one-family version of stages 1–2 (used for the first forward search) |
| `corr.c` | Statistics of the correct key (method A, one 6-round approximation, 8-bit partitioning); test vector and equivalent-representation checks |
| `null.c` | Distribution of the statistics for wrong keys |
| `refstat.c` | Stage 1–2 statistics of the correct key for both families and directions (`KEY=... ./refstat P.bin 1 0` for a given key) |
| `null2.c` | Independence of the two families for wrong hypotheses |
| `bias.c` | Correlations of the eight 6-round approximations |
| `mlc.c` | Reproduction of the first stage of Sakikoyama et al. (rank of the correct key among 2^37) |
| `analyze_ref.py`, `analyze_mlc.py` | Summaries of the `refstat` and `mlc` outputs (numpy, scipy) |
| `cmp_one.sh` | Compares `attackC` and `attackC2` on one known key |

Build them with `make` (they are part of the default target); run them from `src/`.
