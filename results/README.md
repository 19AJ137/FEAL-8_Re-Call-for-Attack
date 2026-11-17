# Logs and results / ログと結果

Raw logs and scripts from the search servers and from the experiments, kept as they were produced
(only collected into this tree). Times are UTC unless stated otherwise. Scripts refer to the
server's working directory (`/root/feal`) and to binaries built from [`../src`](../src)
(`./sk/attackC2_gpu` = `attackC2_gpu`). Comments in the scripts are in Japanese.

サーバと実験の生ログ・スクリプトを、出力されたままの形で集めたもの。時刻は特記なければ UTC。スクリプト中のパス（`/root/feal`）や実行ファイル名（`./sk/attackC2_gpu` など）はサーバ上のもの。

## Timeline / 経過

| Time (UTC) | Event | Folder |
|---|---|---|
| 11-16 16:03 – 11-17 00:54 | Forward, one family (\|z_B\| ≥ 3.5), RTX 3070; 126k candidates to stage 3, no key | `earlier_runs/fwd1_3070` |
| 11-17 00:40 – 03:22 | Reverse, two families (r1 = 3.9175), RTX 3070, DK7 0–702,463; 14 candidates | `earlier_runs/rev2_3070` |
| 03:23 – 03:42 | Reverse, two families (r1 = 3.9175), rest of the space on 8 × RTX 5090; no key | `earlier_runs/rev2_5090` |
| 03:42 – 04:04 | Forward, two families (r1 = 3.9175), 8 × RTX 5090; no key | `earlier_runs/fwd2_5090` |
| 04:24 – 04:26 | All 138 two-family candidates through stage 3 with relaxed thresholds; no key | `earlier_runs/relaxed_stage3` |
| 04:27 | Reverse, two families (r1 = 3.44), 8 × RTX 5090 started (`run9.sh`) | `key_recovery` |
| **05:06:31** | **Key `a5ea6ec2fb0cc741` recovered (14:06:31 JST)** | `key_recovery` |
| 05:15 | Reverse at r1 = 3.0 on a second server, stopped about a minute later | `earlier_runs/r1_3.0_prep` |

## `key_recovery/` — the run that found the key / 鍵を出した探索

| File | Contents |
|---|---|
| `run9.sh` | Driver: stages 1–2 on 8 GPUs (r1 = 3.44, rpre = 2.486, r2 = 8.5), every 2 min new candidates to stage 3 (T1 = 1.0, T2 = 4.5, `REV=1`) |
| `rev3_start.txt` | Start time |
| `rev3_gpu_0..7.log` | Stage 1–2 progress of each GPU (2048 DK7 values per line). GPU 6 covered DK7 0x300000–0x37ffff |
| `rev3_cands_0..7.txt` | Stage-2 candidates (`t e kappa e8 r  r-stat z0 z1`); the correct one is `34fd3b 6f7a 6b aa 45 10.1657 7.4391 6.9283` in `rev3_cands_6.txt` |
| `rev3_s3.log` | Stage-3 log. In the batch of 05:04:39, candidate 5 passed the 3-round test with z3 = −16.25 and the key was recovered |
| `FOUND_KEY.txt` | Key, log name and time written by `run9.sh` |

The literal `\n` after `KEY a5ea6ec2fb0cc741` and after `recovering from cand 5:` comes from an
escaping mistake in `printf` that was fixed in [`../src/stage3.c`](../src/stage3.c) afterwards.
`run9.sh` skipped candidates already tested at r1 = 3.9175 (the list `seen_rev.txt` is not included).

`KEY a5ea6ec2fb0cc741` や `recovering from cand 5:` の後ろの `\n` は printf のエスケープの誤りによるもので、公開版の `stage3.c` では直してある。

## `earlier_runs/` — searches that did not find the key / 鍵が出なかった探索

| Folder | Contents |
|---|---|
| `fwd1_3070` | `run_real.sh`, `run_chain.sh`, GPU and stage-3 logs, and `tier1_all.txt.gz` (the 125,847 candidates with \|z\| ≥ 7.3 that went to stage 3). The full candidate list (19.6 M lines, 167 MB compressed) is not included |
| `rev2_3070` | `run_chain2.sh`, logs and the 14 candidates found on the RTX 3070 |
| `rev2_5090` | `run8.sh` (reverse then forward, r1 = 3.9175), the 14 earlier candidates re-tested (`old_rev2_*`), logs and candidates (`rev2b_*`), all 73 reverse candidates (`all_rev2.txt`) |
| `fwd2_5090` | Forward part of `run8.sh`: logs, candidates, all 65 (`all_fwd2.txt`) |
| `relaxed_stage3` | Stage 3 with T1 = 1.0, T2 = 4.5 on the 138 candidates. `relax2_*.log` (04:24–04:26) is the valid run; `relax_*.log` (04:06) used a stage 3 that could drop the correct hypothesis when more than 4096 passed (fixed at 04:20). `f_fwd.txt`, `f_rev.txt`: stage-2 candidates of known-key test data (seed 11, key 39287fc26939a7df) in the forward and reverse directions |
| `r1_3.0_prep` | Preparation of r1 = 3.0 on a second server: `run10.sh`, known-key tests in `t30/` and `cap_test/` (build with a larger survivor cap), `chain10.*` |

## `experiments/` — measurements quoted in the documentation / 文書で引用した測定

| Folder / file | Contents |
|---|---|
| `actual_key_statistics.txt` | Stage 1–2 statistics of the recovered key in both directions, and why the earlier runs missed it |
| `feasibility/results.txt` | Output of `corr` and `null`: statistics of correct and wrong keys with the challenge plaintexts vs. random plaintexts |
| `two_family/` | Outputs of `bias`, `mlc`, `refstat` (30,000 keys, gzip), `null2`, `analyze_ref.py`, `cmp_one.sh` (60 known keys) |
| `stage12_tests/` | One-family stages 1–2 on known keys (`attackC self`, thresholds 3.0 and 3.5, 2400 DK7 each) |
| `stage3_tests/` | Stage 3 on 240 known keys with default (`known_keys_default_*`, 215 recovered) and relaxed thresholds (`known_keys_relaxed_*`, 239 recovered); wrong-candidate tests (`s3_wrong.log`, `s3_topwrong.log`, `s3_truecl.log`) |
| `reverse_mapping_check/` | Check that the attack also works if the data were C = D_K(P) (test data from `mkdt`) |
| `bruteforce/` | `bs_gpu` benchmarks: 7.6 × 10^10 keys/s on RTX 3070, 3.8 × 10^11 keys/s on RTX 5090 |
