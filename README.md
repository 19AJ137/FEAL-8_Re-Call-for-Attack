# FEAL-8_Re-Call-for-Attack

**Re: Call for Attack — the preliminary problem of the 1989 FEAL-8 prize challenge, solved**

[日本語](README.ja.md)

```
KEY = A5EA6EC2FB0CC741
```

In August 1989 NTT announced "The FEAL-8 Cryptosystem and a Call for Attack" (S. Miyaguchi,
CRYPTO '89). The preliminary problem gives 8192 bytes of plaintext and the corresponding 8192 bytes
of ciphertext, encrypted with FEAL-8 in ECB mode under one secret 64-bit key, and asks for the key.
The call closed on 31 August 1991. In 2012 M. Matsui described the problem as having remained
unsolved, and we are not aware of any published solution.

This repository contains the key, the challenge data, the attack that found it, the source code,
and the logs of the search. The key was recovered on 17 November 2026.

## Check the key

```sh
python3 src/verify_key.py A5EA6EC2FB0CC741 data/plaintext.bin data/ciphertext.bin
# key a5ea6ec2fb0cc741: 1024/1024 blocks match
# block0 0909814081408265 -> 0b2ca544e0c77e30 file 0b2ca544e0c77e30
```

`verify_key.py` is a self-contained FEAL-8 implementation that first checks the standard test vector
(key 0123456789ABCDEF, plaintext 0 → CEEF2C86F2490752). The neighbouring key A5EA6EC2FB0CC740 matches
no block.

## How it was found

Known-plaintext linear cryptanalysis with partitioning, combining Biham–Carmeli (SAC 2014) and
Sakikoyama et al. (IEICE 2016), both originally aimed at FEAL-8X. Full description:
[docs/attack.md](docs/attack.md).

1. **Data**: 1024 blocks, 981 distinct pairs (2^10 known plaintexts; chosen-plaintext attacks do not apply).
2. **Stage 1** — for all 2^45 hypotheses on EK0 (15 bits), DK7 (22 bits) and 8 bits of DK6, an
   8-bit partitioning statistic of the round-6 approximation, for two families of 5-round linear
   approximations at once; keep z0² + z1² ≥ r1².
3. **Stage 2** — complete the remaining 15 bits and compute the full 5-round statistic
   (correlation 1/4, about 8σ per family for the correct key) with XOR and popcount; keep r ≥ 8.5.
4. **Stage 3** — shift the approximation one round in each direction to get EK0, EK1, DK6, DK5, test
   pairs with a 3-round approximation (correct key ≈ 16σ), and invert the key schedule.
5. **Both directions** — decryption has the structure of encryption, so the same attack also runs
   with plaintext and ciphertext swapped.

The key came out of the reverse direction with r1 = 3.44 after about 40 minutes on 8 × RTX 5090
(about 1 hour for a full direction). An average brute-force search would take about 35 days on the
same machine. Earlier runs with stricter thresholds had narrowly missed it; the statistics of the
correct key show why ([results/experiments/actual_key_statistics.txt](results/experiments/actual_key_statistics.txt)).

## Reproduce

```sh
make -C src && make -C src check        # build (gcc with OpenMP) and verify key and data
bash scripts/reproduce_quick.sh         # CPU, about 5-10 minutes: re-runs stages 1-3 on the DK7 value
                                        # that contained the key and recovers it

make -C src gpu ARCH=sm_120             # CUDA (sm_86 for RTX 30xx, sm_89 for RTX 40xx)
bash scripts/run_search.sh rev 3.44     # full search of one direction on all GPUs (about 1 hour on 8 x RTX 5090)
```

## Contents

```
data/      challenge plaintext and ciphertext (.bin, hex, decoded text), printed CBC-MAC values, key
docs/      description of the attack (English / Japanese), references
src/       attack (stages 1-3, CPU and CUDA), key verification, test-data tools, supporting experiments
scripts/   full search driver and quick reproduction
results/   logs of all search runs (including the one that found the key) and experimental results
```

## Data source

The challenge data were transcribed by hand from the appendix figures 6–9 of
松井甲子雄『コンピュータによる暗号解読法入門』(Kineo Matsui, "Introduction to cryptanalysis by
computer", in Japanese), Morikita Publishing, 1990, ISBN 4-627-83350-4, and checked against the 64 CBC-MAC
values printed there (all match). See [data/README.md](data/README.md).

## License

The code and documents written for this repository are released into the public domain under the
[Unlicense](LICENSE). The challenge data in `data/` are reproduced from the publication above for
verification and research; the text contained in the plaintext is the beginning of the FEAL-8 paper
by Miyaguchi, Shiraishi and Shimizu (NTT). The Unlicense does not apply to them.
