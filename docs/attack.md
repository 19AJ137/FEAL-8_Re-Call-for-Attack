# How the key was recovered

[日本語版](attack.ja.md)

This document describes the attack that recovered the key of the FEAL-8 preliminary problem
from the 1989 "Call for Attack", how it was run, and how much it cost. The source code is in
[`../src`](../src) and the raw logs are in [`../results`](../results).

- Key: `A5EA6EC2FB0CC741`
- Data: 1024 known plaintext/ciphertext blocks, ECB mode, one secret 64-bit key (981 distinct pairs)
- Attack: known-plaintext linear cryptanalysis with partitioning, run in both directions
- Cost of the successful search: about 40 minutes on 8 × RTX 5090 (one direction, about 67% of the
  space searched before the key appeared). A complete search of one direction takes about 1 hour on the
  same machine; an average brute-force search (2^63 keys) would take about 35 days.

## 1. The problem

The preliminary problem of the call (Miyaguchi, CRYPTO '89) gives 8192 bytes of plaintext and the
8192 bytes of ciphertext obtained by encrypting it with FEAL-8 in ECB mode under a secret fixed key.
The task is to find the key.

- 8192 bytes = 1024 blocks of 64 bits, i.e. 2^10 known plaintexts. Chosen-plaintext attacks
  (differential cryptanalysis and its relatives) do not apply.
- 43 blocks are repeats, which leaves N = 981 distinct pairs. Only distinct pairs are used.
- 725 of the distinct blocks are Shift-JIS/ASCII text and 256 are random-looking bytes. The text
  does not hurt the attack: with these plaintexts, the statistics of the correct and of wrong keys
  are indistinguishable from those obtained with 981 uniformly random plaintexts
  ([`../results/experiments/feasibility`](../results/experiments/feasibility)).

Published known-plaintext attacks on 8-round FEAL that work with as few as 2^10 pairs exist only for
FEAL-8X (FEAL-8 with a 128-bit key schedule): the estimate of Biham and Carmeli (SAC 2014, Section 5),
and the multiple linear cryptanalysis of Sakikoyama et al. (IEICE 2016), which succeeds with
non-negligible probability. FEAL-8 and FEAL-8X differ only in the key schedule, so their core
applies to FEAL-8. Neither method alone was practical here (Section 7), so the attack below
combines their ingredients.

## 2. Notation

- 32-bit words are written with byte 0 as the most significant byte. Rounds are numbered 0 to 7,
  and (L_i, R_i) is the state entering round i.
- Equivalent representation (Biham–Carmeli, Section 2): F(A, K) = N(M(A) ⊕ (0, K_hi, K_lo, 0)) with
  M(A) = (A0, A0⊕A1, A2⊕A3, A3), where N consists of the four S-boxes. The input and output whitening
  keys can be absorbed into 32-bit effective round keys:
  - plaintext side: EK_i = M(K89 ⊕ Kab) ⊕ (0, K_i, 0) for even i, M(K89) ⊕ (0, K_i, 0) for odd i
  - ciphertext side: DK_i = M(Kcd) ⊕ (0, K_i, 0) for even i, M(Kcd ⊕ Kef) ⊕ (0, K_i, 0) for odd i
- With these keys, the first and last rounds can be computed from the data alone, up to constants
  that only flip the sign of the linear approximations:
  - L1 = PL ⊕ PR, R1 = PL ⊕ N(M(L1) ⊕ EK0)
  - R7 = CL ⊕ CR, L7 = CL ⊕ N(M(R7) ⊕ DK7)

## 3. Linear approximations

The attack uses the 5-round approximation over rounds 1–5 of Biham and Carmeli (the first five rounds
of their "Approximation 1"):

| Round | Output mask ΓY | Input mask ΓX | Correlation |
|---|---|---|---|
| 1 | 00010504 | 00000100 | 1 |
| 2 | 04031104 | 00010504 | 1/2 |
| 3 | 0 | 0 | 1 |
| 4 | 04031104 | 00010504 | 1/2 |
| 5 | 00010504 | 00000100 | 1 |

The total correlation is 1/4 (bias 2^-3). On the state, the masks are (00010504, 04031004) on
(L1, R1) and (04031004, 00010504) on (L6, R6).

The appendix of Sakikoyama et al. lists eight 6-round approximations in two families. The second
family replaces the round-2 and round-4 output mask 04031104 by 04011904; its 5-round part also has
correlation 1/4, with masks (00010504, 04011804) on (L1, R1) and (04011804, 00010504) on (L6, R6). Both
families use exactly the same key bits. For the correct key their noise is almost independent (the
correlation of the two statistics' absolute values is −0.06), so using both roughly doubles the
signal-to-noise ratio in squared terms. Which family a statistic refers to is written as a subscript:
z0 for the first family, z1 for the second.

Each family gives one parity bit per pair on the plaintext side,
p_f = ⟨00010504, L1⟩ ⊕ ⟨a_f, R1⟩, and one on the ciphertext side that also needs the round-6
F output, which depends on DK6. With N = 981 pairs, the correct key gives a statistic of about
√N / 4 ≈ 7.8 standard deviations for each family.

Key bits involved:

| Key | Bits | How they are found |
|---|---|---|
| EK0 | 15 (bytes 1 and 2; the top bit of byte 2 is normalised) | stage 1 (all 2^15) |
| DK7 | 22 of 30 (mask 03FFFF0F) | stage 1 (all 2^22) |
| DK7 | remaining 8 (byte 0 bits 2–6, byte 3 bits 4–6) | stage 2 (e8) |
| DK6 | 8 (κ, the partitioning key) | stage 1, all 256 at once with a Walsh–Hadamard transform |
| DK6 | remaining 7 (r) | stage 2 |

The top bits of bytes 0 and 3 of DK7 only act linearly and are not needed until the key schedule is
solved. In total a full hypothesis has 15 + 30 + 15 = 60 bits, the hypothesis space of
Biham–Carmeli's Section 5 attack. Testing all 2^60 directly costs about 2^64 word operations;
stages 1 and 2 avoid this.

## 4. The attack

### Stage 1: partitioning filter (2^45 hypotheses)

Once DK7 (22 bits) is guessed, L7 is known up to constants, so the round-6 F input is known up to
the corresponding bits of DK6. Eight of those input bits (bits 2–4 of the second byte and bits 2–6 of
the third byte of M(L7) ⊕ DK6) determine a partition of the pairs into 256 groups. For each group g
the correlation e(g) of the round-6 output parity is computed exactly once and for all. This is the
partitioning technique of Biham and Carmeli (1 bit) generalised by Leurent (4 bits) to 8 bits. Its
capacity Σ e(g)^2 / 256 is 0.200, 12.8 times that of a single 6-round approximation (1/64).

For each (DK7, EK0) the signed sums of p ⊕ c within each group are collected, and the weighted
statistic

    z_B(κ) = Σ_g e(g ⊕ κ) S_g / sqrt(Σ_g e(g ⊕ κ)^2 N_g)

is obtained for all 256 κ at once with two Walsh–Hadamard transforms. This is done for both
families, and a hypothesis survives if

    z_B0^2 + z_B1^2 ≥ r1^2.

For the correct key E|z_B| ≈ 3.5 per family. For wrong keys z_B0 and z_B1 are independent standard
normal variables (measured correlation 0.0001), so the survival rate is exp(−r1^2 / 2):
r1 = 3.9175 keeps 4.6 × 10^-4 of the 2^45 hypotheses (the same number as |z| ≥ 3.5 with one family),
and r1 = 3.44 keeps 2.7 × 10^-3 (the same as |z| ≥ 3.0).

### Stage 2: full 5-round statistic (2^15 completions per survivor)

For each surviving hypothesis, the remaining 8 bits of DK7 (e8) and 7 bits of DK6 (r) are completed.
The round-6 output parity is then known exactly, and the method-A statistic of Biham–Carmeli
Section 5 is computed for both families:

    z_Af = (number of pairs where the approximation holds − number where it fails) / sqrt(N)

The plaintext-side bits for all 2^15 EK0 values are precomputed as 1024-bit vectors; the
ciphertext-side vectors are built from lookup tables and bit-matrix transposition; the comparison is
XOR plus popcount. Hypotheses are first tested on the first 256 pairs
(z4_0^2 + z4_1^2 ≥ 2.486^2 · 256, which keeps the same fraction as 2σ with one family), and
candidates with

    r = sqrt(z_A0^2 + z_A1^2) ≥ r2 = 8.5

are passed to stage 3. Over the whole space, about 120 wrong candidates per direction are expected
(measured: 65 and 73 at r1 = 3.9175), while the correct key gives r ≈ 11.

### Stage 3: remaining subkeys and the 64-bit key

Stage 3 (`stage3.c`) shifts the 5-round approximation by one round in each direction and uses the
3-round probability-one-half approximation to tie the two sides together:

| Step | Approximation | Guessed | Result |
|---|---|---|---|
| (a) | rounds 2–6 | bytes 0 and 3 of EK0 (14 bits), a normalisation bit, middle bytes of EK1 (15 bits) | all of EK0, middle of EK1 |
| (c) | rounds 0–4 | bytes 0 and 3 of DK6 (14 bits), a normalisation bit, middle bytes of DK5 (15 bits) | all of DK6, middle of DK5 |
| (t) | 3-round approximation of rounds 2–4 (correlation 1/2) | nothing | test of each (a) × (c) pair; the correct pair scores about 15.7σ |
| (e) | key schedule | 7 ambiguous top bits and three unknown bytes, exhaustively (2^7 × 2^24) | 64-bit key, checked on two blocks |

Steps (a) and (c) use the same two-step method as stages 1 and 2 (partitioning, then completion),
keep up to 4096 hypotheses with |z| ≥ T2, and pass pairs with |z3| ≥ T3 = 8 to step (e), in decreasing
order of z_a^2 + z_c^2 + z_3^2. Step (e) solves the key schedule byte by byte (fK can be inverted in
its second argument) and takes less than 0.5 s per pair.

- Default thresholds T1 = 2.0, T2 = 5.0: 215 of 240 known keys recovered (89.6%), about 40 s of thread
  time per candidate.
- Relaxed thresholds T1 = 1.0, T2 = 4.5 (used in the final search): 239 of 240 (99.6%), about 55 s per
  candidate.
- Wrong candidates almost never pass (t): 3 of 126,000 in the first search, all with z3 ≈ 8, whereas the
  correct key gives |z3| ≈ 15–17.

### Both directions

FEAL-8 decryption has the same structure as encryption with the round keys in reverse order and the
whitening keys swapped (K89 ↔ Kcd, Kab ↔ Kef). Feeding (ciphertext, plaintext) as (plaintext,
ciphertext) therefore runs the same attack on the "reverse cipher". Only step (e) changes
(`REV=1` in `stage3`). The statistics in the two directions come from different intermediate values,
so the two runs fail almost independently.

## 5. Success probability

The probability that the correct key passes stages 1 and 2 was measured with 30,000 random keys
encrypting the challenge plaintexts, using the true intermediate values
([`../results/experiments/two_family`](../results/experiments/two_family)):

| Variant | Forward | Reverse | Either |
|---|---|---|---|
| One family, \|z_B\| ≥ 3.5, \|z_A\| ≥ 7.3 | 42.0% | 42.2% | 66.5% |
| Two families, r1 = 3.9175, r2 = 8.5 | 87.8% | 88.0% | 98.6% |
| Two families, r1 = 3.44, r2 = 8.5 | 95.0% | 95.0% | 99.8% |

Stage 3 with the relaxed thresholds adds a factor of about 0.996.

For comparison, the method of Sakikoyama et al. on its own (eight 6-round approximations, keeping the
top 2^24 of 2^37 candidates) placed the correct key in the top 2^24 for 16 of 100 keys with these
plaintexts, about the same as the 12% in their Table 3.

## 6. The actual search

| Run | Stages 1–2 | Candidates to stage 3 | Passed (t) | Result |
|---|---|---|---|---|
| Forward, one family, \|z_B\| ≥ 3.5 | RTX 3070, 8.6 h | 126,000 (\|z\| ≥ 7.3) | 3 (z3 ≈ 8.0) | no key |
| Reverse, two families, r1 = 3.9175 | RTX 3070 for 0–702,463, then 8 × RTX 5090, 19 min | 73 | 0 | no key |
| Forward, two families, r1 = 3.9175 | 8 × RTX 5090, 22 min | 65 | 0 | no key |
| Same 138 candidates, relaxed stage 3 | CPU, 2 min | 138 | 0 | no key |
| **Reverse, two families, r1 = 3.44** | **8 × RTX 5090, key after 40 min** | **new candidates every 2 min** | **3 (correct: z3 = −16.25)** | **A5EA6EC2FB0CC741** |

Details of the successful run ([`../results/key_recovery`](../results/key_recovery)):

- Started 2026-11-17 04:27 UTC; the key was printed at 05:06:31 UTC (14:06 JST).
- Stage-2 candidate: DK7 index `34fd3b`, EK0 `6f7a`, κ `6b`, e8 `aa`, r `45`, with
  r = 10.17 (z0 = 7.44, z1 = 6.93).
- Stage 3: (a) best z = 7.12, (c) best z = −8.72, (t) z3 = −16.25 (the largest value from a wrong
  candidate in all runs was about 11.5), (e) key found from the first pair.

Statistics of the correct key, computed afterwards from the true intermediate values
([`../results/experiments/actual_key_statistics.txt`](../results/experiments/actual_key_statistics.txt)):

| | Forward | Reverse |
|---|---|---|
| Stage 1, sqrt(z_B0² + z_B1²) | 3.80 | 3.61 |
| Stage 1, \|z_B0\| (one family) | 3.47 | 2.31 |
| Stage 2, sqrt(z_A0² + z_A1²) | 10.30 | 10.17 |

This explains the earlier misses: the one-family run failed by 0.03 at stage 1 (3.47 < 3.5), and both
two-family runs at r1 = 3.9175 failed at stage 1 as well. With r1 = 3.44 both directions would have
found the key.

## 7. Cost

Measured speeds of stages 1–2 (two families):

| Hardware | r1 = 3.9175 | r1 = 3.44 |
|---|---|---|
| RTX 3070, per DK7 | 0.0137 s (one direction ≈ 16 h) | — |
| RTX 5090, per DK7 per GPU | 0.0023 s (8 GPUs: ≈ 20 min per direction) | 0.0067 s (8 GPUs: ≈ 59 min per direction) |
| CPU, thread time per DK7 | ≈ 3 s (≈ 145 core-days per direction) | ≈ 6 s (≈ 290 core-days per direction) |

Stage 3 costs about 1 minute of thread time per candidate. In total the project used about 12 GPU-hours
on an RTX 3070 and about 11 GPU-hours on RTX 5090s.

For comparison, a bitsliced brute-force search (`bs_gpu.cu`) runs at 7.6 × 10^10 keys/s on an RTX 3070
and 3.8 × 10^11 keys/s on an RTX 5090, i.e. 2^63 keys (the average) in 3.8 years and 0.77 years per
GPU, or about 35 days on 8 × RTX 5090.

Other approaches that were evaluated and discarded:

- Biham–Carmeli Section 5 as published (all 2^60 hypotheses): about 2^64 word operations, roughly
  300 core-years, with a 91% chance that the correct key reaches 6.5σ.
- Sakikoyama et al. alone: 16% chance of the correct key being in the top 2^24, and 2^24 candidates
  at their reported 0.2 s each is about 39 days.
- Differential-type attacks: the 480,690 plaintext pairs contain no repeated differences useful for
  any characteristic (the most frequent 64-bit difference appears twice).

## 8. Verification

- `src/verify_key.py`, written independently of the C code and checked against the FEAL-8 test vector
  (key 0123456789ABCDEF, plaintext 0 → CEEF2C86F2490752), encrypts all 1024 plaintext blocks under
  A5EA6EC2FB0CC741 and reproduces all 1024 ciphertext blocks. The neighbouring key A5EA6EC2FB0CC740
  matches none.
- The C implementation (`src/mkct.c`, built on `src/feal8.h`) produces the same ciphertext file.
- The transcription of the data was checked against the 64 CBC-MAC values printed with the problem
  (`src/check_mac.py`, see [`../data/README.md`](../data/README.md)).
- The key satisfies C = E_K(P), as stated in the problem. (After the first searches failed, the
  alternative C = D_K(P) was also investigated; the attack works under either mapping, see
  [`../results/experiments/reverse_mapping_check`](../results/experiments/reverse_mapping_check).)

## References

See [references.md](references.md).
