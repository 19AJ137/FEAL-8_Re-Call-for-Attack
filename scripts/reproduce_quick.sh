#!/usr/bin/env bash
# Quick CPU reproduction of the step that found the key (about 5-10 minutes on a 2-core machine).
#
# Re-runs stages 1-2 for the single DK7 value in which the key was found (index 0x34fd3b,
# reverse direction, r1=3.44), passes the resulting candidates to stage 3, and verifies the
# recovered key against all 1024 blocks with an independent Python implementation.
#
#   bash scripts/reproduce_quick.sh            (needs gcc with OpenMP, make, python3)
#   OUT=dir bash scripts/reproduce_quick.sh    (keep the intermediate files in dir)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/src
DATA=$ROOT/data
OUT=${OUT:-$(mktemp -d)}
mkdir -p "$OUT"

make -C "$SRC" attackC2 stage3 >/dev/null

echo "== stages 1-2: reverse direction, DK7 index 0x34fd3b, r1=3.44 rpre=2.486 r2=8.5"
"$SRC/attackC2" real "$DATA/ciphertext.bin" "$DATA/plaintext.bin" 0x34fd3b 1 "$OUT/cands.txt" 3.44 2.486 8.5
echo "-- candidates (t e kappa e8 r  r-stat z0 z1):"
cat "$OUT/cands.txt"

echo "== stage 3: REV=1 T1=1.0 T2=4.5"
REV=1 RECOVER=1 T1=1.0 T2=4.5 "$SRC/stage3" real "$DATA/ciphertext.bin" "$DATA/plaintext.bin" "$OUT/cands.txt" \
    | tee "$OUT/stage3.log"
KEY=$(awk '/^KEY /{print $2; exit}' "$OUT/stage3.log")
if [ -z "$KEY" ]; then echo "key not recovered"; exit 1; fi

echo "== verification (independent implementation)"
python3 "$SRC/verify_key.py" "$KEY" "$DATA/plaintext.bin" "$DATA/ciphertext.bin"
echo "intermediate files: $OUT"
