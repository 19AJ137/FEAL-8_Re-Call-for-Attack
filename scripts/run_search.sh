#!/usr/bin/env bash
# Full key search with Method C (two-family version), the procedure that recovered the key.
#
#   bash scripts/run_search.sh DIR [R1] [NGPU]
#     DIR   rev = reverse direction (ciphertext.bin is fed in as the "plaintext"; the key was found this way)
#           fwd = forward direction
#     R1    stage-1 threshold on sqrt(z0^2 + z1^2)   (default 3.44; 3.9175 is about 3x cheaper but
#           passes fewer correct keys, see docs/attack.md)
#     NGPU  number of GPUs (default: all GPUs listed by nvidia-smi)
#
# Environment variables:
#   START, END       DK7 index range [START, END)   (default 0 .. 2^22)
#   WORK             output directory                (default ./search_<DIR>_<R1>)
#   POLL             seconds between candidate sweeps (default 120)
#   ATTACK           stage 1-2 program (default src/attackC2_gpu; src/attackC2 also works for small ranges)
#   OMP_NUM_THREADS  threads used by stage 3
#
# Stages 1-2 run on the GPUs, one contiguous slice of the 2^22 DK7 values per GPU. Every POLL seconds
# the new candidates (r >= 8.5) are passed to stage 3 (T1=1.0, T2=4.5). When stage 3 prints
# "KEY ...", the key is written to $WORK/FOUND_KEY.txt and the GPU jobs are stopped.
# On 8x RTX 5090, one direction at R1=3.44 takes about one hour.
#
# This is a cleaned-up version of results/key_recovery/run9.sh (the script that actually ran).
set -u
DIR=${1:?usage: run_search.sh rev|fwd [R1] [NGPU]}
R1=${2:-3.44}
NG=${3:-$(nvidia-smi -L 2>/dev/null | wc -l)}
[ "${NG:-0}" -ge 1 ] 2>/dev/null || NG=1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/src
DATA=$ROOT/data
ATTACK=${ATTACK:-$SRC/attackC2_gpu}
START=$((${START:-0}))
END=$((${END:-4194304}))
POLL=${POLL:-120}
WORK=${WORK:-$PWD/search_${DIR}_${R1}}

case $DIR in
  rev) PF=$DATA/ciphertext.bin; CF=$DATA/plaintext.bin; REVFLAG=1 ;;
  fwd) PF=$DATA/plaintext.bin;  CF=$DATA/ciphertext.bin; REVFLAG= ;;
  *) echo "DIR must be rev or fwd" >&2; exit 2 ;;
esac
[ -x "$ATTACK" ] || { echo "$ATTACK not found (run 'make gpu' in src/)" >&2; exit 2; }
[ -x "$SRC/stage3" ] || { echo "$SRC/stage3 not found (run 'make' in src/)" >&2; exit 2; }

mkdir -p "$WORK" && cd "$WORK" || exit 2
date -u > start.txt
per=$(( (END - START + NG - 1) / NG ))
pids=()
declare -a last
for ((g = 0; g < NG; g++)); do
  s=$((START + g * per)); c=$per
  ((s + c > END)) && c=$((END - s))
  ((c <= 0)) && continue
  rm -f "cands_$g.txt"
  CUDA_VISIBLE_DEVICES=$g "$ATTACK" real "$PF" "$CF" "$s" "$c" "cands_$g.txt" "$R1" 2.486 8.5 > "gpu_$g.log" 2>&1 &
  pids+=($!)
  last[$g]=0
done

stop() { kill "${pids[@]}" 2>/dev/null; }
trap 'stop; exit 130' INT TERM

while :; do
  sleep "$POLL"
  alive=0
  for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && alive=1; done
  : > new.txt
  for ((g = 0; g < NG; g++)); do
    f=cands_$g.txt
    [ -f "$f" ] || continue
    n=$(wc -l < "$f"); l=${last[$g]:-0}
    ((n > l)) && sed -n "$((l + 1)),${n}p" "$f" >> new.txt
    last[$g]=$n
  done
  if [ -s new.txt ]; then
    cat new.txt >> cands_all.txt
    echo "== $(date -u) new candidates: $(wc -l < new.txt)" >> stage3.log
    # stage3 treats any value of REV as "reverse", so the variable is set only for DIR=rev
    env ${REVFLAG:+REV=1} T1=1.0 T2=4.5 RECOVER=1 "$SRC/stage3" real "$PF" "$CF" new.txt >> stage3.log 2>&1
    if grep -q '^KEY ' stage3.log; then
      { grep -m1 '^KEY ' stage3.log; date -u; } > FOUND_KEY.txt
      stop
      cat FOUND_KEY.txt
      exit 0
    fi
  fi
  ((alive)) || break
done
date -u > done.txt
echo "search of [$START, $END) finished without recovering a key (see $WORK/stage3.log)"
exit 1
