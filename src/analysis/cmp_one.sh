#!/bin/sh
# Compare attackC (one family) and attackC2 (two families) on one known key (seed). Run from src/:
#   seq 1001 1060 | xargs -P2 -n1 sh analysis/cmp_one.sh
s=$1; P=${P:-../data/plaintext.bin}
a=$(./analysis/attackC self "$P" "$s" 1 3.5 2.0 6.5 | grep "true key in candidates")
b=$(./attackC2 self "$P" "$s" 1 3.9175 2.486 8.5 | grep "true key in candidates")
echo "$s | $a | $b"
