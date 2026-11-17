#!/bin/bash
cd /root/feal
seq 1001 1240 | xargs -P 120 -I{} sh -c "T2=5.0 EFAST=1 OMP_NUM_THREADS=1 ./stage3 self plaintext.bin {} > s3true5/{}.txt 2>&1"
echo done > runB.done
CANDOUT=cands_201.txt OMP_NUM_THREADS=240 ./attackC self plaintext.bin 201 2400 3.5 2.0 6.0 > c3.5_201_cands.log 2>&1
echo done > runB2.done
