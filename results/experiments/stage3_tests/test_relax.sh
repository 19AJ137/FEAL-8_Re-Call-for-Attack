#!/bin/bash
cd /root/feal
seq 1001 1240 | xargs -P 160 -I{} sh -c "T1=1.0 T2=4.5 EFAST=1 OMP_NUM_THREADS=1 nice -n 10 ./stage3 self plaintext.bin {} > s3relax2/{}.txt 2>&1"
date -u > s3relax2.done
