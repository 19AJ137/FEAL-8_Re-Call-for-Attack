#!/bin/bash
cd /root/feal
seq 1001 1240 | xargs -P 120 -I{} sh -c "EFAST=1 OMP_NUM_THREADS=1 ./stage3 self plaintext.bin {} > s3true/{}.txt 2>&1"
echo done > runA.done
