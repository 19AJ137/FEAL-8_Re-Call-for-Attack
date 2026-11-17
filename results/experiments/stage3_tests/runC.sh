#!/bin/bash
cd /root/feal
awk "\$7==0 && (\$6>=6.5 || \$6<=-6.5)" cands_201.txt | shuf -n 1200 --random-source=cands_201.txt > wrong_sample.txt
awk "\$7==0" cands_201.txt | head -300 > top_wrong.txt
awk "\$7==1" cands_201.txt | head -40 > true_cluster.txt
OMP_NUM_THREADS=240 ./stage3 self plaintext.bin 201 wrong_sample.txt > s3_wrong.log 2>&1
OMP_NUM_THREADS=240 ./stage3 self plaintext.bin 201 top_wrong.txt > s3_topwrong.log 2>&1
RECOVER=1 OMP_NUM_THREADS=240 ./stage3 self plaintext.bin 201 true_cluster.txt > s3_truecl.log 2>&1
echo done > runC.done
rm -rf s3true6; mkdir s3true6
seq 1001 1240 | xargs -P 120 -I{} sh -c "EFAST=1 OMP_NUM_THREADS=1 ./stage3 self plaintext.bin {} > s3true6/{}.txt 2>&1"
echo done > runC2.done
