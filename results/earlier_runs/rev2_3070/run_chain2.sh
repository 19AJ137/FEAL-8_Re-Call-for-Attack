#!/bin/bash
# 順方向（現行、run_real.sh）の GPU が終わったら、2族版の逆方向（r1=3.9175, rpre=2.486, r2=8.5）を回す。
# 候補は全部（1方向 100 件程度）段階3（REV=1）へ。どちらかで鍵が出たら止める。
cd /root/feal
FP=$(cat real_fwd_gpu.pid)
while kill -0 $FP 2>/dev/null; do sleep 60; [ -f FOUND_KEY.txt ] && exit 0; done
[ -f FOUND_KEY.txt ] && exit 0
OUT=real_rev2_cands.txt
rm -f $OUT real_rev2_s3.log run_rev2.done
date -u > real_rev2_start.txt
./sk/attackC2_gpu real ciphertext.bin plaintext.bin 0 4194304 $OUT 3.9175 2.486 8.5 > real_rev2_gpu.log 2>&1 &
GPID=$!
echo $GPID > real_rev2_gpu.pid
LAST=0; DONE=0
while true; do
  sleep 300
  if [ -f FOUND_KEY.txt ]; then kill $GPID 2>/dev/null; break; fi
  kill -0 $GPID 2>/dev/null || DONE=1
  N=$(wc -l < $OUT 2>/dev/null || echo 0)
  if [ "$N" -gt "$LAST" ]; then
    tail -n +$((LAST+1)) $OUT | head -n $((N-LAST)) > rev2_new.txt
    LAST=$N
    echo "== $(date -u) lines=$N new=$(wc -l < rev2_new.txt)" >> real_rev2_s3.log
    REV=1 RECOVER=1 OMP_NUM_THREADS=120 ./stage3 real ciphertext.bin plaintext.bin rev2_new.txt >> real_rev2_s3.log 2>&1
    if grep -q "^KEY " real_rev2_s3.log; then grep "^KEY " real_rev2_s3.log > FOUND_KEY.txt; echo "reverse (2-family)" >> FOUND_KEY.txt; date -u >> FOUND_KEY.txt; kill $GPID; break; fi
  fi
  [ $DONE = 1 ] && break
done
date -u > run_rev2.done
