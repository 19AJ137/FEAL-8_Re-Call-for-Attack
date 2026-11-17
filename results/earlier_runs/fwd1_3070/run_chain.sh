#!/bin/bash
# 順方向（run_real.sh）が鍵を見つけずに終わったら、逆方向（平文と暗号文を入れ替え）を同じ手順で回す
cd /root/feal
while [ ! -f run_real.done ]; do sleep 60; done
[ -f FOUND_KEY.txt ] && exit 0
OUT=real_rev_cands.txt
rm -f $OUT tier1_rev_all.txt real_rev_s3.log run_rev.done
date -u > real_rev_start.txt
./attackC_gpu real ciphertext.bin plaintext.bin 0 4194304 $OUT 3.5 6.5 > real_rev_gpu.log 2>&1 &
GPID=$!
LAST=0; DONE=0
while true; do
  sleep 300
  kill -0 $GPID 2>/dev/null || DONE=1
  N=$(wc -l < $OUT 2>/dev/null || echo 0)
  if [ "$N" -gt "$LAST" ]; then
    tail -n +$((LAST+1)) $OUT | head -n $((N-LAST)) | awk "{z=\$6; if (z<0) z=-z; if (z>=7.3) print}" > tier_rev_new.txt
    LAST=$N
    if [ -s tier_rev_new.txt ]; then
      cat tier_rev_new.txt >> tier1_rev_all.txt
      echo "== $(date -u) lines=$N tier_new=$(wc -l < tier_rev_new.txt)" >> real_rev_s3.log
      REV=1 RECOVER=1 OMP_NUM_THREADS=200 ./stage3 real ciphertext.bin plaintext.bin tier_rev_new.txt >> real_rev_s3.log 2>&1
      if grep -q "^KEY " real_rev_s3.log; then grep "^KEY " real_rev_s3.log > FOUND_KEY.txt; echo reverse >> FOUND_KEY.txt; date -u >> FOUND_KEY.txt; kill $GPID; break; fi
    fi
  fi
  [ $DONE = 1 ] && break
done
date -u > run_rev.done
