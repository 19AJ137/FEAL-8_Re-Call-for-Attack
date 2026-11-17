#!/bin/bash
# 本番：順方向・段階1しきい値 3.5。GPU で段階1・2、CPU で |z|>=7.3 の候補に段階3を逐次実行
cd /root/feal
OUT=real_fwd_cands.txt
rm -f $OUT FOUND_KEY.txt tier1_all.txt real_fwd_s3.log run_real.done
date -u > real_fwd_start.txt
./attackC_gpu real plaintext.bin ciphertext.bin 0 4194304 $OUT 3.5 6.5 > real_fwd_gpu.log 2>&1 &
GPID=$!
echo $GPID > real_fwd_gpu.pid
LAST=0; DONE=0
while true; do
  sleep 300
  kill -0 $GPID 2>/dev/null || DONE=1
  N=$(wc -l < $OUT 2>/dev/null || echo 0)
  if [ "$N" -gt "$LAST" ]; then
    tail -n +$((LAST+1)) $OUT | head -n $((N-LAST)) | awk "{z=\$6; if (z<0) z=-z; if (z>=7.3) print}" > tier_new.txt
    LAST=$N
    if [ -s tier_new.txt ]; then
      cat tier_new.txt >> tier1_all.txt
      echo "== $(date -u) lines=$N tier_new=$(wc -l < tier_new.txt)" >> real_fwd_s3.log
      RECOVER=1 OMP_NUM_THREADS=200 ./stage3 real plaintext.bin ciphertext.bin tier_new.txt >> real_fwd_s3.log 2>&1
      if grep -q "^KEY " real_fwd_s3.log; then grep "^KEY " real_fwd_s3.log > FOUND_KEY.txt; date -u >> FOUND_KEY.txt; kill $GPID; break; fi
    fi
  fi
  [ $DONE = 1 ] && break
done
date -u > run_real.done
