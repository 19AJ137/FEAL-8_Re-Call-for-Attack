#!/bin/bash
# 8 GPU で 2族版（r1=3.9175, rpre=2.486, r2=8.5）を回す。
# 1) 旧サーバの逆方向の候補 14 件に段階3  2) 逆方向の残り [702464, 2^22)  3) 順方向 [0, 2^22)
# 候補は 2 分ごとに全件を段階3へ。鍵が出たら FOUND_KEY.txt に書いて全部止める。
cd /root/feal
G=8
s3() { # $1=候補ファイル $2=REV(0/1) $3=ログ
  if [ "$2" = 1 ]; then REV=1 RECOVER=1 OMP_NUM_THREADS=300 ./stage3 real ciphertext.bin plaintext.bin "$1" >> "$3" 2>&1
  else RECOVER=1 OMP_NUM_THREADS=300 ./stage3 real plaintext.bin ciphertext.bin "$1" >> "$3" 2>&1; fi
  if grep -q "^KEY " "$3"; then { grep "^KEY " "$3"; echo "$3"; date -u; } > FOUND_KEY.txt; return 1; fi
  return 0
}
phase() { # $1=名前 $2=REV $3=開始 $4=終了
  name=$1; rev=$2; S=$3; E=$4
  if [ "$rev" = 1 ]; then PF=ciphertext.bin; CF=plaintext.bin; else PF=plaintext.bin; CF=ciphertext.bin; fi
  per=$(( (E - S + G - 1) / G )); pids=""
  date -u > ${name}_start.txt
  for g in $(seq 0 $((G-1))); do
    s=$((S + g*per)); c=$per; [ $((s + c)) -gt $E ] && c=$((E - s)); [ $c -le 0 ] && continue
    rm -f ${name}_cands_$g.txt
    CUDA_VISIBLE_DEVICES=$g ./sk/attackC2_gpu real $PF $CF $s $c ${name}_cands_$g.txt 3.9175 2.486 8.5 > ${name}_gpu_$g.log 2>&1 &
    pids="$pids $!"
  done
  echo "$pids" > ${name}_pids.txt
  declare -A last
  while true; do
    sleep 120
    alive=0; for p in $pids; do kill -0 $p 2>/dev/null && alive=1; done
    : > ${name}_new.txt
    for g in $(seq 0 $((G-1))); do
      f=${name}_cands_$g.txt; [ -f $f ] || continue
      n=$(wc -l < $f); l=${last[$g]:-0}
      [ $n -gt $l ] && tail -n +$((l+1)) $f | head -n $((n-l)) >> ${name}_new.txt
      last[$g]=$n
    done
    if [ -s ${name}_new.txt ]; then
      cat ${name}_new.txt >> ${name}_cands_all.txt
      echo "== $(date -u) new=$(wc -l < ${name}_new.txt)" >> ${name}_s3.log
      if ! s3 ${name}_new.txt $rev ${name}_s3.log; then kill $pids 2>/dev/null; return 1; fi
    fi
    [ $alive = 0 ] && break
  done
  date -u > ${name}_done.txt
  return 0
}
echo "== $(date -u) old reverse candidates" >> old_rev2_s3.log
s3 old_rev2_cands.txt 1 old_rev2_s3.log || exit 0
phase rev2b 1 702464 4194304 || exit 0
phase fwd2 0 0 4194304 || exit 0
date -u > run8.done
