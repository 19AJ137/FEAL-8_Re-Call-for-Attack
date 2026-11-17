#!/bin/bash
# 2族版 r1=3.44（現行 3.0 相当、1方向の段階1・2 通過率 95%）。逆方向 → 順方向、各 8 GPU。
# 候補は 2 分ごとに、前回（r1=3.9175）で調べ済みのものを除いて、緩めた段階3（T1=1.0, T2=4.5）へ。鍵が出たら全部止める。
cd /root/feal
G=8; R1=3.44
s3() { # $1=候補 $2=REV $3=ログ
  if [ "$2" = 1 ]; then T1=1.0 T2=4.5 REV=1 RECOVER=1 OMP_NUM_THREADS=240 nice -n 5 ./stage3 real ciphertext.bin plaintext.bin "$1" >> "$3" 2>&1
  else T1=1.0 T2=4.5 RECOVER=1 OMP_NUM_THREADS=240 nice -n 5 ./stage3 real plaintext.bin ciphertext.bin "$1" >> "$3" 2>&1; fi
  if grep -q "^KEY " "$3"; then { grep "^KEY " "$3"; echo "$3"; date -u; } > FOUND_KEY.txt; return 1; fi
  return 0
}
phase() { # $1=名前 $2=REV $3=調べ済み一覧
  name=$1; rev=$2; seen=$3; S=0; E=4194304
  if [ "$rev" = 1 ]; then PF=ciphertext.bin; CF=plaintext.bin; else PF=plaintext.bin; CF=ciphertext.bin; fi
  per=$(( (E - S + G - 1) / G )); pids=""
  date -u > ${name}_start.txt
  for g in $(seq 0 $((G-1))); do
    s=$((S + g*per)); c=$per; [ $((s + c)) -gt $E ] && c=$((E - s))
    rm -f ${name}_cands_$g.txt
    CUDA_VISIBLE_DEVICES=$g ./sk/attackC2_gpu real $PF $CF $s $c ${name}_cands_$g.txt $R1 2.486 8.5 > ${name}_gpu_$g.log 2>&1 &
    pids="$pids $!"
  done
  declare -A last
  while true; do
    sleep 120
    alive=0; for p in $pids; do kill -0 $p 2>/dev/null && alive=1; done
    : > ${name}_raw.txt
    for g in $(seq 0 $((G-1))); do
      f=${name}_cands_$g.txt; [ -f $f ] || continue
      n=$(wc -l < $f); l=${last[$g]:-0}
      [ $n -gt $l ] && tail -n +$((l+1)) $f | head -n $((n-l)) >> ${name}_raw.txt
      last[$g]=$n
    done
    if [ -s ${name}_raw.txt ]; then
      cat ${name}_raw.txt >> ${name}_cands_all.txt
      awk "NR==FNR {s[\$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5]=1; next} !((\$1\" \"\$2\" \"\$3\" \"\$4\" \"\$5) in s)" $seen ${name}_raw.txt > ${name}_new.txt
      echo "== $(date -u) raw=$(wc -l < ${name}_raw.txt) new=$(wc -l < ${name}_new.txt)" >> ${name}_s3.log
      if [ -s ${name}_new.txt ]; then
        if ! s3 ${name}_new.txt $rev ${name}_s3.log; then kill $pids 2>/dev/null; return 1; fi
      fi
    fi
    [ $alive = 0 ] && break
  done
  date -u > ${name}_done.txt
  return 0
}
phase rev3 1 seen_rev.txt || exit 0
phase fwd3 0 seen_fwd.txt || exit 0
date -u > run9.done
