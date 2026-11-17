#!/bin/bash
# 2族版の 1 方向を全空間で回す（8 GPU）。生存仮説の上限を広げた attackC2_gpu_c18（CAP=262144）を使う。
#   ./run10.sh 名前 REV r1 調べ済み一覧
# 候補は 2 分ごとに、調べ済み一覧（t e kap e8 r）にないものだけ緩めた段階3（T1=1.0, T2=4.5、両写像）へ。
# 鍵が出たら FOUND_KEY.txt に書いて GPU を止める。終われば ${名前}_done.txt。
cd /root/feal
name=$1; rev=$2; R1=$3; seen=$4; G=8; BIN=./sk/attackC2_gpu_c18
s3() { # $1=候補 $2=REV $3=ログ
  if [ "$2" = 1 ]; then T1=1.0 T2=4.5 REV=1 RECOVER=1 OMP_NUM_THREADS=240 nice -n 5 ./stage3 real ciphertext.bin plaintext.bin "$1" >> "$3" 2>&1
  else T1=1.0 T2=4.5 RECOVER=1 OMP_NUM_THREADS=240 nice -n 5 ./stage3 real plaintext.bin ciphertext.bin "$1" >> "$3" 2>&1; fi
  if grep -q "^KEY " "$3"; then { grep "^KEY " "$3"; echo "$3"; date -u; } > FOUND_KEY.txt; return 1; fi
  return 0
}
S=0; E=4194304
if [ "$rev" = 1 ]; then PF=ciphertext.bin; CF=plaintext.bin; else PF=plaintext.bin; CF=ciphertext.bin; fi
per=$(( (E - S + G - 1) / G )); pids=""
date -u > ${name}_start.txt
for g in $(seq 0 $((G-1))); do
  s=$((S + g*per)); c=$per; [ $((s + c)) -gt $E ] && c=$((E - s))
  rm -f ${name}_cands_$g.txt
  CUDA_VISIBLE_DEVICES=$g $BIN real $PF $CF $s $c ${name}_cands_$g.txt $R1 2.486 8.5 > ${name}_gpu_$g.log 2>&1 &
  pids="$pids $!"
done
echo "$pids" > ${name}_pids.txt
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
      cut -d" " -f1-5 ${name}_new.txt >> $seen
      if ! s3 ${name}_new.txt $rev ${name}_s3.log; then kill $pids 2>/dev/null; exit 0; fi
    fi
  fi
  [ $alive = 0 ] && break
done
date -u > ${name}_done.txt
