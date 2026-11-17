#!/bin/bash
# run9.sh が鍵なしで終わったら（run9.done）、順方向 r1=3.0 を run10.sh で回す。鍵が出ていれば何もしない。
cd /root/feal
while [ ! -f run9.done ]; do
  [ -f FOUND_KEY.txt ] && exit 0
  if ! pgrep -f "^/bin/bash ./run9.sh" > /dev/null; then sleep 5; [ -f run9.done ] || exit 0; fi
  sleep 30
done
[ -f FOUND_KEY.txt ] && exit 0
cat all_fwd2.txt fwd3_cands_all.txt 2> /dev/null | cut -d" " -f1-5 > seen_fwd30.txt
exec ./run10.sh fwd30 0 3.0 seen_fwd30.txt
