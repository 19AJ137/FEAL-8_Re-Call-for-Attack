#!/bin/bash
export OMP_NUM_THREADS=240
for t in 3.0 3.5; do for s in 201 202 203; do
  ./attackC self plaintext.bin $s 2400 $t 2.0 6.0 > c${t}_$s.txt 2>&1
done; done
echo DONE > run1.done
