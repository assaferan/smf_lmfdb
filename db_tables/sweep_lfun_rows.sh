#!/bin/bash
# sweep_lfun_rows.sh INDIR OUTDIR: certified L-function rows for every in_*.txt written by
# `sweep_lfun.py inputs`, with the settings of the loaded pilots (40 bits).  On lava.
# Writes OUTDIR/<name>/{e.json,d.log,rows/...} and OUTDIR/summary.txt (one line per form).
IN=$1; OUT=$2; L=~/lfun; mkdir -p $OUT
for f in $IN/in_*.txt; do
  n=$(basename $f .txt); n=${n#in_}; D=$OUT/$n
  [ -s $D/rows/lfunc_lfunctions.txt ] && continue
  mkdir -p $D
  read lab k j fr < <(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['smf_label'], *d['weight'], d['fricke'])" $f.factors.json)
  $L/smf_driver_pr40.exe $f $D/e.json --self-dual --target-prec 40 > $D/d.log 2>&1
  if python3 $L/sweep/lfunc_rows.py rows --driver $D/e.json --factors $f.factors.json --smf-label $lab \
       --weight $k $j --aut-rep-type G --fricke $fr --accuracy 40 --outdir $D/rows > $D/rows.log 2>&1; then
    echo "$lab ok $(tail -1 $D/rows.log)" >> $OUT/summary.txt
  else
    echo "$lab FAILED $(tail -1 $D/rows.log)" >> $OUT/summary.txt
  fi
done
