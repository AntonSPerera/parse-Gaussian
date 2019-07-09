#!/bin/sh
BEGIN=0
END=120
export GAUSS_EXEDIR=/cm/shared/gaussian/g09.e01/intel64-sandybridge/g09
export GAUSS_SCRDIR=$SCRATCH
mkdir -p $SCRATCH
for NUMBER in $(seq $BEGIN $END); do
  # IMPORTANT - Change job name, the variable $NUMBER goes from $BEGIN to $END
  sed -n '/Normal/p' out_"$NUMBER".log
  echo gau$NUMBER
done
