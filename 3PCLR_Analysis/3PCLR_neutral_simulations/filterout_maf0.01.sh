#!/bin/bash
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
  freq = $8 / $9
  if (freq > 0.01 && freq < 0.99)
    print
}' $1 > $2
