#!/bin/bash

# Usage: ./script.sh input_file heterozygosity_threshold (e.g. 0.05)

input_file=$1
threshold=$2

awk -v threshold="$threshold" '{
  sum = 0
  for (i = 1; i <= NF; i++) if ($i > 0) sum += $i
  for (i = 1; i <= NF; i++) {
    if ($i > 0) printf "%.2f ", $i / sum
    else printf "0 "
  }
  printf "\n"
}' "$input_file" | awk -v t="$threshold" '{
  homo = 1 - t
  if ($1 > homo) print $0" AA";
  else if ($2 > homo) print $0" CC";
  else if ($3 > homo) print $0" GG";
  else if ($4 > homo) print $0" TT";
  else if ($1 > t && $2 > t) print $0" AC HET_trv";
  else if ($1 > t && $3 > t) print $0" AG HET_trn";
  else if ($1 > t && $4 > t) print $0" AT HET_trv";
  else if ($2 > t && $3 > t) print $0" CG HET_trv";
  else if ($2 > t && $4 > t) print $0" CT HET_trn";
  else if ($3 > t && $4 > t) print $0" GT HET_trv";
  else print $0" UNKNOWN"
}'

