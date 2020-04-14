#!/bin/bash
ls parallel.confirmed.*.bed | awk '{print substr($0,19,5)}' | sort | uniq
while read sample; do normalize_factor=$(wc -l parallel.confirmed.${sample}.bed | awk '{print $1/100000}'); awk -v N=$normalize_factor '{sum+=$3} NR%100==0 {print sum/100/N; sum=0}' coverage.parallel.confirmed.${sample}.bed > binned_normalized.coverage.SRs.${sample}.bed; done < mapfile