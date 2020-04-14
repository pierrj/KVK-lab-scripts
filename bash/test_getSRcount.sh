#!/bin/bash
pgbt9=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pGBT9)
puc19=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pUC19_tweaked)
pbr322=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pBR322_tweaked)
touch SR_countable_${1}

echo -e '0.01''\t'$pgbt9 >> SR_countable_${1}
echo -e '0.1''\t'$puc19 >> SR_countable_${1}
echo -e '1''\t'$pbr322 >> SR_countable_${1}

slope=$(awk '{ x[NR] = $1; y[NR] = $2;
 sx += x[NR]; sy += y[NR]; 
 sxx += x[NR]*x[NR];
 sxy += x[NR]*y[NR];
}
END{
 det = NR*sxx - sx*sx;
 a = (NR*sxy - sx*sy)/det;
 b = (-sx*sxy+sxx*sy)/det;
 print a/100000;
# for(i=1;i<=NR;i++) print x[i],a*x[i]+b;
}' SR_countable_${1})

echo -e ${1}'\t'$slope >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/SRcount_slopes