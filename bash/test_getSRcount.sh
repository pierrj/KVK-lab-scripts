#!/bin/bash
pgbt9=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pGBT9)
puc19=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pUC19_tweaked)
pbr322=$(/global/home/users/pierrj/git/bash/get_SR_count.sh ${1} pBR322_tweaked)
echo -e ${1}'\t'$pgbt9'\t'$puc19'\t'$pbr322 >> /global/scratch/users/pierrj/eccDNA/pipeline_tests/normalization/SR_counts