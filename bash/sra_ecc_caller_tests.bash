outdir=SRR6315407
accession=SRR6315407
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/prefetch ${accession} -O ${outdir}
/global/home/users/pierrj/scripts/sratoolkit.2.10.4-centos_linux64/bin/fasterq-dump -e 24 -O ${outdir} -t ${outdir}/tmp ${outdir}/${accession}.sra