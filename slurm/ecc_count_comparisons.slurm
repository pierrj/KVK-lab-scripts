#!/bin/bash
#SBATCH --job-name=ecc_count_comparisons
#SBATCH --account=fc_kvkallow
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/

ecc_count_output=ecc_count_output_moryzae.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/

cd /global/scratch/users/pierrj/eccDNA/stress_experiments/rice_control/

ecc_count_output=ecc_count_output_myrice.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/ORSA_w_organelles.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run

ecc_count_output=ecc_count_output_human_muscle.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/GRCh37.p13.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_muscle

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2018_moller/full_run

ecc_count_output=ecc_count_output_human_leukocytes.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/GRCh37.p13.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_leukocytes

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run

ecc_count_output=ecc_count_output_ara_wt.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/TAIR10.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_ENA_ara_wt

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run

ecc_count_output=ecc_count_output_ara_epi.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/TAIR10.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_ENA_ara_epi

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run

ecc_count_output=ecc_count_output_ory_callus.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/ORSA_w_organelles.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_ENA_ory_callus

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run

ecc_count_output=ecc_count_output_ory_leaf.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/ORSA_w_organelles.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_ENA_ory_leaf

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2017_lanciano/full_run

ecc_count_output=ecc_count_output_ory_seed.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/ORSA_w_organelles.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_ENA_ory_seed

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run

ecc_count_output=ecc_count_output_yeast_d.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/Scer_S288C.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_d

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run

ecc_count_output=ecc_count_output_yeast_z.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/Scer_S288C.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_z

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run

ecc_count_output=ecc_count_output_yeast_g.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/Scer_S288C.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_g

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/


cd /global/scratch/users/pierrj/eccDNA/2015_moller/full_run

ecc_count_output=ecc_count_output_yeast_s.txt

echo $ecc_count_output

genome_file=/global/scratch/users/pierrj/references/Scer_S288C.fasta

if [ -f "${ecc_count_output}" ]; then
    rm ${ecc_count_output}
fi

while read sample;
do
    cd ${sample}
    echo $sample
    ecc_count=$(awk ' $6 != "lowq" ' ecccaller_output.${sample}.renamed.details.tsv | wc -l | awk '{print $1}')
    read_count=$(samtools view -c filtered.sorted.${sample}.bam | awk '{print $1/1000000}') ## million reads
    genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum/1000000;}' ${genome_file} ) ## megabasepairs
    total=$(awk -v ECC="$ecc_count" -v READ="$read_count" -v SIZE="$genome_size" 'BEGIN{print ECC/READ/SIZE}') ## ecc regions/millionreads/megabasepair
    echo -e $sample'\t'$total >> ../${ecc_count_output}
    cd ..
done < mapfile_s

cp ${ecc_count_output} /global/scratch/users/pierrj/eccDNA/comparative/ecc_counts/
