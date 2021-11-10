#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#!/bin/bash
#SBATCH --job-name=Pacbio_v_illumina_eccs
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio

output_table=/global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/illumina_v_pacbio_out

if [ -f "${output_table}" ]; then
    rm ${output_table}
fi

while read sample; do
    # input illumina called eccdnas, illumina split reads (with all false positives) and pacbio eccdnas and compare them
    python /global/home/users/pierrj/git/python/Pacbio_v_illumina_eccs.py /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample}_pacbio_eccs.bed \
        /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample}/${sample}.ecc_caller_out.details.nolowq.txt \
        /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/all_false_positives/${sample}/${sample}.mapq0_and_uniq_splitreads.allfalsepositives.bed \
        ${sample}.Pacbio_v_illumina ${sample}
        cat ${sample}.Pacbio_v_illumina >> ${output_table}
done < /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/mapfile