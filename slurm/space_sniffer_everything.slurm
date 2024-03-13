#!/bin/bash
#SBATCH --job-name=space_sniffer_everything
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/

# du -h EffectorP-3.0/ >> everything_du_03_05_24.txt
# du -h SURVIVOR/ >> everything_du_03_05_24.txt
# du -h SeqPrep/ >> everything_du_03_05_24.txt
# du -h TMHMM_29306/ >> everything_du_03_05_24.txt
# du -h art_bin_MountRainier/ >> everything_du_03_05_24.txt
# du -h bismark/ >> everything_du_03_05_24.txt
# du -h ciri2/ >> everything_du_03_05_24.txt
# du -h citrus_rnaseq/ >> everything_du_03_05_24.txt
# du -h conda_env/ >> everything_du_03_05_24.txt
# du -h conda_envs/ >> everything_du_03_05_24.txt
# du -h delly/ >> everything_du_03_05_24.txt
# du -h diamond/ >> everything_du_03_05_24.txt
# du -h ecc_finder/ >> everything_du_03_05_24.txt
# du -h effector_prediction/ >> everything_du_03_05_24.txt
# du -h foldseek/ >> everything_du_03_05_24.txt
# du -h fungap_runs/ >> everything_du_03_05_24.txt
# du -h genome/ >> everything_du_03_05_24.txt
# du -h genomes_and_plant_health_summer_2021/ >> everything_du_03_05_24.txt
# du -h lumpy-sv/ >> everything_du_03_05_24.txt
# du -h lumpy-sv_old/ >> everything_du_03_05_24.txt
# du -h magnaporthe_eccdnas/ >> everything_du_03_05_24.txt
# du -h manta-1.6.0.centos6_x86_64/ >> everything_du_03_05_24.txt
# du -h methylpy/ >> everything_du_03_05_24.txt
# du -h minimap2/ >> everything_du_03_05_24.txt
# du -h moryzae_virulence_code/ >> everything_du_03_05_24.txt
# du -h moryzae_virulence_project/ >> everything_du_03_05_24.txt
# du -h mummer_4/ >> everything_du_03_05_24.txt
# du -h mummer_runs/ >> everything_du_03_05_24.txt
# du -h new_guy11_reference/ >> everything_du_03_05_24.txt
# du -h pannzer2/ >> everything_du_03_05_24.txt
# du -h parse_pfam_domains/ >> everything_du_03_05_24.txt
# du -h pfam/ >> everything_du_03_05_24.txt
# du -h plot_profiler/ >> everything_du_03_05_24.txt
# du -h primer3/ >> everything_du_03_05_24.txt
# du -h ragoo_output/ >> everything_du_03_05_24.txt
# du -h ragoo_runs/ >> everything_du_03_05_24.txt
# du -h random_forest_memory_test/ >> everything_du_03_05_24.txt
# du -h raxml_ng/ >> everything_du_03_05_24.txt
# du -h raxml_ng_savio1/ >> everything_du_03_05_24.txt
# du -h references/ >> everything_du_03_05_24.txt
# du -h repbase/ >> everything_du_03_05_24.txt
# du -h repeatmasker/ >> everything_du_03_05_24.txt
# du -h rnaseq/ >> everything_du_03_05_24.txt
# du -h scooping_emile/ >> everything_du_03_05_24.txt
# du -h script_tests/ >> everything_du_03_05_24.txt
# du -h scripts/ >> everything_du_03_05_24.txt
# du -h signalp-4.1/ >> everything_du_03_05_24.txt
# du -h software/ >> everything_du_03_05_24.txt
# du -h space_snif/ >> everything_du_03_05_24.txt
# du -h spaced_anchor/ >> everything_du_03_05_24.txt
# du -h spaced_anchor_runs/ >> everything_du_03_05_24.txt
# du -h sra/ >> everything_du_03_05_24.txt
# du -h sv-callers-old/ >> everything_du_03_05_24.txt
# du -h sv-callers/ >> everything_du_03_05_24.txt
# du -h sv-callers_old/ >> everything_du_03_05_24.txt
# du -h sv-callers_runs/ >> everything_du_03_05_24.txt
# du -h sv_calling_moryzae/ >> everything_du_03_05_24.txt
# du -h sv_calling_workshop/ >> everything_du_03_05_24.txt
# du -h syri/ >> everything_du_03_05_24.txt
# du -h tmhmm-2.0c/ >> everything_du_03_05_24.txt
# du -h tmhmm_runs/ >> everything_du_03_05_24.txt
# du -h tmp_signalp/ >> everything_du_03_05_24.txt
# du -h wham/ >> everything_du_03_05_24.txt


rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/sv_calling_moryzae
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/sra
rsync -a --delete /global/scratch/users/pierrj/PAV_SV/PAV/empty_dir/ /global/scratch/users/pierrj/scooping_emile/

tar -czvf /global/scratch/projects/vector_kvklab/pierrj/everything_else.tar.gz --exclude='PAV_SV' --exclude='moryzae_virulence_project' /global/scratch/users/pierrj/