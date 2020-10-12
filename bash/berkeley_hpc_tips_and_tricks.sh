#!/bin/bash


module load imagemagick
display test.png

module load ${module} ## load a new module
module list ## list currently loaded modules
module purge ## remove all loaded modules
module avail ## see all available modules
module remove ${module} ## remove a specific module
echo "module load ${module}" >> ~/.bash_profile ## add loading module to bash profile so it is always loaded when you run a job or log in
emacs ~/.bash_profile ## access .bash_profile to edit, can use any text editor here

conda create --name=$ENV_NAME python=3.7 ## pick python version
source activate $ENV_NAME
conda install $package_name

export R_LIBS_USER=/global/scratch/users/yourusername/scripts/R_libs

install.packages('fields', lib = '/scratch/users/myusername/R')

export SLURM_ACCOUNT=fc_kvkallow
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT

git add $directory
git update-index --chmod +x $script ## make bash scripts executable which can be useful
git commit -m "useful commit message"
git push

git pull
sbatch ~/git/$directory/$script ## example

git clone https://github.com/${your_gitrepo}

git config credential.helper store
git push https://github.com/${your_gitrepo}
## type Username and Password here

srun -A fc_kvkallow -p savio --qos=savio_normal -t 30 --pty bash ## start an interactive job

sbatch $yourjob.slurm ## submit a job
sbatch -p savio2 $yourjob.slum ## submit a job on a specific node, can pass many arguments here and they all overwrite the SLURM header
sbatch --export=variable1=$variable1,variable2=$variable2,variable3=$variable3 $yourjob.slurm ## export variables to your job, useful for submitting multiple jobs with a loop

sinfo ## get info on all nodes to see which nodes are available
squeue -p savio ## check on a specific node to see how full it is

squeue -u pierrj ## check the status of your jobs
scancel $JOBID ## cancel a job
scancel -u pierrj ## cancel all jobs you submitted

sacct -j $JOBID --format=${FORMATOPTIONS} ## get information from a job, see SLURM sacct docs for list of options, MaxRSS gives memory usage which is useful

check_usage.sh -E -u yourusername ## check how many SUs you have used
check_usage.sh -E -a fc_kvkallow ## check how many SUs the lab has used





${SLURM_JOBID} ## get jobid from within the slurm script
${SLURM_NTASKS} ## get number of tasks from within slurm script, useful to get thread number for parallelization
