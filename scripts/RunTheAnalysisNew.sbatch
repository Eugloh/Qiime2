#!/usr/bin/bash
#SBATCH -J "runall-MethylomeRun"
#SBATCH -A "runall-MethylomeRun"
#SBATCH -c 32			   ### number of files to treat 
#SBATCH --mem=80G 				   ### blastn purpose
#SBATCH --mail-user="lohmanne@students.iarc.fr"
#SBATCH --mail-type=ALL
#SBATCH	--output=sbatch/runall-%J.out                    ### %j will be replaced by the batch job id
#SBATCH --error=sbatch/runall-%J.err                     

uname -a
echo starting job
date

#source activate qiime2-2020.8
source activate qiime2-2019.10


/home/lohmanne/Qiime2/scripts/Runtheanalysis-dev.sh
