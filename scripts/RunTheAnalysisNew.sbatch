#!/usr/bin/bash
#SBATCH -J "runall-MethylomeRun"
#SBATCH -A "runall-MethylomeRun"
#SBATCH -c 16				   ### number of files to treat 
#SBATCH --mem=80G 				   ### blastn purpose
#SBATCH --mail-user="lohmanne@students.iarc.fr"
#SBATCH --mail-type=ALL
#SBATCH	--output=sbatch/runall-%J-%j.out                    ### %j will be replaced by the batch job id
#SBATCH --error=sbatch/runall-%J-%j.err                     

uname -a
echo starting job
date

#source activate qiime2-2020.8
source activate qiime2-2019.10


/home/lohmanne/Qiime2/scripts/RuntheanalysisNew.sh
