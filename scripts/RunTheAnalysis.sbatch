#!/usr/bin/bash
#SBATCH -J "MethylomeRun"
#SBATCH -A "MethylomeRun"
#SBATCH -c 16				   ### number of files to treat 
#SBATCH --mem=30G 				   ### blastn purpose
#SBATCH --mail-user="lohmanne@students.iarc.fr"
#SBATCH --mail-type=ALL
#SBATCH	--output=sbatch/%J-%j.out                    ### %j will be replaced by the batch job id
#SBATCH --error=sbatch/%J-%j.err                     

uname -a
echo starting job
date

source activate qiime2-2020.8


RunTheAnalysis.sh 

