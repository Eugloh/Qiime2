# Script folder

.
├── apparence.sh
├── PATHWAY.R
├── README.md # current file 
├── Runtheanalysis-dev.sh
├── RunTheAnalysisNew.sbatch # github import
└── RuntheanalysisNew.sh # will be removed

## Details on scripts

### Runtheanalysis-dev.sh 

Bash script dedicaced to the Qiime2 analysis. The current version should be launched with the conda environment qiime2-2019.10 and will need further extention for the picrust2 run. 

#### beware: not compatible with qiime2-2020.8 but with 2019.10 due to some dependancies

```bash
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.10-py36-linux-conda.yml

# then DO : 
conda activate qiime2-2019.10 
# and install q2-picrust2 with conda

conda install -c gavinmdouglas q2-picrust2 
```

### apparence.sh

script dedicaced to output display. 

### PATHWAY.R

Picrust2 results analysis with ALDEx2 R package. 

R script, takes long time to run and a certain RAM capacity on your computer.
Once the script has been launch once you should retrieve the data opening the .RData files in the dedicaced folder ../Rdata.

### RunTheAnalysisNew.sbatch 

Slurm call 


