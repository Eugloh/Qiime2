# <div align="center">Qiime2 analysis</div>

Development under the conda environment qiime2-2019.10

```bash
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml --no-check-certificate
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.10-py36-linux-conda.yml
conda activate qiime2-2019.10 
## then install q2-picrust2 with conda
conda install -c gavinmdouglas q2-picrust2 
```

## Oral and Gut microbiota in obese patients subjected to bariatric surgery : A Longitudinal study 

### All in one script 
**1** : Importing and demultiplex data, summarize the results, and examing quality of the reads.

**2** : Quality controlling sequences and building Feature Table and Feature Data.

**3** : Summarizing Feature Table and Feature Data.

**4** : Assigning Taxonomy.

**5** : Generating a phylogenetic tree.

**6** : Analyzing Alpha and Beta diversities.


__Beware__ : To launch the main script on a beegfs dependant environment you might faces some troubles and have to change your bash_profile file in a way that the tmp don't create itself on beegfs scratch. 

### Detailled steps of analysis :

    1 - qiime tools import
    2 - qiime cutadapt trim-paired on demultiplexed data
    3 - qiime dada2 denoise-paired
    4 - qiime merge
    5 - qiime feature-table filter-samples
    6 - qiime feature-table filter-seqs
    7 - qiime tools import otus fasta
    8 - qiime tools import taxonomy otus
    9 - qiime feature-classifier extract-reads 
    10 - qiime feature-classifier fit-classifier-naive-bayes
    11 - qiime feature-classifier classify-sklearn
    12 - qiime taxa collapse 
    13 - qiime feature-table relative-frequency
    14 - qiime alignment mafft 
    15 - qiime alignment mask
    16 - qiime phylogeny fasttree 
    17 - qiime phylogeny midpoint-root
    18 - qiime diversity alpha-rarefaction 
    19 - qiime diversity beta
    20 - qiime diversity pcoa
    21 - qiime diversity core-metrics-phylogenetic
    22 - qiime gneiss correlation-clustering
    23 - qiime gneiss gradient-clustering
    24 - qiime gneiss ilr-hierarchical
    25 - qiime picrust2 full-pipeline
    26 - qiime diversity core-metrics 
    27 - qiime feature-table relative-frequency 

### Environment setup 
    Qiime2
      .
      ├── data
      │   ├── 1
      │   ├── 2
      │   ├── 3
      │   └── 4
      ├── database
      │   ├── C+
      │   └── greengenes
      ├── files
      │   ├── manifest
      │   │   └── old
      │   ├── metadata
      │   ├── method
      │   └── old
      ├── log
      ├── output
      │   ├── 1
      │   ├── 2
      │   ├── 3
      │   ├── 4
      │   └── merged
      │       ├── export
      │       ├── F
      │       │   ├── core_metrics
      │       │   ├── diversity
      │       │   ├── export_picrust2
      │       │   ├── path_rel-freq_picrust2
      │       │   └── picrust2
      │       ├── FOB
      │       │   ├── core_metrics
      │       │   ├── diversity
      │       │   ├── export_picrust2
      │       │   ├── path_rel-freq_picrust2
      │       │   └── picrust2
      │       ├── S
      │       │   ├── core_metrics
      │       │   ├── diversity
      │       │   ├── export_picrust2
      │       │   ├── path_rel-freq_picrust2
      │       │   └── picrust2
      │       └── SOB
      │           ├── core_metrics
      │           ├── diversity
      │           ├── export_picrust2
      │           ├── path_rel-freq_picrust2
      │           └── picrust2
      ├── sbatch
      └── scripts
