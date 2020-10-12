#!/usr/bin/bash
##########################################mise en page##########################################
#
function print_centered {
     [[ $# == 0 ]] && return 1

     declare -i TERM_COLS="$(tput cols)"
     declare -i str_len="${#1}"
     [[ $str_len -ge $TERM_COLS ]] && {
          echo "$1";
          return 0;
     }

     declare -i filler_len="$(( (TERM_COLS - str_len) / 2 ))"
     [[ $# -ge 2 ]] && ch="${2:0:1}" || ch=" "
     filler=""
     for (( i = 0; i < filler_len; i++ )); do
          filler="${filler}${ch}"
     done

     printf "%s%s%s" "$filler" "$1" "$filler"
     [[ $(( (TERM_COLS - str_len) % 2 )) -ne 0 ]] && printf "%s" "${ch}"
     printf "\n"

     return 0
}
#
#################################################################################################
set -e 
<< ////

Tips:
The .qza and .qzv files can be unpacked using the unix command unzip or the qiime commands qiime tools extract or qiime tools export.

QIIME2 uses two different file types that contain the data and metadata from an analysis: .qza files are data files while .qzv files are visualizations.

The authors of QIIME2 call these data files “data artifacts” to indicate that they are objects containing data and metadata about an experiment.

The raw data in these files can be accessed using the command qiime tools export.

////


print_centered "QIIME2 ANALYSIS rerun"
print_centered "don't forget to activate conda environment"
print_centered "Development under qiime2-2020.8"
date 
start=$SECONDS

DIR=$(pwd)
Manifest=${DIR}"/Manifest-localisation-filessmall.txt"
MetaNames=${DIR}"/Metadata18_OB_F.txt"
OUT=${DIR}"/output"
LOG=${DIR}"/log"

#########################################   BEGINNING   #########################################
#########################################       1       #########################################

print_centered "1 - qiime tools import"
<< ////
Usage:
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/manifest_CRC_OB.txt --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
Utility:   
tools               Tools for working with QIIME 2 files.
////

if [ ! -e ${OUT}/paired-end-demux.qza ]; then

    begin=$SECONDS
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${Manifest} \
        --output-path ${OUT}/paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2 
        echo "DONE"
else echo "ALREADY DONE"
fi

#########################################       2       #########################################
print_centered "2 - qiime cutadapt trim-paired on demultiplexed data"

<< ////
Usage:
qiime cutadapt trim-paired --i-demultiplexed-sequences paired-end-demux.qza --p-cores 15 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --o-trimmed-sequences paired-end-demux.trim.qza
Utility: 
  cutadapt            Plugin for removing adapter sequences, primers, and
                      other unwanted sequence from sequence data.
////
if [ ! -e ${OUT}/paired-end-demux.trim.qza ]; then
    begin=$SECONDS
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${OUT}/paired-end-demux.qza \
        --p-cores 4 \
        --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
        --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
        --o-trimmed-sequences ${OUT}/paired-end-demux.trim.qza
        echo "DONE"
else echo "ALREADY DONE"
fi


#########################################       3       #########################################
print_centered "3 - qiime dada2 denoise-paired"

<< ////
Usage:
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.trim.qza --p-trim-left-f 20 --p-trim-left-r 0 --p-trunc-len-f 290 --p-trunc-len-r 260 --p-no-hashed-feature-ids --o-representative-sequences okrep-seqs-dada2.qza --o-denoising-stats okden-seqs-dada2.qza --o-table oktable-dada2.qza --verbose
Utility: 
  dada2               Plugin for sequence quality control with DADA2.
  denoise-paired  Denoise and dereplicate paired-end sequences
////

if [ ! -e ${OUT}/oktable-dada2.qza ]; then
    begin=$SECONDS
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${OUT}/paired-end-demux.trim.qza \
    --p-trim-left-f 20 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 290 \
    --p-trunc-len-r 260 \
    --p-no-hashed-feature-ids \
    --o-representative-sequences ${OUT}/okrep-seqs-dada2.qza \
    --o-denoising-stats ${OUT}/okden-seqs-dada2.qza \
    --o-table ${OUT}/oktable-dada2.qza --verbose > ${LOG}/okall-dada2.log
        echo "DONE"
else echo "ALREADY DONE"
fi


#########################################       4       #########################################
print_centered "4 - qiime feature-table filter-samples"

<< ////
Usage:
qiime feature-table filter-samples --i-table oktable-dada2.qza --m-metadata-file MetadataRun2_OB_F.txt --o-filtered-table Run2_OB_F-table-dada2.qza
Utility: 
  feature-table       Plugin for working with sample by feature tables.
  filter-samples      Filter samples from table
////

if [ ! -e ${OUT}/Run2_OB_F-table-dada2.qza ]; then
qiime feature-table filter-samples \
     --i-table ${OUT}/oktable-dada2.qza \
     --m-metadata-file ${MetaNames} \
     --o-filtered-table ${OUT}/Run2_OB_F-table-dada2.qza
        echo "DONE"
else echo "ALREADY DONE"
fi


#########################################       5       #########################################
print_centered "5 - qiime feature-table filter-seqs"

<< ////
Usage:
qiime feature-table filter-seqs --i-data okrep-seqs-dada2.qza --i-table Run2_OB_F-table-dada2.qza --o-filtered-data Run2_OB_F-rep-seqs-dada2.qza
Utility: 
  feature-table       Plugin for working with sample by feature tables.
  filter-seqs         Filter features from sequences
////

if [ ! -e ${OUT}/Run2_OB_F-rep-seqs-dada2.qza ]; then
qiime feature-table filter-seqs \
     --i-data ${OUT}/okrep-seqs-dada2.qza \
     --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
     --o-filtered-data ${OUT}/Run2_OB_F-rep-seqs-dada2.qza
        echo "DONE"
else echo "ALREADY DONE"
fi

#########################################       6       #########################################
print_centered "5 - qiime tools export the dada2 table to readable"

<< ////
Usage:
qiime tools export --input-path table-dada2.qza --output-path exported-table-dada2.tsv
Utility: 
  tools               Tools for working with QIIME 2 files.
  export            Export data from a QIIME 2 Artifact or a Visualization
////

qiime tools export --input-path table-dada2.qza \
     --output-path exported-table-dada2.tsv



#########################################       7       #########################################
#########################################       8       #########################################
#########################################       9       #########################################
#########################################      10       #########################################
#########################################      11       #########################################
#########################################      12       #########################################
#########################################      13       #########################################
#########################################      14       #########################################
#########################################      15       #########################################
#########################################      16       #########################################
#########################################      17       #########################################
#########################################      18       #########################################
#########################################      19       #########################################
#########################################      20       #########################################














print_centered "X - end of the analysis"
date 