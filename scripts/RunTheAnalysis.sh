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

#########################################       old2       #########################################
# echo "assuming already demuplexed data"
# print_centered "2 - qiime demux summarize"

# << ////
# Usage:
# qiime demux summarize --i-data paired-end-demux.qza --o-visualization paired-end-demux.qzv
# Utility: 
# demux               Plugin for demultiplexing & viewing sequence quality.
# ////

# if [ ! -e ${OUT}/paired-end-demux.qzv ]; then
#      begin=$SECONDS
#      qiime demux summarize \
#           --i-data ${OUT}/paired-end-demux.qza \
#           --o-visualization ${OUT}/paired-end-demux.qzv
#      echo "DONE"
# else echo "ALREADY DONE"
# fi

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

qiime dada2 denoise-paired -\
    -i-demultiplexed-seqs ${OUT}/paired-end-demux.trim.qza \
    --p-trim-left-f 20 \
    --p-trim-left-r 0 \
    --p-trunc-len-f 290 \
    --p-trunc-len-r 260 \
    --p-no-hashed-feature-ids \
    --o-representative-sequences ${OUT}/okrep-seqs-dada2.qza \
    --o-denoising-stats ${OUT}/okden-seqs-dada2.qza \
    --o-table ${OUT}/oktable-dada2.qza --verbose > ${LOG}/okall-dada2.log


#########################################       4       #########################################
#########################################       5       #########################################
#########################################       6       #########################################
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