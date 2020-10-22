#!/usr/bin/bash 
##########################################mise en page#########################################
#
set -e #stop the current script if an error occur
trap 'kill $(jobs -p)' EXIT

source apparence.sh 

<< ////

Tips:
The .qza and .qzv files can be unpacked using the unix command unzip or the qiime commands qiime tools extract or qiime tools export.

QIIME2 uses two different file types that contain the data and metadata from an analysis: .qza files are data files while .qzv files are visualizations.

The authors of QIIME2 call these data files “data artifacts” to indicate that they are objects containing data and metadata about an experiment.

The raw data in these files can be accessed using the command qiime tools export.

////

print_centered "QIIME2 ANALYSIS rerun"
print_centered "don't forget to activate conda environment"
print_centered "Development under qiime2-2019.10"
date 
start=$SECONDS

DIR=$(pwd)
Manifest=${DIR}"/files/manifest/Manifest"
MetaNames=${DIR}"/files/metada/Metadata"
OUT=${DIR}"/output"
LOG=${DIR}"/log"
OTU_database=${DIR}"/database/greengenes"
threads=8
pathOUT=${DIR}"/pathabun_core_metrics_out"
picrustOUT=${DIR}"/q2-picrust2_output"

# we have always used the same primer to amplify the V3-V4 sequences of the 16sRNA gene
forwardcut="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG "
reversecut="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

fcutadap=($forwardcut)
rcutadap=($reversecut)


#########################################   BEGINNING   #########################################
#########################################       0       #########################################
print_centered "0 - setting up of the working environment..."

if [ ! -d ${OUT} ]; then
  mkdir ${OUT}
fi
if [ ! -d ${OTU_database} ]; then
  echo "There is no database for the analysis to go through"
fi
if [ ! -d ${LOG} ]; then
  mkdir ${LOG}
fi

print_centered " Everything is well set ! "
print_centered "."
print_centered "."
print_centered "."

#########################################       1       #########################################

print_centered "1 - qiime tools import"

<< ////
Usage:
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/manifest_CRC_OB.txt --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
Utility:   
tools               Tools for working with QIIME 2 files.
////

for el in 1 2 3 4; do
if [ ! -e ${OUT}/${el}/paired-end-demux${el}.qza ]; then
    begin=$SECONDS
    echo ${el} " run processing..." ; 

    exe qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path ${Manifest}${el}".txt" \
     --output-path ${OUT}/${el}/paired-end-demux${el}.qza --input-format PairedEndFastqManifestPhred33V2 ;

     echo ${el} " DONE"
else echo ${el} " ALREADY DONE"
fi 
done
wait

########################################       2       #########################################
print_centered "2 - qiime cutadapt trim-paired on demultiplexed data"

<< ////
Usage:
qiime cutadapt trim-paired --i-demultiplexed-sequences paired-end-demux.qza --p-cores 15 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --o-trimmed-sequences paired-end-demux.trim.qza
Utility: 
  cutadapt            Plugin for removing adapter sequences, primers, and
                      other unwanted sequence from sequence data.
////

for el in 1 2 3 4; do
if [ ! -e ${OUT}/${el}/paired-end-demux${el}_trim.qza ]; then
    begin=$SECONDS; 

    exe qiime cutadapt trim-paired \
        --i-demultiplexed-sequences ${OUT}/${el}/paired-end-demux${el}.qza \
        --p-cores ${threads} \
        --p-front-f ${fcutadap[${el}]} \
        --p-front-r ${rcutadap[${el}]} \
        --o-trimmed-sequences ${OUT}/${el}/paired-end-demux${el}_trim.qza ;
     echo ${el} " DONE"
else echo ${el} " ALREADY DONE"
fi 
done
wait
#########################################       3       #########################################
print_centered "3 - qiime dada2 denoise-paired"

<< ////
Usage:
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.trim.qza --p-trim-left-f 20 --p-trim-left-r 0 --p-trunc-len-f 290 --p-trunc-len-r 260 --p-no-hashed-feature-ids --o-representative-sequences okrep-seqs-dada2.qza --o-denoising-stats okden-seqs-dada2.qza --o-table oktable-dada2.qza --verbose
This method denoises paired-end sequences, dereplicates them, and filters
  chimeras.
Utility: 
  dada2               Plugin for sequence quality control with DADA2.
  denoise-paired  Denoise and dereplicate paired-end sequences
Inputs:
  --i-demultiplexed-seqs ARTIFACT SampleData[PairedEndSequencesWithQuality]
                         The paired-end demultiplexed sequences to be
                         denoised.

Outputs: 
--o-representative-sequences ARTIFACT FeatureData[Sequence]
                         The resulting feature sequences. Each feature in the
                         feature table will be represented by exactly one
                         sequence, and these sequences will be the joined
                         paired-end sequences         

--o-denoising-stats ARTIFACT SampleData[DADA2Stats]

--o-table ARTIFACT FeatureTable[Frequency]
                         The resulting feature table.

Good to know : 
     ${OUT}/okrep-seqs-dada2.qza will be FeatureData[Sequence]
     ${OUT}/okden-seqs-dada2.qza will be SampleData[DADA2Stats]
     ${OUT}/oktable-dada2.qza will be FeatureTable[Frequency]
////

for el in 1 2 3 4; do
if [ ! -e ${OUT}/${el}/oktable-dada2_${el}.qza ]; then
    begin=$SECONDS ;
exe qiime dada2 denoise-paired \
     --i-demultiplexed-seqs ${OUT}/${el}/paired-end-demux${el}_trim.qza \
     --p-trim-left-f 20 \
     --p-trim-left-r 0 \
     --p-trunc-len-f 290 \
     --p-trunc-len-r 260 \
     --p-no-hashed-feature-ids \
     --p-n-threads ${threads} \
     --o-representative-sequences ${OUT}/${el}/okrep-seqs-dada2_${el}.qza \
     --o-denoising-stats ${OUT}/${el}/okden-seqs-dada2_${el}.qza \
     --o-table ${OUT}/${el}/oktable-dada2_${el}.qza --verbose  > ${LOG}/okall-dada2_${el}.log ;
     echo ${el} " DONE"
else echo ${el} " ALREADY DONE"
fi 
done
wait


# #########################################       4       #########################################
# print_centered "4 - qiime feature-table filter-samples"

# << ////
# Usage:
# qiime feature-table filter-samples --i-table oktable-dada2.qza --m-metadata-file MetadataRun2_OB_F.txt --o-filtered-table Run2_OB_F-table-dada2.qza
# Utility: 
#   feature-table       Plugin for working with sample by feature tables.
#   filter-samples      Filter samples from table
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency¹ | RelativeFrequency² |
#     PresenceAbsence³ | Composition⁴]
#                        The feature table from which samples should be
#                        filtered.
# Output:
#   --o-filtered-table ARTIFACT FeatureTable[Frequency¹ | RelativeFrequency²
#     | PresenceAbsence³ | Composition⁴]
#                        The resulting feature table filtered by sample.

# ////

# if [ ! -e ${OUT}/Run2_OB_F-table-dada2.qza ]; then
#      qiime feature-table filter-samples \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --m-metadata-file ${MetaNames} \
#      --o-filtered-table ${OUT}/Run2_OB_F-table-dada2.qza
#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# #########################################       5       #########################################
# print_centered "5 - qiime feature-table filter-seqs"

# << ////
# Usage:
# qiime feature-table filter-seqs --i-data okrep-seqs-dada2.qza --i-table Run2_OB_F-table-dada2.qza --o-filtered-data Run2_OB_F-rep-seqs-dada2.qza
# Utility: 
#   feature-table       Plugin for working with sample by feature tables.
#   filter-seqs         Filter features from sequences
# ////

# if [ ! -e ${OUT}/Run2_OB_F-rep-seqs-dada2.qza ]; then
# qiime feature-table filter-seqs \
#      --i-data ${OUT}/okrep-seqs-dada2.qza \
#      --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
#      --o-filtered-data ${OUT}/Run2_OB_F-rep-seqs-dada2.qza
#         echo "DONE"
# else echo "ALREADY DONE"
# fi


# #########################################       6       #########################################
# print_centered "6 - qiime tools import otus fasta"

# << ////
# Usage:
# qiime tools import --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otus.fasta --output-path 99_otus.qza --type 'FeatureData[Sequence]'
# Utility: 
#   tools               Tools for working with QIIME 2 files.
#   import            Import data into a new QIIME 2 Artifact.
#   --input-path PATH       Path to file or directory that should be imported.
#   --output-path ARTIFACT  Path where output artifact should be written.
#   --type TEXT             The semantic type of the artifact that will be
#                           created upon importing.
# ////
# if [ ! -e ${OUT}/99_otus.qza ]; then
# qiime tools import \
#      --input-path ${OTU_database}/99_otus.fasta \
#      --output-path ${OUT}/99_otus.qza \
#      --type 'FeatureData[Sequence]'
#         echo "DONE"
# else echo "ALREADY DONE"
# fi

# #########################################       7       #########################################
# print_centered "7 - qiime tools import taxonomy otus"

# << ////
# Usage:
# qiime tools import --type FeatureData[Taxonomy] --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otu_taxonomy.txt --input-format HeaderlessTSVTaxonomyFormat --output-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otu_taxonomy.qza
# Utility: 
#   tools               Tools for working with QIIME 2 files.
#   import            Import data into a new QIIME 2 Artifact.
#   --type TEXT             The semantic type of the artifact that will be
#                           created upon importing.
# ////
# if [ ! -e ${OUT}/99_otu_taxonomy.qza ]; then
# qiime tools import \
#      --type FeatureData[Taxonomy] \
#      --input-path ${OTU_database}/99_otu_taxonomy.txt \
#      --input-format HeaderlessTSVTaxonomyFormat \
#      --output-path ${OUT}/99_otu_taxonomy.qza
#         echo "DONE"
# else echo "ALREADY DONE"
# fi

# #########################################                #########################################
# print_centered "8 - qiime feature-classifier extract-reads "
# << ////
# Usage:
# qiime feature-classifier extract-reads --i-sequences 99_otus.qza --p-f-primer CCTACGGGNGGCWGCAG --p-r-primer GACTACHVGGGTATCTAATCC --o-reads 99_otus-ref.seqs.qza

# Utility: 
#   Extract simulated amplicon reads from a reference database.
# Inputs:
#   --i-sequences ARTIFACT FeatureData[Sequence]
# Outputs:
#   --o-reads ARTIFACT FeatureData[Sequence]
# To Do:
# ////

# if [ ! -e ${OUT}/99_otus-ref.seqs.qza ]; then
#      begin=$SECONDS
#      qiime feature-classifier extract-reads \
#      --i-sequences ${OUT}/99_otus.qza \
#      --p-f-primer CCTACGGGNGGCWGCAG \
#      --p-r-primer GACTACHVGGGTATCTAATCC \
#      --o-reads ${OUT}/99_otus-ref.seqs.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# #########################################               #########################################
# print_centered "9 - qiime feature-classifier fit-classifier-naive-bayes"
# << ////
# Usage:
# qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads 99_otus-ref.seqs.qza --i-reference-taxonomy 99_otu_taxonomy.qza --o-classifier classifier.trained.qza

# Utility: 
#   Create a scikit-learn naive_bayes classifier for reads
# Input:
#   --i-reference-reads ARTIFACT FeatureData[Sequence]
#   --i-reference-taxonomy ARTIFACT FeatureData[Taxonomy]
# To Do:
# ////

# if [ ! -e ${OUT}/classifier.trained.qza ]; then
#      begin=$SECONDS
#      qiime feature-classifier fit-classifier-naive-bayes \
#      --i-reference-reads ${OUT}/99_otus-ref.seqs.qza \
#      --i-reference-taxonomy ${OUT}/99_otu_taxonomy.qza \
#      --o-classifier ${OUT}/classifier.trained.qza
#      echo "DONE"
# else echo "ALREADY DONE"
# fi

# #########################################               #########################################
# print_centered "10 - qiime feature-classifier classify-sklearn"
# << ////
# Usage:
# qiime feature-classifier classify-sklearn --i-classifier classifier.trained.qza --i-reads newCpositive-rep-seqs-dada2.qza --o-classification newtaxonomy.sklearn.qza 

# Utility: 
#   Classify reads by taxon using a fitted classifier.

# To Do:
# Notice: I change newCpositive-rep-seqs-dada2.qza to Run2_OB_F-rep-seqs-dada2.qza for tests
# ////
# if [ ! -e ${OUT}/newtaxonomy.sklearn.qza  ]; then
#      begin=$SECONDS
#      qiime feature-classifier classify-sklearn \
#      --i-classifier ${OUT}/classifier.trained.qza \
#      --i-reads ${OUT}/okrep-seqs-dada2.qza \
#      --o-classification ${OUT}/newtaxonomy.sklearn.qza 
#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# # #########################################               #########################################
# print_centered "11  - qiime taxa collapse "
# << ////
# Usage:
# qiime taxa collapse --i-table table-dada2.qza --i-taxonomy taxonomy.sklearn.qza --p-level 2 --o-collapsed-table FOB_freq-2.qza

# Utility: 
#   Collapse groups of features that have the same taxonomic assignment
#   through the specified level. The frequencies of all features will be
#   summed when they are collapsed.

# To Do:
# ////
# if [ ! -e ${OUT}/FOB_freq-2.qza ]; then
#      begin=$SECONDS
#      qiime taxa collapse \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --i-taxonomy ${OUT}/newtaxonomy.sklearn.qza  \
#      --p-level 2 \
#      --o-collapsed-table ${OUT}/FOB_freq-2.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# # # ########################################               #########################################
#  print_centered " 12 - qiime feature-table relative-frequency"
# << ////
# Usage:
# qiime feature-table relative-frequency --i-table FOB_freq-2.qza --o-relative-frequency-table relative-freq-2.qza
# Utility: 
# Convert frequencies to relative frequencies by dividing each frequency in
#   a sample by the sum of frequencies in that sample.

# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency]
#                        The feature table to be converted into relative
#                        frequencies.
# Outputs:
#   --o-relative-frequency-table ARTIFACT FeatureTable[RelativeFrequency]
#                        The resulting relative frequency feature table.
# To Do:
# ////

# if [ ! -e ${OUT}/relative-freq-2.qza ]; then
#      begin=$SECONDS
#      qiime feature-table relative-frequency \
#      --i-table ${OUT}/FOB_freq-2.qza \
#      --o-relative-frequency-table ${OUT}/relative-freq-2.qza \
#      --verbose > ${LOG}/qiime-feature-table-relative-frequency.log
#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# #########################################            #########################################
# print_centered "13 - qiime alignment mafft "
# << ////
# Usage:
# qiime alignment mafft --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza

# Utility: 
#   Perform de novo multiple sequence alignment using MAFFT.
# Input:
#   --i-sequences ARTIFACT FeatureData[Sequence]
#                           The sequences to be aligned.
# Output:
#   --o-alignment ARTIFACT FeatureData[AlignedSequence]

# To Do:
# ////
# if [ ! -e ${OUT}/aligned-rep-seqs.qza ]; then
#      begin=$SECONDS
#      qiime alignment mafft \
#      --i-sequences ${OUT}/okrep-seqs-dada2.qza  \
#      --o-alignment ${OUT}/aligned-rep-seqs.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered "14 - qiime alignment mask"
# << ////
# Usage:
# qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

# Utility: 
#  Mask (i.e., filter) unconserved and highly gapped columns from an
#   alignment. Default min_conservation was chosen to reproduce the mask
#   presented in Lane (1991)

# To Do:
# ////
# if [ ! -e ${OUT}/masked-aligned-rep-seqs.qza ]; then
#      begin=$SECONDS
#      qiime alignment mask \
#      --i-alignment ${OUT}/aligned-rep-seqs.qza \
#      --o-masked-alignment ${OUT}/masked-aligned-rep-seqs.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered "15 - qiime phylogeny fasttree "
# << ////
# Usage:
# qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

# Utility: 
#   Construct a phylogenetic tree with FastTree.
# Inputs:
#   --i-alignment ARTIFACT FeatureData[AlignedSequence]
#                           Aligned sequences to be used for phylogenetic
#                           reconstruction.   

# To Do:
# ////
# if [ ! -e ${OUT}/unrooted-tree.qza ]; then
#      begin=$SECONDS
#      qiime phylogeny fasttree \
#      --i-alignment ${OUT}/masked-aligned-rep-seqs.qza \
#      --o-tree ${OUT}/unrooted-tree.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered "16 - qiime phylogeny midpoint-root"
# << ////
# Usage:
# qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza 

# Utility: 
#   Midpoint root an unrooted phylogenetic tree.
# Inputs:
#   --i-tree ARTIFACT Phylogeny[Unrooted]
#                        The phylogenetic tree to be rooted.          
# Outputs:
#   --o-rooted-tree ARTIFACT
#     Phylogeny[Rooted]  The rooted phylogenetic tree.   

# To Do:
# ////
# if [ ! -e ${OUT}/rooted-tree.qza  ]; then
#      begin=$SECONDS
#      qiime phylogeny midpoint-root \
#      --i-tree ${OUT}/unrooted-tree.qza \
#      --o-rooted-tree ${OUT}/rooted-tree.qza 

#      echo "DONE"
# else echo "ALREADY DONE"
# fi

# # #########################################        BUG      #########################################
# # print_centered " - qiime diversity alpha-rarefaction "
# # << ////
# # Usage:
# # #qiime diversity alpha-rarefaction --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-max-depth 27828 --p-steps 75 --p-iterations 55 --m-metadata-file Metadata_FOB_BS.tsv --p-metrics chao1 --p-metrics simpson_e --p-metrics simpson --p-metrics shannon --p-metrics observed_otus --p-metrics faith_pd --o-visualization rarefaction-curve.qzv
# # Utility: 
# #   Generate interactive alpha rarefaction curves by computing rarefactions
# #   between min_depth and max_depth. The number of intermediate depths to
# #   compute is controlled by the steps parameter, with n iterations being
# #   computed at each rarefaction depth. If sample metadata is provided,
# #   samples may be grouped based on distinct values within a metadata column.
# # Inputs:
# #   --i-table ARTIFACT FeatureTable[Frequency]

# # Outputs:
# # Outputs:
# #   --o-visualization VISUALIZATION

# # To Do:
# # ////
# # if [ ! -e ${OUT}/rarefaction-curve.qzv ]; then
# #      begin=$SECONDS
# #      qiime diversity alpha-rarefaction \
# #      --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
# #      --i-phylogeny ${OUT}/rooted-tree.qza \
# #      --p-iterations 55 \
# #      --m-metadata-file ${MetaNames} \
# #      --p-max-depth 27828 \
# #      --p-steps 75 \
# #      --p-metrics chao1 \
# #      --p-metrics simpson_e \
# #      --p-metrics simpson \
# #      --p-metrics shannon \
# #      --p-metrics observed_otus \
# #      --p-metrics faith_pd \
# #      --o-visualization ${OUT}/rarefaction-curve.qzv

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi
# # #########################################       seems to work        #########################################
# print_centered "17 - qiime diversity beta"
# << ////
# Usage:
# qiime diversity beta --i-table table-dada2.qza --p-metric braycurtis --o-distance-matrix dada2.braycurtis.notNorm.diversity.qza
# Utility: 
#   Computes a user-specified beta diversity metric for all pairs of samples
#   in a feature table.
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency | RelativeFrequency |
#     PresenceAbsence]   The feature table containing the samples over which
#                        beta diversity should be computed. 
# Outputs:
#   --o-distance-matrix ARTIFACT
#     DistanceMatrix     The resulting distance matrix.

# To Do:
# ////
# if [ ! -e ${OUT}/dada2.braycurtis.notNorm.diversity.qza ]; then
#      begin=$SECONDS
#      qiime diversity beta \
#      --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
#      --p-metric braycurtis \
#      --o-distance-matrix ${OUT}/dada2.braycurtis.notNorm.diversity.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi

# #########################################               #########################################
# print_centered "18 - qiime diversity pcoa"
# << ////
# Usage:
# qiime diversity pcoa --i-distance-matrix dada2.braycurtis.notNorm.diversity.qza --o-pcoa dada2.braycurtis.notNorm.diversity.pcoa.qza

# Utility: 
#   Apply principal coordinate analysis.

# To Do:
# ////
# if [ ! -e ${OUT}/dada2.braycurtis.notNorm.diversity.pcoa.qza ]; then
#      begin=$SECONDS
#      qiime diversity pcoa \
#      --i-distance-matrix ${OUT}/dada2.braycurtis.notNorm.diversity.qza \
#      --o-pcoa ${OUT}/dada2.braycurtis.notNorm.diversity.pcoa.qza
     
#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# # #########################################       BUG        #########################################
# # print_centered " - qiime diversity core-metrics-phylogenetic"
# # << ////
# # Usage:
# # qiime diversity core-metrics-phylogenetic --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-sampling-depth 5801 --output-dir dada2-diversity-5801 --m-metadata-file Metadata_FOB_BS.tsv

# # Utility: 
# #  Applies a collection of diversity metrics (both phylogenetic and non-
# #   phylogenetic) to a feature table.
# # Inputs:
# #   --i-table ARTIFACT FeatureTable[Frequency]
# #                           The feature table containing the samples over which
# #                           diversity metrics should be computed.    
# #   --i-phylogeny ARTIFACT  Phylogenetic tree containing tip identifiers that
# #     Phylogeny[Rooted]     correspond to the feature identifiers in the table.
# #                           This tree can contain tip ids that are not present
# #                           in the table, but all feature ids in the table must
# #                           be present in this tree.                  

# # To Do:
# # ////
# # if [ ! -d diversity ]; then
# #      mkdir ${DIR}/diversity
# #      begin=$SECONDS
# #      qiime diversity core-metrics-phylogenetic \
# #      --i-table ${OUT}/oktable-dada2.qza \
# #      --i-phylogeny r ${OUT}/rooted-tree.qza \
# #      --p-sampling-depth 5801 \
# #      --output-dir ${DIR}/diversity \
# #      --m-metadata-file ${MetaNames}
# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi

# #########################################               #########################################
# print_centered "19 - qiime gneiss correlation-clustering"
# << ////
# Usage:
# qiime gneiss correlation-clustering --i-table table-dada2.qza --p-pseudocount 1 --o-clustering hierarchy1.qza         

# Utility: 
#   Build a bifurcating tree that represents a hierarchical clustering of
#   features.  The hiearchical clustering uses Ward hierarchical clustering
#   based on the degree of proportionality between features.
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency]
#                           The feature table containing the samples in which
#                           the columns will be clustered.  

# To Do:
# ////
# if [ ! -e ${OUT}/hierarchy1.qza  ]; then
#      begin=$SECONDS
#      qiime gneiss correlation-clustering \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --p-pseudocount 1 \
#      --o-clustering ${OUT}/hierarchy1.qza         

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# # #########################################               #########################################
# # print_centered " - "
# # << ////
# # Usage:
# # qiime gneiss gradient-clustering --i-table table-dada2.qza --m-gradient-file Metadata_FOB_BS.tsv --m-gradient-column TimePoint --o-clustering gradient-hierarchy.qza 

# # Utility: 

# # To Do:
# # ////
# # if [ ! -e ${OUT}/... ]; then
# #      begin=$SECONDS
# #      qiime gneiss gradient-clustering --i-table table-dada2.qza --m-gradient-file Metadata_FOB_BS.tsv --m-gradient-column TimePoint --o-clustering gradient-hierarchy.qza 

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi

# #########################################               #########################################
# print_centered "20 - qiime gneiss ilr-hierarchical"
# << ////
# Usage:
# qiime gneiss ilr-hierarchical --i-table table-dada2.qza --i-tree hierarchy1.qza --o-balances balances.qza

# Utility: 
#   Calculate balances given a hierarchy.
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency | Composition]
#                           The feature table containing the samples in which
#                           the ilr transform will be performed.      
#   --i-tree ARTIFACT       A hierarchy of feature identifiers that defines the
#     Hierarchy             partitions of features.  Each tip in the
#                           hierarchycorresponds to the feature identifiers in
#                           the table. This tree can contain tip ids that are
#                           not present in the table, but all feature ids in the
#                           table must be present in this tree.  This assumes
#                           that all of the internal nodes in the tree have
#                           labels.                                   

# Outputs:
#   --o-balances ARTIFACT FeatureTable[Balance]
#                           The resulting balances from the ilr transform.
                                                                    

# To Do:
# ////
# if [ ! -e ${OUT}/balances.qza ]; then
#      begin=$SECONDS
#      qiime gneiss ilr-hierarchical \
#      --i-table ${OUT}/oktable-dada2.qza  \
#      --i-tree ${OUT}/hierarchy1.qza \
#      --o-balances ${OUT}/balances.qza

#      echo "DONE"
# else echo "ALREADY DONE"
# fi


# # #########################################    feature table biom that I don't have            #########################################
# # print_centered " - "
# # << ////
# # Usage:
# # qiime tools import --input-path feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table2.qza    

# # Utility: 

# # To Do:
# # ////
# # if [ ! -e ${OUT}/... ]; then
# #      begin=$SECONDS

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi
# #########################################               #########################################
# print_centered "21 - qiime picrust2 full-pipeline"
# << ////
# beware: not compatible with qiime2-2020.8 but with 2019.10 due to some dependancies
# wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml
# conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
# # OPTIONAL CLEANUP
# rm qiime2-2019.10-py36-linux-conda.yml

# DO : conda activate qiime2-2019.10 
# and install q2-picrust2 with conda
# Instruction :  conda install -c gavinmdouglas q2-picrust2 
# Usage:

# qiime picrust2 full-pipeline --i-table table-dada2.qza --i-seq rep-seqs-dada2.qza --output-dir q2-picrust2_output --p-threads 1 --p-hsp-method pic --p-max-nsti 2 --verbose

# Utility: 
#   QIIME 2 plugin for default 16S PICRUSt2 pipeline

# To Do: need to install qiime2 https://forum.qiime2.org/t/qiime2-modulenotfounderror-with-picrust2-plugin-installation/8985/2 
# ////

# if [ ! -d ${picrustOUT} ]; then
#      begin=$SECONDS
#      qiime picrust2 full-pipeline \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --i-seq ${OUT}/okrep-seqs-dada2.qza \
#      --output-dir ${picrustOUT} \
#      --p-threads ${threads} \
#      --p-hsp-method pic \
#      --p-max-nsti 2 --verbose

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered "22 -  BUG - qiime diversity core-metrics "
# # << ////
# # Usage:
# # qiime diversity core-metrics --i-table q2-picrust2_output/pathway_abundance.qza --p-sampling-depth 	548854 --m-metadata-file Metadata_FOB_BS.tsv --output-dir pathabun_core_metrics_out --p-n-jobs 1

# # Inputs:
# #   --i-table ARTIFACT FeatureTable[Frequency]
# #                        The feature table containing the samples over which
# #                        diversity metrics should be computed.

# # Utility: 
# #       Applies a collection of diversity metrics (non-phylogenetic) to a feature
# #         table.
  
# # To Do:
# # Supposition : ${OUT}/pathway_abundance.qza must have been created by the picrust2 command. 
# # ////
# # if [ ! -d ${pathOUT} ]; then
# #      begin=$SECONDS
# #     qiime diversity core-metrics \
# #       --i-table ${picrustOUT}/pathway_abundance.qza \
# #       --p-sampling-depth 	548854 \
# #       --m-metadata-file ${MetaNames} \
# #       --output-dir ${pathOUT} 
# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi
# # #########################################               #########################################
# print_centered "23 - BUG - qiime tools export"
# # << ////
# # Usage:
# # qiime tools export --input-path q2-picrust2_output/pathway_abundance.qza --output-path pathabun_exported 

# # Utility: 

# # To Do:
# # ////
# # if [ ! -d ${picrustOUT}/export ]; then
# #     mkdir ${picrustOUT}/export
# #     begin=$SECONDS
# #     qiime tools export \
# #       --input-path ${picrustOUT}/pathway_abundance.qza \
# #       --output-path ${picrustOUT}/export

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi

# #########################################               #########################################
# print_centered "24 - result BUG -> can't - biom convert pathway abundance picrust2 output"
# # << ////
# # Usage:
# # biom convert -i pathabun_exported/feature-table.biom -o pathabun_exported/feature-table.biom.tsv --to-tsv

# # Utility: 

# # To Do:
# # ////
# # if [ ! -e ${picrustOUT}/export/feature-table.biom.tsv ]; then
# #      begin=$SECONDS
# #       biom convert -i ${picrustOUT}/export/feature-table.biom \
# #       -o ${picrustOUT}/export/feature-table.biom.tsv --to-tsv

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi
# ########################################               #########################################
# print_centered "25 - qiime feature-table relative-frequency "
# << ////
# Usage:
# qiime feature-table relative-frequency --i-table pathway_abundance.qza --o-relative-frequency-table pathway_rel-freq.qza   

# Utility: Convert frequencies to relative frequencies by dividing each frequency in
#   a sample by the sum of frequencies in that sample.

# To Do:
# ////
# if [ ! -e ${picrustOUT}/pathway_rel-freq.qza  ]; then
#      begin=$SECONDS
#     qiime feature-table relative-frequency \
#       --i-table ${picrustOUT}/pathway_abundance.qza \
#       --o-relative-frequency-table ${picrustOUT}/pathway_rel-freq.qza   

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered "26 - qiime tools export pathway_rel-freq.qza "
# << ////
# Usage:
# qiime tools export --input-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq.qza --output-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq

# Utility: 

# To Do:
# ////

# if [ ! -e ${picrustOUT}/pathway_rel-freq   ]; then
#      begin=$SECONDS
#     qiime tools export \
#     --input-path ${picrustOUT}/pathway_rel-freq.qza   \
#     --output-path ${picrustOUT}/pathway_rel-freq  

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# ########################################               #########################################
# print_centered "27 - biom convert  pathway ref freq to tsv"
# << ////
# Usage:
# biom convert -i feature-table.biom -o pathway_rel-freq_tsv.txt --to-tsv

# Utility: 

# To Do:
# ////
# if [ ! -e ${picrustOUT}/pathway_rel-freq_tsv.txt ]; then
#      begin=$SECONDS
#     biom convert \
#     -i ${picrustOUT}/pathway_rel-freq/feature-table.biom \
#     -o ${picrustOUT}/pathway_rel-freq_tsv.txt --to-tsv

#      echo "DONE"
# else echo "ALREADY DONE"
# fi









print_centered "X - end of the analysis"

date 
