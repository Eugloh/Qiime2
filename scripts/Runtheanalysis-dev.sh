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
MetaNames=${DIR}"/files/metadata/Metadata"
OUT=${DIR}"/output"
LOG=${DIR}"/log"
MERGED=${OUT}"/merged"
OTU_database=${DIR}"/database/greengenes"
threads=8

# we have always used the same primer to amplify the V3-V4 sequences of the 16sRNA gene
forwardcut="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG "
reversecut="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

fcutadap=($forwardcut)
rcutadap=($reversecut)

echo -e 'QIIME2 pipeline steps :

Step1 :\t Importing and demultiplex data, summarize the results, and examing quality of the reads.\n
Step2 :\t Quality controlling sequences and building Feature Table and Feature Data.\n
Step3 :\t Summarizing Feature Table and Feature Data.\n
Step4 :\t Assigning Taxonomy.\n
Step5 :\t Generating a phylogenetic tree.\n
Step6 :\t Analyzing Alpha and Beta diversities.\n'


#########################################   BEGINNING   #########################################
#########################################               #########################################
print_centered "0 - setting up of the working environment..."

if [ ! -d ${OUT} ]; then
  mkdir ${OUT}
  mkdir ${OUT}/1
  mkdir ${OUT}/2
  mkdir ${OUT}/3
  mkdir ${OUT}/4
fi
if [ ! -d ${OTU_database} ]; then
  echo "There is no database for the analysis to go through"
fi
if [ ! -d ${LOG} ]; then
  mkdir ${LOG}
fi 

print_centered "Everything is well set !"


#########################################               #########################################

print_centered "1 - qiime tools import"

<< ////
Usage:
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_/qime2try/manifest_CRC_OB.txt --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
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

     echo -e ${el} "\trun DONE"
else echo -e ${el} "\trun ALREADY DONE"
fi 
done
wait
date
########################################               #########################################
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
     echo -e ${el} "\trun DONE"
else echo -e ${el} "\trun ALREADY DONE"
fi 
done
wait
date
#########################################               #########################################
print_centered "3 - qiime dada2 denoise-paired"

<< ////
Usage:
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.trim.qza --p-trim-left-f 20 --p-trim-left-r 0 --p-trunc-len-f 290 --p-trunc-len-r 260 --p-no-hashed-feature-ids --o-representative-sequences rep-seqs-dada2.qza --o-denoising-stats den-seqs-dada2.qza --o-table table-dada2.qza --verbose
This method denoises paired-end sequences, dereplicates them, and filters
  chimeras.
Utility: 
  dada2               Plugin for sequence quality control with DADA2.
  denoise-paired  Denoise and dereplicate paired-end sequences
Inputs:
  --i-demultiplexed-seqs ARTIFACT SampleData[PairedEndSequencesWithQuality]

Outputs: 
--o-representative-sequences ARTIFACT FeatureData[Sequence]
--o-denoising-stats ARTIFACT SampleData[DADA2Stats]
--o-table ARTIFACT FeatureTable[Frequency]

Good to know : 
     ${OUT}/rep-seqs-dada2.qza will be FeatureData[Sequence]
     ${OUT}/den-seqs-dada2.qza will be SampleData[DADA2Stats]
     ${OUT}/table-dada2.qza will be FeatureTable[Frequency]
////

for el in 1 2 3 4; do
if [ ! -e ${OUT}/${el}/table-dada2_${el}.qza ]; then
    begin=$SECONDS ;
exe qiime dada2 denoise-paired \
     --i-demultiplexed-seqs ${OUT}/${el}/paired-end-demux${el}_trim.qza \
     --p-trim-left-f 20 \
     --p-trim-left-r 0 \
     --p-trunc-len-f 290 \
     --p-trunc-len-r 260 \
     --p-no-hashed-feature-ids \
     --p-n-threads ${threads} \
     --o-representative-sequences ${OUT}/${el}/rep-seqs-dada2_${el}.qza \
     --o-denoising-stats ${OUT}/${el}/den-seqs-dada2_${el}.qza \
     --o-table ${OUT}/${el}/table-dada2_${el}.qza --verbose  > ${LOG}/all-dada2_${el}.log 2>&1 ;
     echo -e ${el} "\trun DONE"
else echo -e ${el} "\trun ALREADY DONE"
fi 
done
wait
date
###############################||| I need to merge my tables |||#################################
#########################################               #########################################
print_centered "4 - qiime merge"

if [ ! -d ${MERGED} ];then
  mkdir ${MERGED}
fi ;

if [ ! -e ${MERGED}/table-dada2_merged.qza ]; then
exe qiime feature-table merge \
  --i-tables ${OUT}/1/table-dada2_1.qza \
  --i-tables ${OUT}/2/table-dada2_2.qza \
  --i-tables ${OUT}/3/table-dada2_3.qza \
  --i-tables ${OUT}/4/table-dada2_4.qza \
  --o-merged-table ${MERGED}/table-dada2_merged.qza ;
  echo "DONE"
else echo "ALREADY DONE"
fi 

if [ ! -e ${MERGED}/run-rep-seqs_merged.qza ]; then
exe qiime feature-table merge-seqs \
  --i-data ${OUT}/1/rep-seqs-dada2_1.qza \
  --i-data ${OUT}/2/rep-seqs-dada2_2.qza \
  --i-data ${OUT}/3/rep-seqs-dada2_3.qza \
  --i-data ${OUT}/4/rep-seqs-dada2_4.qza \
  --o-merged-data ${MERGED}/run-rep-seqs_merged.qza ;
  echo "DONE"
else echo "ALREADY DONE"
fi 
date

# #########################################                #########################################
print_centered "5 - qiime feature-table filter-samples"

<< ////
Usage:
qiime feature-table filter-samples --i-table table-dada2.qza --m-metadata-file MetadataRun2_OB_F.txt --o-filtered-table Run2_OB_F-table-dada2.qza
Utility: 
  feature-table       Plugin for working with sample by feature tables.
  filter-samples      Filter samples from table
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency¹ | RelativeFrequency² |
    PresenceAbsence³ | Composition⁴]
                       The feature table from which samples should be
                       filtered.
Output:
  --o-filtered-table ARTIFACT FeatureTable[Frequency¹ | RelativeFrequency²
    | PresenceAbsence³ | Composition⁴]
                       The resulting feature table filtered by sample.

////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/run-rep-dada2_${el}.qza ]; then
    exe qiime feature-table filter-samples \
     --i-table ${MERGED}/table-dada2_merged.qza \
     --m-metadata-file ${MetaNames}${el}".csv" \
     --o-filtered-table ${MERGED}/${el}/run-rep-dada2_${el}.qza ;
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
#########################################                #########################################
print_centered "6 - qiime feature-table filter-seqs"

<< ////
Usage:
qiime feature-table filter-seqs --i-data rep-seqs-dada2.qza --i-table Run2_OB_F-table-dada2.qza --o-filtered-data Run2_OB_F-rep-seqs-dada2.qza
Utility: 
  feature-table       Plugin for working with sample by feature tables.
  filter-seqs         Filter features from sequences
////

for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/data-dada2_${el}_filtered.qza ]; then
exe qiime feature-table filter-seqs \
     --i-data ${MERGED}/run-rep-seqs_merged.qza \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --o-filtered-data ${MERGED}/${el}/data-dada2_${el}_filtered.qza ;
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
#########################################               #########################################
print_centered "7 - qiime tools import otus fasta"

<< ////
Usage:
qiime tools import --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_/qime2try/99_otus.fasta --output-path 99_otus.qza --type 'FeatureData[Sequence]'
Utility: 
  tools               Tools for working with QIIME 2 files.
  import            Import data into a new QIIME 2 Artifact.
  --input-path PATH       Path to file or directory that should be imported.
  --output-path ARTIFACT  Path where output artifact should be written.
  --type TEXT             The semantic type of the artifact that will be
                          created upon importing.
////
if [ ! -e ${OUT}/99_otus.qza ]; then
exe qiime tools import \
     --input-path ${OTU_database}/99_otus.fasta \
     --output-path ${OUT}/99_otus.qza \
     --type 'FeatureData[Sequence]'
      echo "DONE"
else echo "ALREADY DONE"
fi
date
#########################################               #########################################
print_centered "8 - qiime tools import taxonomy otus"

<< ////
Usage:
qiime tools import --type FeatureData[Taxonomy] --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_/qime2try/99_otu_taxonomy.txt --input-format HeaderlessTSVTaxonomyFormat --output-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_/qime2try/99_otu_taxonomy.qza
Utility: 
  tools               Tools for working with QIIME 2 files.
  import            Import data into a new QIIME 2 Artifact.
  --type TEXT             The semantic type of the artifact that will be
                          created upon importing.
////
if [ ! -e ${OUT}/99_otu_taxonomy.qza ]; then
exe qiime tools import \
     --type FeatureData[Taxonomy] \
     --input-path ${OTU_database}/99_otu_taxonomy.txt \
     --input-format HeaderlessTSVTaxonomyFormat \
     --output-path ${OUT}/99_otu_taxonomy.qza
        echo "DONE"
else echo "ALREADY DONE"
fi
date
#########################################                #########################################
print_centered "9 - qiime feature-classifier extract-reads "
<< ////
Usage:
qiime feature-classifier extract-reads --i-sequences 99_otus.qza --p-f-primer CCTACGGGNGGCWGCAG --p-r-primer GACTACHVGGGTATCTAATCC --o-reads 99_otus-ref.seqs.qza

Utility: 
  Extract simulated amplicon reads from a reference database.
Inputs:
  --i-sequences ARTIFACT FeatureData[Sequence]
Outputs:
  --o-reads ARTIFACT FeatureData[Sequence]
To Do:
////

if [ ! -e ${OUT}/99_otus-ref.seqs.qza ]; then
     begin=$SECONDS
    exe qiime feature-classifier extract-reads \
     --i-sequences ${OUT}/99_otus.qza \
     --p-f-primer CCTACGGGNGGCWGCAG \
     --p-r-primer GACTACHVGGGTATCTAATCC \
     --o-reads ${OUT}/99_otus-ref.seqs.qza

     echo "DONE"
else echo "ALREADY DONE"
fi
date

#########################################               #########################################
print_centered "10 - qiime feature-classifier fit-classifier-naive-bayes"
<< ////
Usage:
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads 99_otus-ref.seqs.qza --i-reference-taxonomy 99_otu_taxonomy.qza --o-classifier classifier.trained.qza

Utility: 
  Create a scikit-learn naive_bayes classifier for reads
Input:
  --i-reference-reads ARTIFACT FeatureData[Sequence]
  --i-reference-taxonomy ARTIFACT FeatureData[Taxonomy]
To Do:
////

if [ ! -e ${OUT}/classifier.trained.qza ]; then
     begin=$SECONDS
     exe qiime feature-classifier fit-classifier-naive-bayes \
     --i-reference-reads ${OUT}/99_otus-ref.seqs.qza \
     --i-reference-taxonomy ${OUT}/99_otu_taxonomy.qza \
     --o-classifier ${OUT}/classifier.trained.qza
     echo "DONE"
else echo "ALREADY DONE"
fi
date
#########################################               #########################################
print_centered "11 - qiime feature-classifier classify-sklearn"
<< ////
Usage:
qiime feature-classifier classify-sklearn --i-classifier classifier.trained.qza --i-reads newCpositive-rep-seqs-dada2.qza --o-classification newtaxonomy.sklearn.qza 

Utility: 
  Classify reads by taxon using a fitted classifier.

To Do:
Notice: I change newCpositive-rep-seqs-dada2.qza to Run2_OB_F-rep-seqs-dada2.qza for tests
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/newtaxonomy.sklearn_${el}.qza  ]; then
     begin=$SECONDS
     exe qiime feature-classifier classify-sklearn \
     --i-classifier ${OUT}/classifier.trained.qza \
     --i-reads ${MERGED}/${el}/data-dada2_${el}_filtered.qza \
     --o-classification ${MERGED}/${el}/newtaxonomy.sklearn_${el}.qza 
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # #########################################               #########################################
print_centered "12  - qiime taxa collapse "
<< ////
Usage:
qiime taxa collapse --i-table table-dada2.qza --i-taxonomy taxonomy.sklearn.qza --p-level 2 --o-collapsed-table FOB_freq-2.qza

Utility: 
  Collapse groups of features that have the same taxonomic assignment
  through the specified level. The frequencies of all features will be
  summed when they are collapsed.

To Do:
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/freq_${el}.qza ]; then
     begin=$SECONDS
     exe qiime taxa collapse \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --i-taxonomy ${MERGED}/${el}/newtaxonomy.sklearn_${el}.qza   \
     --p-level 2 \
     --o-collapsed-table ${MERGED}/${el}/freq_${el}.qza

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # # ########################################               #########################################
print_centered " 13 - qiime feature-table relative-frequency"
<< ////
Usage:
qiime feature-table relative-frequency --i-table FOB_freq-2.qza --o-relative-frequency-table relative-freq-2.qza
Utility: 
Convert frequencies to relative frequencies by dividing each frequency in
  a sample by the sum of frequencies in that sample.

Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]
                       The feature table to be converted into relative
                       frequencies.
Outputs:
  --o-relative-frequency-table ARTIFACT FeatureTable[RelativeFrequency]
                       The resulting relative frequency feature table.
To Do:
////

for el in F FOB S SOB; do
if [ ! -e  ${MERGED}/${el}/relative_freq_${el}.qza ]; then
     begin=$SECONDS
     exe qiime feature-table relative-frequency \
     --i-table ${MERGED}/${el}/freq_${el}.qza \
     --o-relative-frequency-table ${MERGED}/${el}/relative_freq_${el}.qza \
     --verbose > ${LOG}/qiime-feature-table-relative-frequency_${el}.log 2>&1
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################            #########################################
print_centered "14 - qiime alignment mafft "
<< ////
Usage:
qiime alignment mafft --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza

Utility: 
  Perform de novo multiple sequence alignment using MAFFT.
Input:
  --i-sequences ARTIFACT FeatureData[Sequence]
                          The sequences to be aligned.
Output:
  --o-alignment ARTIFACT FeatureData[AlignedSequence]

To Do:
////
for el in F FOB S SOB; do
if [ ! -e  ${MERGED}/${el}/aligned-rep-seqs_${el}.qza ]; then
     begin=$SECONDS
     qiime alignment mafft \
     --i-sequences ${MERGED}/${el}/data-dada2_${el}_filtered.qza \
     --o-alignment  ${MERGED}/${el}/aligned-rep-seqs_${el}.qza

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "15 - qiime alignment mask"
<< ////
Usage:
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

Utility: 
 Mask (i.e., filter) unconserved and highly gapped columns from an
  alignment. Default min_conservation was chosen to reproduce the mask
  presented in Lane (1991)

To Do:
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/masked-aligned-rep-seqs_${el}.qza ]; then
     begin=$SECONDS
     exe qiime alignment mask \
     --i-alignment ${MERGED}/${el}/aligned-rep-seqs_${el}.qza \
     --o-masked-alignment ${MERGED}/${el}/masked-aligned-rep-seqs_${el}.qza

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "16 - qiime phylogeny fasttree "
<< ////
Usage:
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

Utility: 
  Construct a phylogenetic tree with FastTree.
Inputs:
  --i-alignment ARTIFACT FeatureData[AlignedSequence]
                          Aligned sequences to be used for phylogenetic
                          reconstruction.   

To Do:
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/unrooted-tree_${el}.qza ]; then
     begin=$SECONDS
     exe qiime phylogeny fasttree \
     --i-alignment ${MERGED}/${el}/masked-aligned-rep-seqs_${el}.qza \
     --o-tree ${MERGED}/${el}/unrooted-tree_${el}.qza

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "17 - qiime phylogeny midpoint-root"
<< ////
Usage:
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza 

Utility: 
  Midpoint root an unrooted phylogenetic tree.
Inputs:
  --i-tree ARTIFACT Phylogeny[Unrooted]
                       The phylogenetic tree to be rooted.          
Outputs:
  --o-rooted-tree ARTIFACT
    Phylogeny[Rooted]  The rooted phylogenetic tree.   

To Do:
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/rooted-tree_${el}.qza  ]; then
     begin=$SECONDS
     exe qiime phylogeny midpoint-root \
     --i-tree ${MERGED}/${el}/unrooted-tree_${el}.qza \
     --o-rooted-tree ${MERGED}/${el}/rooted-tree_${el}.qza
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # #########################################        BUG      #########################################
print_centered "18 - qiime diversity alpha-rarefaction "
<< ////
Usage:
#qiime diversity alpha-rarefaction --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-max-depth 27828 --p-steps 75 --p-iterations 55 --m-metadata-file Metadata_FOB_BS.tsv --p-metrics chao1 --p-metrics simpson_e --p-metrics simpson --p-metrics shannon --p-metrics observed_otus --p-metrics faith_pd --o-visualization rarefaction-curve.qzv
Utility: 
  Generate interactive alpha rarefaction curves by computing rarefactions
  between min_depth and max_depth. The number of intermediate depths to
  compute is controlled by the steps parameter, with n iterations being
  computed at each rarefaction depth. If sample metadata is provided,
  samples may be grouped based on distinct values within a metadata column.
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]

Outputs:
Outputs:
  --o-visualization VISUALIZATION

To Do:
////
for el in F FOB S SOB; do
if [ ! -e ${MERGED}/${el}/rarefaction-curve_${el}.qzv ]; then
     begin=$SECONDS
     qiime diversity alpha-rarefaction \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --i-phylogeny ${MERGED}/${el}/rooted-tree_${el}.qza \
     --p-iterations 55 \
     --m-metadata-file ${MetaNames}${el}".csv" \
     --p-max-depth 27828 \
     --p-steps 75 \
     --p-metrics chao1 \
     --p-metrics simpson_e \
     --p-metrics simpson \
     --p-metrics shannon \
     --p-metrics observed_otus \
     --p-metrics faith_pd \
     --o-visualization ${MERGED}/${el}/rarefaction-curve_${el}.qzv

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################       seems to work        #########################################
print_centered "19 - qiime diversity beta"
<< ////
Usage:
qiime diversity beta --i-table table-dada2.qza --p-metric braycurtis --o-distance-matrix dada2.braycurtis.notNorm.diversity.qza
Utility: 
  Computes a user-specified beta diversity metric for all pairs of samples
  in a feature table.
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency | RelativeFrequency |
    PresenceAbsence]   The feature table containing the samples over which
                       beta diversity should be computed. 
Outputs:
  --o-distance-matrix ARTIFACT
    DistanceMatrix     The resulting distance matrix.

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/dada2_braycurtis_notNorm_diversity_${el}.qza ]; then
     begin=$SECONDS
     exe qiime diversity beta \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --p-metric braycurtis \
     --o-distance-matrix ${MERGED}/${el}/dada2_braycurtis_notNorm_diversity_${el}.qza

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "20 - qiime diversity pcoa"
<< ////
Usage:
qiime diversity pcoa --i-distance-matrix dada2.braycurtis.notNorm.diversity.qza --o-pcoa dada2.braycurtis.notNorm.diversity.pcoa.qza

Utility: 
  Apply principal coordinate analysis.

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/dada2_braycurtis_notNorm_diversity_pcoa_${el}.qza ]; then
     begin=$SECONDS
     exe qiime diversity pcoa \
     --i-distance-matrix ${MERGED}/${el}/dada2_braycurtis_notNorm_diversity_${el}.qza \
     --o-pcoa ${MERGED}/${el}/dada2_braycurtis_notNorm_diversity_pcoa_${el}.qza
     
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # #########################################       BUG        #########################################
print_centered "21 - qiime diversity core-metrics-phylogenetic"
<< ////
Usage:
qiime diversity core-metrics-phylogenetic --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-sampling-depth 5801 --output-dir dada2-diversity-5801 --m-metadata-file Metadata_FOB_BS.tsv

Utility: 
 Applies a collection of diversity metrics (both phylogenetic and non-
  phylogenetic) to a feature table.
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]
                          The feature table containing the samples over which
                          diversity metrics should be computed.    
  --i-phylogeny ARTIFACT  Phylogenetic tree containing tip identifiers that
    Phylogeny[Rooted]     correspond to the feature identifiers in the table.
                          This tree can contain tip ids that are not present
                          in the table, but all feature ids in the table must
                          be present in this tree.                  

To Do:
////
for el in F FOB S SOB; do

if [ ! -d ${MERGED}/${el}/diversity ]; then
     begin=$SECONDS
     exe qiime diversity core-metrics-phylogenetic \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --m-metadata-file ${MetaNames}${el}".csv" \
     --i-phylogeny ${MERGED}/${el}/rooted-tree_${el}.qza \
     --p-sampling-depth 5801 \
     --output-dir ${MERGED}/${el}/diversity
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "22 - qiime gneiss correlation-clustering"
<< ////
Usage:
qiime gneiss correlation-clustering --i-table table-dada2.qza --p-pseudocount 1 --o-clustering hierarchy1.qza         

Utility: 
  Build a bifurcating tree that represents a hierarchical clustering of
  features.  The hiearchical clustering uses Ward hierarchical clustering
  based on the degree of proportionality between features.
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]
                          The feature table containing the samples in which
                          the columns will be clustered.  

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/hierarchy_${el}.qza  ]; then
     begin=$SECONDS
     exe qiime gneiss correlation-clustering \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --p-pseudocount 1 \
     --o-clustering ${MERGED}/${el}/hierarchy_${el}.qza         

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # #########################################               #########################################
print_centered "23 - qiime gneiss gradient-clustering"
<< ////
Usage:
qiime gneiss gradient-clustering --i-table table-dada2.qza --m-gradient-file Metadata_FOB_BS.tsv --m-gradient-column TimePoint --o-clustering gradient-hierarchy.qza 

Utility: 

To Do:
////
for el in F FOB SOB; do
if [ ! -e ${MERGED}/${el}/gradient-hierarchy_${el}.qza ]; then
     begin=$SECONDS
     exe qiime gneiss gradient-clustering \
      --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
      --m-gradient-file ${MetaNames}${el}".csv" \
      --m-gradient-column BMI \
      --o-clustering ${MERGED}/${el}/gradient-hierarchy_${el}.qza
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "24 - qiime gneiss ilr-hierarchical"
<< ////
Usage:
qiime gneiss ilr-hierarchical --i-table table-dada2.qza --i-tree hierarchy1.qza --o-balances balances.qza

Utility: 
  Calculate balances given a hierarchy.
Inputs:
  --i-table ARTIFACT FeatureTable[Frequency | Composition]
                          The feature table containing the samples in which
                          the ilr transform will be performed.      
  --i-tree ARTIFACT       A hierarchy of feature identifiers that defines the
    Hierarchy             partitions of features.  Each tip in the
                          hierarchycorresponds to the feature identifiers in
                          the table. This tree can contain tip ids that are
                          not present in the table, but all feature ids in the
                          table must be present in this tree.  This assumes
                          that all of the internal nodes in the tree have
                          labels.                                   

Outputs:
  --o-balances ARTIFACT FeatureTable[Balance]
                          The resulting balances from the ilr transform.
                                                                    

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/balances_${el}.qza ]; then
     begin=$SECONDS
     exe qiime gneiss ilr-hierarchical \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza  \
     --i-tree ${MERGED}/${el}/hierarchy_${el}.qza  \
     --o-balances ${MERGED}/${el}/balances_${el}.qza
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "qiime tools export"
<< ////
Usage:
qiime tools export --input-path table-dada2.qza --output-path exported-table-dada2.tsv
Utility: 

To Do:
////
if [ ! -e ${MERGED}/export ]; then
     begin=$SECONDS
    exe qiime tools export \
      --input-path ${MERGED}/table-dada2_merged.qza \
      --output-path ${MERGED}/export
     echo "DONE"
else echo "ALREADY DONE"
fi

# #########################################               #########################################
print_centered "biom convert"
<< ////
Usage:
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

Utility: 

To Do:
////
if [ ! -e ${MERGED}/export/feature-table.tsv ]; then
     begin=$SECONDS
    exe biom convert \
    -i ${MERGED}/export/feature-table.biom -o ${MERGED}/export/feature-table.tsv --to-tsv
     echo "DONE"
else echo "ALREADY DONE"
fi

# #########################################               #########################################
print_centered "summerize"
<< ////
Usage:
qiime feature-table summarize --i-table table-dada2.qza --o-visualization table-dada2.qzv --m-sample-metadata-file metadata.tsv

Utility: 

To Do:
////
if [ ! -e ${MERGED}/table-dada2_merged.qzv  ]; then
     begin=$SECONDS
    exe qiime feature-table summarize \
      --i-table ${MERGED}/table-dada2_merged.qza \
      --o-visualization ${MERGED}/table-dada2_merged.qzv \
      --m-sample-metadata-file ${MetaNames}"_all.csv"
    
    echo "DONE"
else echo "ALREADY DONE"
fi


# #########################################               #########################################
print_centered "25 - qiime picrust2 full-pipeline"
<< ////
beware: not compatible with qiime2-2020.8 but with 2019.10 due to some dependancies
wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml --no-check-certificate
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.10-py36-linux-conda.yml

DO : conda activate qiime2-2019.10 
and install q2-picrust2 with conda
Instruction :  conda install -c gavinmdouglas q2-picrust2 
Usage:

qiime picrust2 full-pipeline --i-table table-dada2.qza --i-seq rep-seqs-dada2.qza --output-dir q2-picrust2_output --p-threads 1 --p-hsp-method pic --p-max-nsti 2 --verbose

Utility: 
  QIIME 2 plugin for default 16S PICRUSt2 pipeline

To Do: need to install qiime2 https://forum.qiime2.org/t/qiime2-modulenotfounderror-with-picrust2-plugin-installation/8985/2 
////

for el in F FOB S SOB; do

if [ ! -d ${MERGED}/${el}/picrust2 ]; then
     begin=$SECONDS
     exe qiime picrust2 full-pipeline \
     --i-table ${MERGED}/${el}/run-rep-dada2_${el}.qza \
     --i-seq ${MERGED}/${el}/data-dada2_${el}_filtered.qza \
     --output-dir ${MERGED}/${el}/picrust2 \
     --p-threads ${threads} \
     --p-hsp-method pic \
     --p-max-nsti 2 --verbose > ${LOG}/picrust2_${el}.log 2>&1 

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "26 -  qiime diversity core-metrics "
<< ////
Usage:
qiime diversity core-metrics --i-table q2-picrust2_output/pathway_abundance.qza --p-sampling-depth 	548854 --m-metadata-file Metadata_FOB_BS.tsv --output-dir pathabun_core_metrics_out --p-n-jobs 1

Inputs:
  --i-table ARTIFACT FeatureTable[Frequency]
                       The feature table containing the samples over which
                       diversity metrics should be computed.

Utility: 
      Applies a collection of diversity metrics (non-phylogenetic) to a feature
        table.

exout: pathabun_core_metrics_out
To Do:
Supposition : ${OUT}/pathway_abundance.qza must have been created by the picrust2 command. 
////
for el in F FOB S SOB; do

if [ ! -d ${MERGED}/${el}/core_metrics ]; then
    begin=$SECONDS
    exe qiime diversity core-metrics \
      --i-table ${MERGED}/${el}/picrust2/pathway_abundance.qza \
      --p-sampling-depth 	548854 \
      --m-metadata-file ${MetaNames}${el}".csv" \
      --output-dir ${MERGED}/${el}/core_metrics 
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# # #########################################               #########################################
print_centered "qiime tools export"
print_centered "creation of feature-table.biom files"
<< ////
Usage:
qiime tools export --input-path q2-picrust2_output/pathway_abundance.qza --output-path pathabun_exported 

Utility: 

To Do:
////
for el in F FOB S SOB; do

if [ ! -d ${MERGED}/${el}/export_picrust2 ]; then
    begin=$SECONDS
    exe qiime tools export \
      --input-path ${MERGED}/${el}/picrust2/pathway_abundance.qza \
      --output-path ${MERGED}/${el}/export_picrust2

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "result"
<< ////
Usage:
biom convert -i pathabun_exported/feature-table.biom -o pathabun_exported/feature-table.biom.tsv --to-tsv

Utility: 

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/export_picrust2/feature-table_biom_${el}.tsv ]; then
     begin=$SECONDS
      exe biom convert -i ${MERGED}/${el}/export_picrust2/feature-table.biom \
      -o ${MERGED}/${el}/export_picrust2/feature-table_biom_${el}.tsv --to-tsv
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# ########################################               #########################################
print_centered "27 - qiime feature-table relative-frequency "
<< ////
Usage:
qiime feature-table relative-frequency --i-table pathway_abundance.qza --o-relative-frequency-table pathway_rel-freq.qza   

Utility: Convert frequencies to relative frequencies by dividing each frequency in
  a sample by the sum of frequencies in that sample.

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/picrust2/pathway_rel-freq_${el}.qza ]; then
     begin=$SECONDS
     exe qiime feature-table relative-frequency \
      --i-table ${MERGED}/${el}/picrust2/pathway_abundance.qza \
      --o-relative-frequency-table ${MERGED}/${el}/picrust2/pathway_rel-freq_${el}.qza   
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# #########################################               #########################################
print_centered "qiime tools export pathway_rel-freq.qza "
print_centered "creation of feature-table.biom files"
<< ////
Usage:
qiime tools export --input-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq.qza --output-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq

Utility: 

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/path_rel-freq_picrust2    ]; then
     begin=$SECONDS
    exe qiime tools export \
    --input-path ${MERGED}/${el}/picrust2/pathway_rel-freq_${el}.qza    \
    --output-path ${MERGED}/${el}/path_rel-freq_picrust2  

     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date
# ########################################               #########################################
print_centered "biom convert pathway ref freq to tsv"
<< ////
Usage:
biom convert -i feature-table.biom -o pathway_rel-freq_tsv.txt --to-tsv

Utility: 

To Do:
////
for el in F FOB S SOB; do

if [ ! -e ${MERGED}/${el}/path_rel-freq_picrust2/pathway_rel-freq_${el}_tsv.txt  ]; then
     begin=$SECONDS
    exe biom convert \
    -i ${MERGED}/${el}/path_rel-freq_picrust2/feature-table.biom \
    -o ${MERGED}/${el}/path_rel-freq_picrust2/pathway_rel-freq_${el}_tsv.txt --to-tsv
     echo -e ${el} "\tDONE"
else echo -e ${el} "\tALREADY DONE"
fi 
done
date







print_centered "X - end of the analysis"

date 
