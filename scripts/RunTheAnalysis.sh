#!/usr/bin/bash
##########################################mise en page#########################################
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
OTU_database=${DIR}"/database/greengenes"
threads=4


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
        --p-cores ${threads} \
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

if [ ! -e ${OUT}/oktable-dada2.qza ]; then
    begin=$SECONDS
qiime dada2 denoise-paired \
     --i-demultiplexed-seqs ${OUT}/paired-end-demux.trim.qza \
     --p-trim-left-f 20 \
     --p-trim-left-r 0 \
     --p-trunc-len-f 290 \
     --p-trunc-len-r 260 \
     --p-no-hashed-feature-ids \
     --p-n-threads ${threads} \
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
print_centered "6 - qiime tools import otus fasta"

<< ////
Usage:
qiime tools import --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otus.fasta --output-path 99_otus.qza --type 'FeatureData[Sequence]'
Utility: 
  tools               Tools for working with QIIME 2 files.
  import            Import data into a new QIIME 2 Artifact.
  --input-path PATH       Path to file or directory that should be imported.
  --output-path ARTIFACT  Path where output artifact should be written.
  --type TEXT             The semantic type of the artifact that will be
                          created upon importing.
////
if [ ! -e ${OUT}/99_otus.qza ]; then
qiime tools import \
     --input-path ${OTU_database}/99_otus.fasta \
     --output-path ${OUT}/99_otus.qza \
     --type 'FeatureData[Sequence]'
        echo "DONE"
else echo "ALREADY DONE"
fi

#########################################       7       #########################################
print_centered "7 - qiime tools import taxonomy otus"

<< ////
Usage:
qiime tools import --type FeatureData[Taxonomy] --input-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otu_taxonomy.txt --input-format HeaderlessTSVTaxonomyFormat --output-path /data/icb/16sRNA/work/CRC_obese_microbiome_run1_ok/qime2try/99_otu_taxonomy.qza
Utility: 
  tools               Tools for working with QIIME 2 files.
  import            Import data into a new QIIME 2 Artifact.
  --type TEXT             The semantic type of the artifact that will be
                          created upon importing.
////
if [ ! -e ${OUT}/99_otu_taxonomy.qza ]; then
qiime tools import \
     --type FeatureData[Taxonomy] \
     --input-path ${OTU_database}/99_otu_taxonomy.txt \
     --input-format HeaderlessTSVTaxonomyFormat \
     --output-path ${OUT}/99_otu_taxonomy.qza
        echo "DONE"
else echo "ALREADY DONE"
fi

#########################################                #########################################
print_centered " - qiime feature-classifier extract-reads "
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
     qiime feature-classifier extract-reads \
     --i-sequences ${OUT}/99_otus.qza \
     --p-f-primer CCTACGGGNGGCWGCAG \
     --p-r-primer GACTACHVGGGTATCTAATCC \
     --o-reads ${OUT}/99_otus-ref.seqs.qza

     echo "DONE"
else echo "ALREADY DONE"
fi


#########################################               #########################################
print_centered " - qiime feature-classifier fit-classifier-naive-bayes"
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
     qiime feature-classifier fit-classifier-naive-bayes \
     --i-reference-reads ${OUT}/99_otus-ref.seqs.qza \
     --i-reference-taxonomy ${OUT}/99_otu_taxonomy.qza \
     --o-classifier ${OUT}/classifier.trained.qza
     echo "DONE"
else echo "ALREADY DONE"
fi

#########################################               #########################################
print_centered " - qiime feature-classifier classify-sklearn"
<< ////
Usage:
qiime feature-classifier classify-sklearn --i-classifier classifier.trained.qza --i-reads newCpositive-rep-seqs-dada2.qza --o-classification newtaxonomy.sklearn.qza 

Utility: 
  Classify reads by taxon using a fitted classifier.

To Do:
Notice: I change newCpositive-rep-seqs-dada2.qza to Run2_OB_F-rep-seqs-dada2.qza for tests
////
if [ ! -e ${OUT}/newtaxonomy.sklearn.qza  ]; then
     begin=$SECONDS
     qiime feature-classifier classify-sklearn \
     --i-classifier ${OUT}/classifier.trained.qza \
     --i-reads ${OUT}/Run2_OB_F-table-dada2.qza \
     --o-classification ${OUT}/newtaxonomy.sklearn.qza 
     echo "DONE"
else echo "ALREADY DONE"
fi


# #########################################    BUG           #########################################
# # print_centered " bug - qiime taxa collapse "
# # << ////
# # Usage:
# # qiime taxa collapse --i-table table-dada2.qza --i-taxonomy taxonomy.sklearn.qza --p-level 2 --o-collapsed-table FOB_freq-2.qza

# # Utility: 
# #   Collapse groups of features that have the same taxonomic assignment
# #   through the specified level. The frequencies of all features will be
# #   summed when they are collapsed.

# # To Do:
# # ////
# # if [ ! -e ${OUT}/FOB_freq-2.qza ]; then
# #      begin=$SECONDS
# #      qiime taxa collapse \
# #      --i-table ${OUT}/oktable-dada2.qza \
# #      --i-taxonomy ${OUT}/newtaxonomy.sklearn.qza  \
# #      --p-level 2 \
# #      --o-collapsed-table ${OUT}/FOB_freq-2.qza

# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi
# # ########################################        BUG       #########################################
# #  print_centered " bug - qiime feature-table relative-frequency"
# # << ////
# # Usage:
# # qiime feature-table relative-frequency --i-table FOB_freq-2.qza --o-relative-frequency-table relative-freq-2.qza
# # Utility: 
# # Convert frequencies to relative frequencies by dividing each frequency in
# #   a sample by the sum of frequencies in that sample.

# # Inputs:
# #   --i-table ARTIFACT FeatureTable[Frequency]
# #                        The feature table to be converted into relative
# #                        frequencies.
# # Outputs:
# #   --o-relative-frequency-table ARTIFACT FeatureTable[RelativeFrequency]
# #                        The resulting relative frequency feature table.
# # To Do:
# # ////

# # if [ ! -e ${OUT}/relative-freq-2.qza ]; then
# #      begin=$SECONDS
# #      qiime feature-table relative-frequency \
# #      --i-table ${OUT}/FOB_freq-2.qza \
# #      --o-relative-frequency-table ${OUT}/relative-freq-2.qza \
# #      --verbose > ${LOG}/qiime-feature-table-relative-frequency.log
# #      echo "DONE"
# # else echo "ALREADY DONE"
# # fi


#########################################       TODO        #########################################
print_centered " to do - biom convert "  
<< ////
Usage:
biom convert -i feature-table.biom -o FOB_relative_freq5_tsv.txt --to-tsv

Utility: 
  Convert to/from the BIOM table format.


To Do: do it :-( 
////

#########################################      should work       #########################################
print_centered " - qiime alignment mafft "
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
if [ ! -e ${OUT}/aligned-rep-seqs.qza ]; then
     begin=$SECONDS
     qiime alignment mafft \
     --i-sequences ${OUT}/okrep-seqs-dada2.qza  \
     --o-alignment ${OUT}/aligned-rep-seqs.qza

     echo "DONE"
else echo "ALREADY DONE"
fi
#########################################               #########################################
print_centered " - qiime alignment mask"
<< ////
Usage:
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

Utility: 
 Mask (i.e., filter) unconserved and highly gapped columns from an
  alignment. Default min_conservation was chosen to reproduce the mask
  presented in Lane (1991)

To Do:
////
if [ ! -e ${OUT}/masked-aligned-rep-seqs.qza ]; then
     begin=$SECONDS
     qiime alignment mask \
     --i-alignment ${OUT}/aligned-rep-seqs.qza \
     --o-masked-alignment ${OUT}/masked-aligned-rep-seqs.qza

     echo "DONE"
else echo "ALREADY DONE"
fi
#########################################               #########################################
print_centered " - qiime phylogeny fasttree "
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
if [ ! -e ${OUT}/unrooted-tree.qza ]; then
     begin=$SECONDS
     qiime phylogeny fasttree \
     --i-alignment ${OUT}/masked-aligned-rep-seqs.qza \
     --o-tree ${OUT}/unrooted-tree.qza

     echo "DONE"
else echo "ALREADY DONE"
fi
#########################################               #########################################
print_centered " - qiime phylogeny midpoint-root"
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
if [ ! -e ${OUT}/rooted-tree.qza  ]; then
     begin=$SECONDS
     qiime phylogeny midpoint-root \
     --i-tree ${OUT}/unrooted-tree.qza \
     --o-rooted-tree ${OUT}/rooted-tree.qza 

     echo "DONE"
else echo "ALREADY DONE"
fi
# #########################################        BUG      #########################################
# print_centered " - qiime diversity alpha-rarefaction "

# << ////
# Usage:
# qiime diversity alpha-rarefaction --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-max-depth 27828 --p-steps 75 --p-iterations 55 --m-metadata-file Metadata_FOB_BS.tsv --p-metrics chao1 --p-metrics simpson_e --p-metrics simpson --p-metrics shannon --p-metrics observed_otus --p-metrics faith_pd --o-visualization rarefaction-curve.qzv

# Utility: 
#   compute is controlled by the `steps` parameter, with n `iterations` being
#   computed at each rarefaction depth. If sample metadata is provided,
#   samples may be grouped based on distinct values within a metadata column.
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency]
#                           Feature table to compute rarefaction curves from.
#   --i-phylogeny ARTIFACT  Optional phylogeny for phylogenetic metrics.
#     Phylogeny[Rooted]  


# To Do:
#      --p-max-depth 27828 \
#      --p-steps 75 \

# ////

# if [ ! -e ${OUT}/rarefaction-curve.qzv ]; then
#      begin=$SECONDS
#      qiime diversity alpha-rarefaction \
#      --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
#      --i-phylogeny ${OUT}/rooted-tree.qza \
#      --p-iterations 55 \
#      --m-metadata-file ${MetaNames} \
#      --p-max-depth 27828 \
#      --p-steps 75 \
#      --p-metrics chao1 \
#      --p-metrics simpson_e \
#      --p-metrics simpson \
#      --p-metrics shannon \
#      --p-metrics observed_otus \
#      --p-metrics faith_pd \
#      --o-visualization ${OUT}/rarefaction-curve.qzv

#      echo "DONE"
# else echo "ALREADY DONE"
# fi

#########################################       seems to work        #########################################
print_centered " - qiime diversity beta"
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
if [ ! -e ${OUT}/dada2.braycurtis.notNorm.diversity.qza ]; then
     begin=$SECONDS
     qiime diversity beta \
     --i-table ${OUT}/Run2_OB_F-table-dada2.qza \
     --p-metric braycurtis \
     --o-distance-matrix ${OUT}/dada2.braycurtis.notNorm.diversity.qza

     echo "DONE"
else echo "ALREADY DONE"
fi

#########################################               #########################################
print_centered " - qiime diversity pcoa"
<< ////
Usage:
qiime diversity pcoa --i-distance-matrix dada2.braycurtis.notNorm.diversity.qza --o-pcoa dada2.braycurtis.notNorm.diversity.pcoa.qza

Utility: 
  Apply principal coordinate analysis.

To Do:
////
if [ ! -e ${OUT}/dada2.braycurtis.notNorm.diversity.pcoa.qza ]; then
     begin=$SECONDS
     qiime diversity pcoa \
     --i-distance-matrix ${OUT}/dada2.braycurtis.notNorm.diversity.qza \
     --o-pcoa ${OUT}/dada2.braycurtis.notNorm.diversity.pcoa.qza
     
     echo "DONE"
else echo "ALREADY DONE"
fi


# #########################################       BUG        #########################################
# print_centered " - qiime diversity core-metrics-phylogenetic"
# << ////
# Usage:
# qiime diversity core-metrics-phylogenetic --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-sampling-depth 5801 --output-dir dada2-diversity-5801 --m-metadata-file Metadata_FOB_BS.tsv

# Utility: 
#  Applies a collection of diversity metrics (both phylogenetic and non-
#   phylogenetic) to a feature table.
# Inputs:
#   --i-table ARTIFACT FeatureTable[Frequency]
#                           The feature table containing the samples over which
#                           diversity metrics should be computed.    
#   --i-phylogeny ARTIFACT  Phylogenetic tree containing tip identifiers that
#     Phylogeny[Rooted]     correspond to the feature identifiers in the table.
#                           This tree can contain tip ids that are not present
#                           in the table, but all feature ids in the table must
#                           be present in this tree.                  

# To Do:
# ////
# if [ ! -d diversity ]; then
#      mkdir ${DIR}/diversity
#      begin=$SECONDS
#      qiime diversity core-metrics-phylogenetic \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --i-phylogeny r ${OUT}/rooted-tree.qza \
#      --p-sampling-depth 5801 \
#      --output-dir ${DIR}/diversity \
#      --m-metadata-file ${MetaNames}
#      echo "DONE"
# else echo "ALREADY DONE"
# fi

#########################################               #########################################
print_centered " - qiime gneiss correlation-clustering"
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
if [ ! -e ${OUT}/hierarchy1.qza  ]; then
     begin=$SECONDS
     qiime gneiss correlation-clustering \
     --i-table ${OUT}/oktable-dada2.qza \
     --p-pseudocount 1 \
     --o-clustering ${OUT}/hierarchy1.qza         

     echo "DONE"
else echo "ALREADY DONE"
fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime gneiss gradient-clustering --i-table table-dada2.qza --m-gradient-file Metadata_FOB_BS.tsv --m-gradient-column TimePoint --o-clustering gradient-hierarchy.qza 

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS
#      qiime gneiss gradient-clustering --i-table table-dada2.qza --m-gradient-file Metadata_FOB_BS.tsv --m-gradient-column TimePoint --o-clustering gradient-hierarchy.qza 

#      echo "DONE"
# else echo "ALREADY DONE"
# fi

#########################################               #########################################
print_centered " - qiime gneiss ilr-hierarchical"
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
if [ ! -e ${OUT}/balances.qza ]; then
     begin=$SECONDS
     qiime gneiss ilr-hierarchical \
     --i-table ${OUT}/oktable-dada2.qza  \
     --i-tree ${OUT}/hierarchy1.qza \
     --o-balances ${OUT}/balances.qza

     echo "DONE"
else echo "ALREADY DONE"
fi


# #########################################    feature table biom that I don't have            #########################################
# print_centered " - "
# << ////
# Usage:
# qiime tools import --input-path feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table2.qza    

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
#########################################               #########################################
print_centered " - qiime picrust2 full-pipeline"
# << ////
# Usage:
# qiime picrust2 full-pipeline --i-table table-dada2.qza --i-seq rep-seqs-dada2.qza --output-dir q2-picrust2_output --p-threads 1 --p-hsp-method pic --p-max-nsti 2 --verbose

# Utility: 

# To Do: need to install qiime2 https://forum.qiime2.org/t/qiime2-modulenotfounderror-with-picrust2-plugin-installation/8985/2 
# ////
# if [ ! -d ${DIR}/q2-picrust2_output ]; then
#      mkdir ${DIR}/q2-picrust2_output
#      begin=$SECONDS
#      qiime picrust2 full-pipeline \
#      --i-table ${OUT}/oktable-dada2.qza \
#      --i-seq ${OUT}/okrep-seqs-dada2.qza \
#      --output-dir q2-picrust2_output \
#      --p-threads ${threads} \
#      --p-hsp-method pic \
#      --p-max-nsti 2 --verbose

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime feature-table summarize --i-table q2-picrust2_output/pathway_abundance.qza --o-visualization q2-picrust2_output/pathway_abundance.qzv

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime diversity core-metrics --i-table q2-picrust2_output/pathway_abundance.qza --p-sampling-depth 	548854 --m-metadata-file Metadata_FOB_BS.tsv --output-dir pathabun_core_metrics_out --p-n-jobs 1

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime tools export --input-path q2-picrust2_output/pathway_abundance.qza --output-path pathabun_exported 

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi

# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# biom convert -i pathabun_exported/feature-table.biom -o pathabun_exported/feature-table.biom.tsv --to-tsv

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime feature-table relative-frequency --i-table pathway_abundance.qza --o-relative-frequency-table pathway_rel-freq.qza   

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# qiime tools export --input-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq.qza --output-path /data/icb/16sRNA_DAVID/work/ob_Bariatric_Surgery/feces/q2-picrust2_output/pathway_rel-freq

# Utility: 

# To Do:
# ////

# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi
# #########################################               #########################################
# print_centered " - "
# << ////
# Usage:
# biom convert -i feature-table.biom -o pathway_rel-freq_tsv.txt --to-tsv

# Utility: 

# To Do:
# ////
# if [ ! -e ${OUT}/... ]; then
#      begin=$SECONDS

#      echo "DONE"
# else echo "ALREADY DONE"
# fi









print_centered "X - end of the analysis"
date 