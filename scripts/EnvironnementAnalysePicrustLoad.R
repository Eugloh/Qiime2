#!/usr/bin/Rscript
############################################################################################################################################
# Set the environment for the script dedicated to picrust analysis 
#
# I inspired myself from the setup_microbiome_analysis.R of Microbial-bioinformatics-introductory-course-Material-2018
# a online course available at : https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html
# Sudarshan A. Shetty, Leo Lahti, Gerben DA. Hermes
# last uptaded 2020-04-11
# github : https://github.com/mibwurrepo/Microbial-bioinformatics-introductory-course-Material-2018
#
############################################################################################################################################

options(repos = "https://pbil.univ-lyon1.fr/CRAN/")

setup_ALDEx2_analysis <- function(){
  
  .packages = c("dplyr" , "devtools","optparse")
  
  .bioc_packages <- c("phyloseq","graphite", "SBGNview", "ALDEx2")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	
  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

  .inst <- .bioc_packages %in% installed.packages()
  if(any(!.inst)) {
    BiocManager::install(.bioc_packages[!.inst])
  }
  
  message("If there was no error then you are ready to launch the analysis")
  
}

setup_ALDEx2_analysis()
