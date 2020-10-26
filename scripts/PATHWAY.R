#!/usr/bin/Rscript
##############################Info########################
### PICRUST output analysis 
### Jully 2020 
###
### work for automatisation: October 2020
##########################################################

#############################Packages#####################
# check for installation and packages 
source("./scripts/EnvironnementAnalysePicrustLoad.R")

library(SBGNview)
library(phyloseq)
library(ALDEx2)
library(dplyr)
library(optparse)

##########################################################

directory <- "/media/lohmanne/E125-D7651/DAVID_rerun/output/tsv" # directory to set here 
setwd(directory)

option_list = list(
  make_option(c("-p", "--pattern"), type="character", default="tsv", help="pattern for file names [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sampleFiles  <- grep(opt$pattern,list.files(directory),value=TRUE) # find names

TABLES <- list()
for (el in sampleFiles){
  el_name <- sapply(strsplit(as.character(el), "\\."), "[[", 1)
  TABLES[[el_name]] <- read.delim(paste(directory,el, sep = "/"), row.names=1)
}

 names(TABLES) ## retrieve the table of interest with TABLES[["<name>"]]

#  [1] "ec_metagenome_F"       "ec_metagenome_FOB"     "ec_metagenome_S"      
#  [4] "ec_metagenome_SOB"     "ko_metagenome_F"       "ko_metagenome_FOB"    
#  [7] "ko_metagenome_S"       "ko_metagenome_SOB"     "pathway_abundance_F"  
# [10] "pathway_abundance_FOB" "pathway_abundance_S"   "pathway_abundance_SOB"
# [13] "pathway_rel-freq_F"    "pathway_rel-freq_FOB"  "pathway_rel-freq_S"   
# [16] "pathway_rel-freq_SOB" 

EC <- cbind(TABLES[["ec_metagenome_F"]],TABLES[["ec_metagenome_FOB"]],TABLES[["ec_metagenome_S"]],TABLES[["ec_metagenome_SOB"]])
summary(EC)
### take the colnames to determine sous groups, from the same cohort or at the same time point. 
Sample_name <-c(sort(unique(sapply(strsplit(as.character(colnames(EC)), "\\."), "[[", 1))))

### the information is sorted into lists, I create table of specifics sous groups. 

# creation of the EC tables, concat in a le list. 


le<- list()
for (time in c("T0","T6","T12")){
  le[[time]] <- list()
  for (statut in c("FOB","SOB","F","S")){
    le[[time]][[statut]] <- EC[,grep(paste0(statut, ".*", time, "$"),colnames(EC))]
  }
}


## visualization and statistics
# ALDEx2
# modelling the data as a log-ratio transformed probability distribution rather than as counts

conds <- list()
ECselex.sub <- list()
x.all <- list()

# differential(relative) abundance analysis of proportional data
ECselex.sub[["06FOB"]] <- round(cbind(le[["T0"]][["FOB"]],le[["T6"]][["FOB"]]))
ECselex.sub[["012FOB"]]<- round(cbind(le[["T0"]][["FOB"]],le[["T12"]][["FOB"]]))
ECselex.sub[["612FOB"]] <- round(cbind(le[["T6"]][["FOB"]],le[["T12"]][["FOB"]]))
# differential(relative) abundance analysis of proportional data
ECselex.sub[["06SOB"]] <- round(cbind(le[["T0"]][["SOB"]],le[["T6"]][["SOB"]]))
ECselex.sub[["012SOB"]]<- round(cbind(le[["T0"]][["SOB"]],le[["T12"]][["SOB"]]))
ECselex.sub[["612SOB"]] <- round(cbind(le[["T6"]][["SOB"]],le[["T12"]][["SOB"]]))
# differential(relative) abundance analysis of proportional data
ECselex.sub[["06S"]] <- round(cbind(le[["T0"]][["S"]],le[["T6"]][["S"]]))
ECselex.sub[["012S"]]<- round(cbind(le[["T0"]][["S"]],le[["T12"]][["S"]]))
ECselex.sub[["612S"]] <- round(cbind(le[["T6"]][["S"]],le[["T12"]][["S"]]))
# differential(relative) abundance analysis of proportional data
ECselex.sub[["06F"]] <- round(cbind(le[["T0"]][["F"]],le[["T6"]][["F"]]))
ECselex.sub[["012F"]]<- round(cbind(le[["T0"]][["F"]],le[["T12"]][["F"]]))
ECselex.sub[["612F"]] <- round(cbind(le[["T6"]][["F"]],le[["T12"]][["F"]]))



conds[["06FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["FOB"]],le[["T6"]][["FOB"]]))), "\\."), "[[", 2)
conds[["012FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["FOB"]],le[["T12"]][["FOB"]]))), "\\."), "[[", 2)
conds[["612FOB"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["FOB"]],le[["T12"]][["FOB"]]))), "\\."), "[[", 2)

conds[["06SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["SOB"]],le[["T6"]][["SOB"]]))), "\\."), "[[", 2)
conds[["012SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["SOB"]],le[["T12"]][["SOB"]]))), "\\."), "[[", 2)
conds[["612SOB"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["SOB"]],le[["T12"]][["SOB"]]))), "\\."), "[[", 2)


conds[["06S"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["S"]],le[["T6"]][["S"]]))), "\\."), "[[", 2)
conds[["012S"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["S"]],le[["T12"]][["S"]]))), "\\."), "[[", 2)
conds[["612S"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["S"]],le[["T12"]][["S"]]))), "\\."), "[[", 2)


conds[["06F"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["F"]],le[["T6"]][["F"]]))), "\\."), "[[", 2)
conds[["012F"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["F"]],le[["T12"]][["F"]]))), "\\."), "[[", 2)
conds[["612F"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["F"]],le[["T12"]][["F"]]))), "\\."), "[[", 2)



x.all[["06FOB"]] <- aldex(ECselex.sub[["06FOB"]], conds[["06FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["012FOB"]] <- aldex(ECselex.sub[["012FOB"]], conds[["012FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["612FOB"]] <- aldex(ECselex.sub[["612FOB"]], conds[["612FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

x.all[["06SOB"]] <- aldex(ECselex.sub[["06SOB"]], conds[["06SOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012SOB"]] <- aldex(ECselex.sub[["012SOB"]], conds[["012SOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612SOB"]] <- aldex(ECselex.sub[["612SOB"]], conds[["612SOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

x.all[["06S"]] <- aldex(ECselex.sub[["06S"]], conds[["06S"]],mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012S"]] <- aldex(ECselex.sub[["012S"]], conds[["012S"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612S"]] <- aldex(ECselex.sub[["612S"]], conds[["612S"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

x.all[["06F"]] <- aldex(ECselex.sub[["06F"]], conds[["06F"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012F"]] <- aldex(ECselex.sub[["012F"]], conds[["012F"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612F"]] <- aldex(ECselex.sub[["612F"]], conds[["612F"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)



pdf(paste(directory,"pdf","FOB_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(x.all[["06FOB"]], type="MW", test="welch")+title("FOB T0 against T6")
aldex.plot(x.all[["012FOB"]], type="MW", test="welch")+title("FOB T0 against T12")
par(mfrow=c(1,1))
aldex.plot(x.all[["612FOB"]], type="MW", test="welch")+title("FOB T6 against T12")
dev.off() 

pdf(paste(directory,"pdf","SOB_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(x.all[["06SOB"]], type="MW", test="welch")+title("SOB T0 against T6")
aldex.plot(x.all[["012SOB"]], type="MW", test="welch")+title("SOB T0 against T12")
par(mfrow=c(1,1))
aldex.plot(x.all[["612SOB"]], type="MW", test="welch")+title("SOB T6 against T12")
dev.off() 

pdf(paste(directory,"pdf","S_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(x.all[["06S"]], type="MW", test="welch")+title("S T0 against T6")
aldex.plot(x.all[["012S"]], type="MW", test="welch")+title("S T0 against T12")
par(mfrow=c(1,1))
aldex.plot(x.all[["612S"]], type="MW", test="welch")+title("S T6 against T12")
dev.off() 

pdf(paste(directory,"pdf","F_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(x.all[["06F"]], type="MW", test="welch")+title("F T0 against T6")
aldex.plot(x.all[["012F"]], type="MW", test="welch")+title("F T0 against T12")
par(mfrow=c(1,1))
aldex.plot(x.all[["612F"]], type="MW", test="welch")+title("F T6 against T12")

dev.off() 
 

########################################################
########################################################
##### KO ANALYSIS 
##### September 2020 
########################################################
########################################################
names(TABLES)
KO <- cbind(TABLES[["ko_metagenome_S"]],TABLES[["ko_metagenome_SOB"]],TABLES[["ko_metagenome_F"]],TABLES[["ko_metagenome_FOB"]])
summary(KO)
### take the colnames to determine sous groups, from the same cohort or at the same time point. 
Sample_name <-c(sort(unique(sapply(strsplit(as.character(colnames(KO)), "\\."), "[[", 1))))
names(KO)

lk<- list()
for (time in c("T0","T6","T12")){
  lk[[time]] <- list()
  for (statut in c("FOB","SOB","F","S")){
    lk[[time]][[statut]] <- KO[,grep(paste0(statut, ".*", time, "$"),colnames(KO))]
  }
}


KOconds <- list()
KOselex.sub <- list()
KOx.all <- list()




# differential(relative) abundance analysis of proportional data
KOselex.sub[["06FOB"]] <- round(cbind(lk[["T0"]][["FOB"]],lk[["T6"]][["FOB"]]))
KOselex.sub[["012FOB"]]<- round(cbind(lk[["T0"]][["FOB"]],lk[["T12"]][["FOB"]]))
KOselex.sub[["612FOB"]] <- round(cbind(lk[["T6"]][["FOB"]],lk[["T12"]][["FOB"]]))

# differential(relative) abundance analysis of proportional data
KOselex.sub[["06SOB"]] <- round(cbind(lk[["T0"]][["SOB"]],lk[["T6"]][["SOB"]]))
KOselex.sub[["012SOB"]]<- round(cbind(lk[["T0"]][["SOB"]],lk[["T12"]][["SOB"]]))
KOselex.sub[["612SOB"]] <- round(cbind(lk[["T6"]][["SOB"]],lk[["T12"]][["SOB"]]))

# differential(relative) abundance analysis of proportional data
KOselex.sub[["06S"]] <- round(cbind(lk[["T0"]][["S"]],lk[["T6"]][["S"]]))
KOselex.sub[["012S"]]<- round(cbind(lk[["T0"]][["S"]],lk[["T12"]][["S"]]))
KOselex.sub[["612S"]] <- round(cbind(lk[["T6"]][["S"]],lk[["T12"]][["S"]]))

# differential(relative) abundance analysis of proportional data
KOselex.sub[["06F"]] <- round(cbind(lk[["T0"]][["F"]],lk[["T6"]][["F"]]))
KOselex.sub[["012F"]]<- round(cbind(lk[["T0"]][["F"]],lk[["T12"]][["F"]]))
KOselex.sub[["612F"]] <- round(cbind(lk[["T6"]][["F"]],lk[["T12"]][["F"]]))




KOconds[["06FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["FOB"]],lk[["T6"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["012FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["FOB"]],lk[["T12"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["612FOB"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["FOB"]],lk[["T12"]][["FOB"]]))), "\\."), "[[", 2)

KOconds[["06SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["SOB"]],lk[["T6"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["012SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["SOB"]],lk[["T12"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["612SOB"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["SOB"]],lk[["T12"]][["SOB"]]))), "\\."), "[[", 2)


KOconds[["06S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["S"]],lk[["T6"]][["S"]]))), "\\."), "[[", 2)
KOconds[["012S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["S"]],lk[["T12"]][["S"]]))), "\\."), "[[", 2)
KOconds[["612S"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["S"]],lk[["T12"]][["S"]]))), "\\."), "[[", 2)


KOconds[["06F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["F"]],lk[["T6"]][["F"]]))), "\\."), "[[", 2)
KOconds[["012F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["F"]],lk[["T12"]][["F"]]))), "\\."), "[[", 2)
KOconds[["612F"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["F"]],lk[["T12"]][["F"]]))), "\\."), "[[", 2)


KOx.all[["06FOB"]] <- aldex(KOselex.sub[["06FOB"]], KOconds[["06FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["012FOB"]] <- aldex(KOselex.sub[["012FOB"]], KOconds[["012FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["612FOB"]] <- aldex(KOselex.sub[["612FOB"]], KOconds[["612FOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

KOx.all[["06SOB"]] <- aldex(KOselex.sub[["06SOB"]], KOconds[["06SOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012SOB"]] <- aldex(KOselex.sub[["012SOB"]], KOconds[["012SOB"]],mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612SOB"]] <- aldex(KOselex.sub[["612SOB"]], KOconds[["612SOB"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)


KOx.all[["06S"]] <- aldex(KOselex.sub[["06S"]], KOconds[["06S"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012S"]] <- aldex(KOselex.sub[["012S"]], KOconds[["012S"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612S"]] <- aldex(KOselex.sub[["612S"]], KOconds[["612S"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

KOx.all[["06F"]] <- aldex(KOselex.sub[["06F"]], KOconds[["06F"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012F"]] <- aldex(KOselex.sub[["012F"]], KOconds[["012F"]],mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612F"]] <- aldex(KOselex.sub[["612F"]], KOconds[["612F"]], mc.samples=1000, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

pdf(paste(directory,"pdf","FOB_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(KOx.all[["06FOB"]], type="MW", test="welch")+title("FOB T0 against T6")
aldex.plot(KOx.all[["012FOB"]], type="MW", test="welch")+title("FOB T0 against T12")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["612FOB"]], type="MW", test="welch")+title("FOB T6 against T12")
dev.off() 

pdf(paste(directory,"pdf","SOB_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(KOx.all[["06SOB"]], type="MW", test="welch")+title("SOB T0 against T6")
aldex.plot(KOx.all[["012SOB"]], type="MW", test="welch")+title("SOB T0 against T12")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["612SOB"]], type="MW", test="welch")+title("SOB T6 against T12")
dev.off() 

pdf(paste(directory,"pdf","S_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(KOx.all[["06S"]], type="MW", test="welch")+title("S T0 against T6")
aldex.plot(KOx.all[["012S"]], type="MW", test="welch")+title("S T0 against T12")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["612S"]], type="MW", test="welch")+title("S T6 against T12")
dev.off() 


pdf(paste(directory,"pdf","F_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(2,1))
aldex.plot(KOx.all[["06F"]], type="MW", test="welch")+title("F T0 against T6")
aldex.plot(KOx.all[["012F"]], type="MW", test="welch")+title("F T0 against T12")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["612F"]], type="MW", test="welch")+title("F T6 against T12")
dev.off() 


# sigKO <- list()
# for (el in names(KOx.all)){
#   sigKO[[el]] <- KOx.all[[el]][which(KOx.all[[el]]$we.eBH<0.01),]
# }

# sigEC <- list()
# for (el in names(x.all)){
#   sigEC[[el]] <- x.all[[el]][which(x.all[[el]]$we.eBH<0.01),]
# }


# for (el in names(sigKO)){
#   if (nrow(sigKO[[el]])!=0){
#   write.csv2(sigKO[el] , paste(directory,"/KO/","KO_",el,"_BH0.01.csv",sep="")) 
#   }
# }

# for (el in names(sigEC)){
#   if (nrow(sigEC[[el]])!=0){
#   write.csv2(sigEC[el] , paste(directory,"/EC/","EC_",el,"_BH0.01.csv",sep="")) 
#   }
# }
