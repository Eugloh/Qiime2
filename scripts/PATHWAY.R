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

directory <- "./output/merged/picrust2" # directory to set here 
#setwd(directory)

option_list = list(
  make_option(c("-p", "--pattern"), type="character", default="tsv", help="pattern for file names [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sampleFiles  <- grep(opt$pattern,list.files(directory),value=TRUE) # find names

TABLES <- list()
for (el in sampleFiles){
  el_name <- sapply(strsplit(as.character(el), "\\."), "[[", 1)
  TABLES[[el_name]] <- read.delim(paste(dir,el, sep = ""), row.names=1)
}

 names(TABLES) ## retrieve the table of interest with TABLES[["<name>"]]

EC <- cbind(TABLES[["ec_metagenome"]],TABLES[["OB-S_EC_feature-table"]],TABLES[["TW-F_EC_feature-table"]],TABLES[["TW-S_EC_feature-table"]])
summary(EC)
### take the colnames to determine sous groups, from the same cohort or at the same time point. 
Sample_name <-c(sort(unique(sapply(strsplit(as.character(colnames(EC)), "\\."), "[[", 1))))

### the information is sorted into lists, I create table of specifics sous groups. 

# creation of the EC tables, concat in a le list. 


le<- list()
for (time in c("T0","T6","T12","T24")){
  le[[time]] <- list()
  for (statut in c("FOB","SOB","F","S")){
    le[[time]][[statut]] <- EC[,grep(paste0(statut, ".*", time, "$"),colnames(EC))]
  }
}


## visualization and statistics
# ALDEx2
# modelling the data as a log-ratio transformed probability distribution rather than as counts

conds <- list()
KOselex.sub <- list()
x.all <- list()

# differential(relative) abundance analysis of proportional data
KOselex.sub[["06FOB"]] <- round(cbind(le[["T0"]][["FOB"]],le[["T6"]][["FOB"]]))
KOselex.sub[["012FOB"]]<- round(cbind(le[["T0"]][["FOB"]],le[["T12"]][["FOB"]]))
KOselex.sub[["024FOB"]]  <- round(cbind(le[["T0"]][["FOB"]],le[["T24"]][["FOB"]]))
KOselex.sub[["612FOB"]] <- round(cbind(le[["T6"]][["FOB"]],le[["T12"]][["FOB"]]))
KOselex.sub[["624FOB"]] <- round(cbind(le[["T6"]][["FOB"]],le[["T24"]][["FOB"]]))
KOselex.sub[["1224FOB"]] <- round(cbind(le[["T12"]][["FOB"]],le[["T24"]][["FOB"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06SOB"]] <- round(cbind(le[["T0"]][["SOB"]],le[["T6"]][["SOB"]]))
KOselex.sub[["012SOB"]]<- round(cbind(le[["T0"]][["SOB"]],le[["T12"]][["SOB"]]))
KOselex.sub[["024SOB"]]  <- round(cbind(le[["T0"]][["SOB"]],le[["T24"]][["SOB"]]))
KOselex.sub[["612SOB"]] <- round(cbind(le[["T6"]][["SOB"]],le[["T12"]][["SOB"]]))
KOselex.sub[["624SOB"]] <- round(cbind(le[["T6"]][["SOB"]],le[["T24"]][["SOB"]]))
KOselex.sub[["1224SOB"]] <- round(cbind(le[["T12"]][["SOB"]],le[["T24"]][["SOB"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06S"]] <- round(cbind(le[["T0"]][["S"]],le[["T6"]][["S"]]))
KOselex.sub[["012S"]]<- round(cbind(le[["T0"]][["S"]],le[["T12"]][["S"]]))
KOselex.sub[["024S"]]  <- round(cbind(le[["T0"]][["S"]],le[["T24"]][["S"]]))
KOselex.sub[["612S"]] <- round(cbind(le[["T6"]][["S"]],le[["T12"]][["S"]]))
KOselex.sub[["624S"]] <- round(cbind(le[["T6"]][["S"]],le[["T24"]][["S"]]))
KOselex.sub[["1224S"]] <- round(cbind(le[["T12"]][["S"]],le[["T24"]][["S"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06F"]] <- round(cbind(le[["T0"]][["F"]],le[["T6"]][["F"]]))
KOselex.sub[["012F"]]<- round(cbind(le[["T0"]][["F"]],le[["T12"]][["F"]]))
KOselex.sub[["024F"]]  <- round(cbind(le[["T0"]][["F"]],le[["T24"]][["F"]]))
KOselex.sub[["612F"]] <- round(cbind(le[["T6"]][["F"]],le[["T12"]][["F"]]))
KOselex.sub[["624F"]] <- round(cbind(le[["T6"]][["F"]],le[["T24"]][["F"]]))
KOselex.sub[["1224F"]] <- round(cbind(le[["T12"]][["F"]],le[["T24"]][["F"]]))



conds[["06FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["FOB"]],le[["T6"]][["FOB"]]))), "\\."), "[[", 2)
conds[["012FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["FOB"]],le[["T12"]][["FOB"]]))), "\\."), "[[", 2)
conds[["024FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["FOB"]],le[["T24"]][["FOB"]]))), "\\."), "[[", 2)
conds[["612FOB"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["FOB"]],le[["T12"]][["FOB"]]))), "\\."), "[[", 2)
conds[["624FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T6"]][["FOB"]],le[["T24"]][["FOB"]]))), "\\."), "[[", 2)
conds[["1224FOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T12"]][["FOB"]],le[["T24"]][["FOB"]]))), "\\."), "[[", 2)

conds[["06SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["SOB"]],le[["T6"]][["SOB"]]))), "\\."), "[[", 2)
conds[["012SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["SOB"]],le[["T12"]][["SOB"]]))), "\\."), "[[", 2)
conds[["024SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["SOB"]],le[["T24"]][["SOB"]]))), "\\."), "[[", 2)
conds[["612SOB"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["SOB"]],le[["T12"]][["SOB"]]))), "\\."), "[[", 2)
conds[["624SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T6"]][["SOB"]],le[["T24"]][["SOB"]]))), "\\."), "[[", 2)
conds[["1224SOB"]] <-sapply(strsplit(as.character(names(cbind(le[["T12"]][["SOB"]],le[["T24"]][["SOB"]]))), "\\."), "[[", 2)


conds[["06S"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["S"]],le[["T6"]][["S"]]))), "\\."), "[[", 2)
conds[["012S"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["S"]],le[["T12"]][["S"]]))), "\\."), "[[", 2)
conds[["024S"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["S"]],le[["T24"]][["S"]]))), "\\."), "[[", 2)
conds[["612S"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["S"]],le[["T12"]][["S"]]))), "\\."), "[[", 2)
conds[["624S"]] <-sapply(strsplit(as.character(names(cbind(le[["T6"]][["S"]],le[["T24"]][["S"]]))), "\\."), "[[", 2)
conds[["1224S"]] <-sapply(strsplit(as.character(names(cbind(le[["T12"]][["S"]],le[["T24"]][["S"]]))), "\\."), "[[", 2)


conds[["06F"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["F"]],le[["T6"]][["F"]]))), "\\."), "[[", 2)
conds[["012F"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["F"]],le[["T12"]][["F"]]))), "\\."), "[[", 2)
conds[["024F"]] <-sapply(strsplit(as.character(names(cbind(le[["T0"]][["F"]],le[["T24"]][["F"]]))), "\\."), "[[", 2)
conds[["612F"]]<-sapply(strsplit(as.character(names(cbind(le[["T6"]][["F"]],le[["T12"]][["F"]]))), "\\."), "[[", 2)
conds[["624F"]] <-sapply(strsplit(as.character(names(cbind(le[["T6"]][["F"]],le[["T24"]][["F"]]))), "\\."), "[[", 2)
conds[["1224F"]] <-sapply(strsplit(as.character(names(cbind(le[["T12"]][["F"]],le[["T24"]][["F"]]))), "\\."), "[[", 2)



x.all[["06FOB"]] <- aldex(selex.sub[["06FOB"]], conds[["06FOB"]], mc.samples=length(conds[["06FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["012FOB"]] <- aldex(selex.sub[["012FOB"]], conds[["012FOB"]], mc.samples=length(conds[["012FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["024FOB"]] <- aldex(selex.sub[["024FOB"]], conds[["024FOB"]], mc.samples=length(conds[["024FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612FOB"]] <- aldex(selex.sub[["612FOB"]], conds[["612FOB"]], mc.samples=length(conds[["612FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["624FOB"]] <- aldex(selex.sub[["624FOB"]], conds[["624FOB"]], mc.samples=length(conds[["624FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
x.all[["1224FOB"]] <- aldex(selex.sub[["1224FOB"]], conds[["1224FOB"]], mc.samples=length(conds[["1224FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

x.all[["06SOB"]] <- aldex(selex.sub[["06SOB"]], conds[["06SOB"]], mc.samples=length(conds[["06SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012SOB"]] <- aldex(selex.sub[["012SOB"]], conds[["012SOB"]], mc.samples=length(conds[["012SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["024SOB"]] <- aldex(selex.sub[["024SOB"]], conds[["024SOB"]], mc.samples=length(conds[["024SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612SOB"]] <- aldex(selex.sub[["612SOB"]], conds[["612SOB"]], mc.samples=length(conds[["612SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["624SOB"]] <- aldex(selex.sub[["624SOB"]], conds[["624SOB"]], mc.samples=length(conds[["624SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["1224SOB"]] <- aldex(selex.sub[["1224SOB"]], conds[["1224SOB"]], mc.samples=length(conds[["1224SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

x.all[["06S"]] <- aldex(selex.sub[["06S"]], conds[["06S"]], mc.samples=length(conds[["06S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012S"]] <- aldex(selex.sub[["012S"]], conds[["012S"]], mc.samples=length(conds[["012S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["024S"]] <- aldex(selex.sub[["024S"]], conds[["024S"]], mc.samples=length(conds[["024S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612S"]] <- aldex(selex.sub[["612S"]], conds[["612S"]], mc.samples=length(conds[["612S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["624S"]] <- aldex(selex.sub[["624S"]], conds[["624S"]], mc.samples=length(conds[["624S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["1224S"]] <- aldex(selex.sub[["1224S"]], conds[["1224S"]], mc.samples=length(conds[["1224S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

x.all[["06F"]] <- aldex(selex.sub[["06F"]], conds[["06F"]], mc.samples=length(conds[["06F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["012F"]] <- aldex(selex.sub[["012F"]], conds[["012F"]], mc.samples=length(conds[["012F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["024F"]] <- aldex(selex.sub[["024F"]], conds[["024F"]], mc.samples=length(conds[["024F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["612F"]] <- aldex(selex.sub[["612F"]], conds[["612F"]], mc.samples=length(conds[["612F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["624F"]] <- aldex(selex.sub[["624F"]], conds[["624F"]], mc.samples=length(conds[["624F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
x.all[["1224F"]] <- aldex(selex.sub[["1224F"]], conds[["1224F"]], mc.samples=length(conds[["1224F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)



pdf(paste(directory,"pdf","FOB_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(x.all[["06FOB"]], type="MW", test="welch")+title("FOB T0 against T6")
aldex.plot(x.all[["012FOB"]], type="MW", test="welch")+title("FOB T0 against T12")
aldex.plot(x.all[["024FOB"]], type="MW", test="welch")+title("FOB T0 against T24")
par(mfrow=c(2,1))
aldex.plot(x.all[["612FOB"]], type="MW", test="welch")+title("FOB T6 against T12")
aldex.plot(x.all[["624FOB"]], type="MW", test="welch")+title("FOB T6 against T24")
par(mfrow=c(1,1))
aldex.plot(x.all[["1224FOB"]], type="MW", test="welch")+title("FOB T12 against T24")
dev.off() 

pdf(paste(directory,"pdf","SOB_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(x.all[["06SOB"]], type="MW", test="welch")+title("SOB T0 against T6")
aldex.plot(x.all[["012SOB"]], type="MW", test="welch")+title("SOB T0 against T12")
aldex.plot(x.all[["024SOB"]], type="MW", test="welch")+title("SOB T0 against T24")
par(mfrow=c(2,1))
aldex.plot(x.all[["612SOB"]], type="MW", test="welch")+title("SOB T6 against T12")
aldex.plot(x.all[["624SOB"]], type="MW", test="welch")+title("SOB T6 against T24")
par(mfrow=c(1,1))
aldex.plot(x.all[["1224SOB"]], type="MW", test="welch")+title("SOB T12 against T24")
dev.off() 

pdf(paste(directory,"pdf","S_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(x.all[["06S"]], type="MW", test="welch")+title("S T0 against T6")
aldex.plot(x.all[["012S"]], type="MW", test="welch")+title("S T0 against T12")
aldex.plot(x.all[["024S"]], type="MW", test="welch")+title("S T0 against T24")
par(mfrow=c(2,1))
aldex.plot(x.all[["612S"]], type="MW", test="welch")+title("S T6 against T12")
aldex.plot(x.all[["624S"]], type="MW", test="welch")+title("S T6 against T24")
par(mfrow=c(1,1))
aldex.plot(x.all[["1224S"]], type="MW", test="welch")+title("S T12 against T24")
dev.off() 

pdf(paste(directory,"pdf","F_EC_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(x.all[["06F"]], type="MW", test="welch")+title("F T0 against T6")
aldex.plot(x.all[["012F"]], type="MW", test="welch")+title("F T0 against T12")
aldex.plot(x.all[["024F"]], type="MW", test="welch")+title("F T0 against T24")
par(mfrow=c(2,1))
aldex.plot(x.all[["612F"]], type="MW", test="welch")+title("F T6 against T12")
aldex.plot(x.all[["624F"]], type="MW", test="welch")+title("F T6 against T24")
par(mfrow=c(1,1))
aldex.plot(x.all[["1224F"]], type="MW", test="welch")+title("F T12 against T24")
dev.off() 
 
library(dplyr)

test <- x.all[["06FOB"]][which(x.all[["06FOB"]]$we.eBH<0.05),] 
test1 <- x.all[["012FOB"]][which(x.all[["012FOB"]]$we.eBH<0.05),]
test2 <- x.all[["024FOB"]][which(x.all[["024FOB"]]$we.eBH<0.05),]
test10 <- x.all[["06S"]][which(x.all[["012S"]]$we.eBH<0.05),]
test4 <- x.all[["012S"]][which(x.all[["012S"]]$we.eBH<0.05),]
test5 <- x.all[["024S"]][which(x.all[["024S"]]$we.eBH<0.05),]
test6 <- x.all[["624S"]][which(x.all[["624S"]]$we.eBH<0.05),]

test7 <- x.all[["012F"]][which(x.all[["012F"]]$we.eBH<0.05),]
test8 <- x.all[["024F"]][which(x.all[["024F"]]$we.eBH<0.05),]
test9 <- x.all[["624F"]][which(x.all[["624F"]]$we.eBH<0.05),]

setwd("/media/eugloh/CL/DAVID")
dataframes2xls::write.xls(c(test,test1,test2,test10,
                            test4,test5,
                            test6,test7,test8,
                            test9), 
                          row.names = TRUE,col.names = TRUE, "table_weeBH05.xls")

########################################################
########################################################
##### KO ANALYSIS 
##### September 2020 
########################################################
########################################################
names(TABLES)
KO <- cbind(TABLES[["ko_metagenome"]],TABLES[["OB_S_KO_feature-table"]],TABLES[["TW-F_KO_feature-table"]],TABLES[["TW-S_KO_feature-table"]])
summary(KO)
### take the colnames to determine sous groups, from the same cohort or at the same time point. 
Sample_name <-c(sort(unique(sapply(strsplit(as.character(colnames(KO)), "\\."), "[[", 1))))
names(KO)

lk<- list()
for (time in c("T0","T6","T12","T24")){
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
KOselex.sub[["024FOB"]]  <- round(cbind(lk[["T0"]][["FOB"]],lk[["T24"]][["FOB"]]))
KOselex.sub[["612FOB"]] <- round(cbind(lk[["T6"]][["FOB"]],lk[["T12"]][["FOB"]]))
KOselex.sub[["624FOB"]] <- round(cbind(lk[["T6"]][["FOB"]],lk[["T24"]][["FOB"]]))
KOselex.sub[["1224FOB"]] <- round(cbind(lk[["T12"]][["FOB"]],lk[["T24"]][["FOB"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06SOB"]] <- round(cbind(lk[["T0"]][["SOB"]],lk[["T6"]][["SOB"]]))
KOselex.sub[["012SOB"]]<- round(cbind(lk[["T0"]][["SOB"]],lk[["T12"]][["SOB"]]))
KOselex.sub[["024SOB"]]  <- round(cbind(lk[["T0"]][["SOB"]],lk[["T24"]][["SOB"]]))
KOselex.sub[["612SOB"]] <- round(cbind(lk[["T6"]][["SOB"]],lk[["T12"]][["SOB"]]))
KOselex.sub[["624SOB"]] <- round(cbind(lk[["T6"]][["SOB"]],lk[["T24"]][["SOB"]]))
KOselex.sub[["1224SOB"]] <- round(cbind(lk[["T12"]][["SOB"]],lk[["T24"]][["SOB"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06S"]] <- round(cbind(lk[["T0"]][["S"]],lk[["T6"]][["S"]]))
KOselex.sub[["012S"]]<- round(cbind(lk[["T0"]][["S"]],lk[["T12"]][["S"]]))
KOselex.sub[["024S"]]  <- round(cbind(lk[["T0"]][["S"]],lk[["T24"]][["S"]]))
KOselex.sub[["612S"]] <- round(cbind(lk[["T6"]][["S"]],lk[["T12"]][["S"]]))
KOselex.sub[["624S"]] <- round(cbind(lk[["T6"]][["S"]],lk[["T24"]][["S"]]))
KOselex.sub[["1224S"]] <- round(cbind(lk[["T12"]][["S"]],lk[["T24"]][["S"]]))
# differential(relative) abundance analysis of proportional data
KOselex.sub[["06F"]] <- round(cbind(lk[["T0"]][["F"]],lk[["T6"]][["F"]]))
KOselex.sub[["012F"]]<- round(cbind(lk[["T0"]][["F"]],lk[["T12"]][["F"]]))
KOselex.sub[["024F"]]  <- round(cbind(lk[["T0"]][["F"]],lk[["T24"]][["F"]]))
KOselex.sub[["612F"]] <- round(cbind(lk[["T6"]][["F"]],lk[["T12"]][["F"]]))
KOselex.sub[["624F"]] <- round(cbind(lk[["T6"]][["F"]],lk[["T24"]][["F"]]))
KOselex.sub[["1224F"]] <- round(cbind(lk[["T12"]][["F"]],lk[["T24"]][["F"]]))



KOconds[["06FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["FOB"]],lk[["T6"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["012FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["FOB"]],lk[["T12"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["024FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["FOB"]],lk[["T24"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["612FOB"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["FOB"]],lk[["T12"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["624FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["FOB"]],lk[["T24"]][["FOB"]]))), "\\."), "[[", 2)
KOconds[["1224FOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T12"]][["FOB"]],lk[["T24"]][["FOB"]]))), "\\."), "[[", 2)

KOconds[["06SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["SOB"]],lk[["T6"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["012SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["SOB"]],lk[["T12"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["024SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["SOB"]],lk[["T24"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["612SOB"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["SOB"]],lk[["T12"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["624SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["SOB"]],lk[["T24"]][["SOB"]]))), "\\."), "[[", 2)
KOconds[["1224SOB"]] <-sapply(strsplit(as.character(names(cbind(lk[["T12"]][["SOB"]],lk[["T24"]][["SOB"]]))), "\\."), "[[", 2)


KOconds[["06S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["S"]],lk[["T6"]][["S"]]))), "\\."), "[[", 2)
KOconds[["012S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["S"]],lk[["T12"]][["S"]]))), "\\."), "[[", 2)
KOconds[["024S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["S"]],lk[["T24"]][["S"]]))), "\\."), "[[", 2)
KOconds[["612S"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["S"]],lk[["T12"]][["S"]]))), "\\."), "[[", 2)
KOconds[["624S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["S"]],lk[["T24"]][["S"]]))), "\\."), "[[", 2)
KOconds[["1224S"]] <-sapply(strsplit(as.character(names(cbind(lk[["T12"]][["S"]],lk[["T24"]][["S"]]))), "\\."), "[[", 2)


KOconds[["06F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["F"]],lk[["T6"]][["F"]]))), "\\."), "[[", 2)
KOconds[["012F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["F"]],lk[["T12"]][["F"]]))), "\\."), "[[", 2)
KOconds[["024F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T0"]][["F"]],lk[["T24"]][["F"]]))), "\\."), "[[", 2)
KOconds[["612F"]]<-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["F"]],lk[["T12"]][["F"]]))), "\\."), "[[", 2)
KOconds[["624F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T6"]][["F"]],lk[["T24"]][["F"]]))), "\\."), "[[", 2)
KOconds[["1224F"]] <-sapply(strsplit(as.character(names(cbind(lk[["T12"]][["F"]],lk[["T24"]][["F"]]))), "\\."), "[[", 2)


KOx.all[["06FOB"]] <- aldex(KOselex.sub[["06FOB"]], KOconds[["06FOB"]], mc.samples=length(KOconds[["06FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["012FOB"]] <- aldex(KOselex.sub[["012FOB"]], KOconds[["012FOB"]], mc.samples=length(KOconds[["012FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["024FOB"]] <- aldex(KOselex.sub[["024FOB"]], KOconds[["024FOB"]], mc.samples=length(KOconds[["024FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612FOB"]] <- aldex(KOselex.sub[["612FOB"]], KOconds[["612FOB"]], mc.samples=length(KOconds[["612FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["624FOB"]] <- aldex(KOselex.sub[["624FOB"]], KOconds[["624FOB"]], mc.samples=length(KOconds[["624FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
KOx.all[["1224FOB"]] <- aldex(KOselex.sub[["1224FOB"]], KOconds[["1224FOB"]], mc.samples=length(KOconds[["1224FOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

KOx.all[["06SOB"]] <- aldex(KOselex.sub[["06SOB"]], KOconds[["06SOB"]], mc.samples=length(KOconds[["06SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012SOB"]] <- aldex(KOselex.sub[["012SOB"]], KOconds[["012SOB"]], mc.samples=length(KOconds[["012SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["024SOB"]] <- aldex(KOselex.sub[["024SOB"]], KOconds[["024SOB"]], mc.samples=length(KOconds[["024SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612SOB"]] <- aldex(KOselex.sub[["612SOB"]], KOconds[["612SOB"]], mc.samples=length(KOconds[["612SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["624SOB"]] <- aldex(KOselex.sub[["624SOB"]], KOconds[["624SOB"]], mc.samples=length(KOconds[["624SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["1224SOB"]] <- aldex(KOselex.sub[["1224SOB"]], KOconds[["1224SOB"]], mc.samples=length(KOconds[["1224SOB"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)


KOx.all[["06S"]] <- aldex(KOselex.sub[["06S"]], KOconds[["06S"]], mc.samples=length(KOconds[["06S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012S"]] <- aldex(KOselex.sub[["012S"]], KOconds[["012S"]], mc.samples=length(KOconds[["012S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["024S"]] <- aldex(KOselex.sub[["024S"]], KOconds[["024S"]], mc.samples=length(KOconds[["024S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612S"]] <- aldex(KOselex.sub[["612S"]], KOconds[["612S"]], mc.samples=length(KOconds[["612S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["624S"]] <- aldex(KOselex.sub[["624S"]], KOconds[["624S"]], mc.samples=length(KOconds[["624S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["1224S"]] <- aldex(KOselex.sub[["1224S"]], KOconds[["1224S"]], mc.samples=length(KOconds[["1224S"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

KOx.all[["06F"]] <- aldex(KOselex.sub[["06F"]], KOconds[["06F"]], mc.samples=length(KOconds[["06F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["012F"]] <- aldex(KOselex.sub[["012F"]], KOconds[["012F"]], mc.samples=length(KOconds[["012F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["024F"]] <- aldex(KOselex.sub[["024F"]], KOconds[["024F"]], mc.samples=length(KOconds[["024F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["612F"]] <- aldex(KOselex.sub[["612F"]], KOconds[["612F"]], mc.samples=length(KOconds[["612F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["624F"]] <- aldex(KOselex.sub[["624F"]], KOconds[["624F"]], mc.samples=length(KOconds[["624F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)
KOx.all[["1224F"]] <- aldex(KOselex.sub[["1224F"]], KOconds[["1224F"]], mc.samples=length(KOconds[["1224F"]]), test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=TRUE)

pdf(paste(directory,"pdf","FOB_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(KOx.all[["06FOB"]], type="MW", test="welch")+title("FOB T0 against T6")
aldex.plot(KOx.all[["012FOB"]], type="MW", test="welch")+title("FOB T0 against T12")
aldex.plot(KOx.all[["024FOB"]], type="MW", test="welch")+title("FOB T0 against T24")
par(mfrow=c(2,1))
aldex.plot(KOx.all[["612FOB"]], type="MW", test="welch")+title("FOB T6 against T12")
aldex.plot(KOx.all[["624FOB"]], type="MW", test="welch")+title("FOB T6 against T24")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["1224FOB"]], type="MW", test="welch")+title("FOB T12 against T24")
dev.off() 

pdf(paste(directory,"pdf","SOB_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(KOx.all[["06SOB"]], type="MW", test="welch")+title("SOB T0 against T6")
aldex.plot(KOx.all[["012SOB"]], type="MW", test="welch")+title("SOB T0 against T12")
aldex.plot(KOx.all[["024SOB"]], type="MW", test="welch")+title("SOB T0 against T24")
par(mfrow=c(2,1))
aldex.plot(KOx.all[["612SOB"]], type="MW", test="welch")+title("SOB T6 against T12")
aldex.plot(KOx.all[["624SOB"]], type="MW", test="welch")+title("SOB T6 against T24")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["1224SOB"]], type="MW", test="welch")+title("SOB T12 against T24")
dev.off() 

pdf(paste(directory,"pdf","S_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(KOx.all[["06S"]], type="MW", test="welch")+title("S T0 against T6")
aldex.plot(KOx.all[["012S"]], type="MW", test="welch")+title("S T0 against T12")
aldex.plot(KOx.all[["024S"]], type="MW", test="welch")+title("S T0 against T24")
par(mfrow=c(2,1))
aldex.plot(KOx.all[["612S"]], type="MW", test="welch")+title("S T6 against T12")
aldex.plot(KOx.all[["624S"]], type="MW", test="welch")+title("S T6 against T24")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["1224S"]], type="MW", test="welch")+title("S T12 against T24")
dev.off() 


pdf(paste(directory,"pdf","F_KO_welch_test_timepoints_comparation.pdf",sep="/"))
par(mfrow=c(3,1))
aldex.plot(KOx.all[["06F"]], type="MW", test="welch")+title("F T0 against T6")
aldex.plot(KOx.all[["012F"]], type="MW", test="welch")+title("F T0 against T12")
aldex.plot(KOx.all[["024F"]], type="MW", test="welch")+title("F T0 against T24")
par(mfrow=c(2,1))
aldex.plot(KOx.all[["612F"]], type="MW", test="welch")+title("F T6 against T12")
aldex.plot(KOx.all[["624F"]], type="MW", test="welch")+title("F T6 against T24")
par(mfrow=c(1,1))
aldex.plot(KOx.all[["1224F"]], type="MW", test="welch")+title("F T12 against T24")
dev.off() 

test <- KOx.all[["06FOB"]][which([["06FOB"]]$we.eBH<0.05),] 
test1 <- KOx.all[["012FOB"]][which(KOx.all[["012FOB"]]$we.eBH<0.05),]
test2 <- KOx.all[["024FOB"]][which(KOx.all[["024FOB"]]$we.eBH<0.05),]
test10 <- KOx.all[["06S"]][which(KOx.all[["012S"]]$we.eBH<0.05),]
test4 <- KOx.all[["012S"]][which(KOx.all[["012S"]]$we.eBH<0.05),]
test5 <- KOx.all[["024S"]][which(KOx.all[["024S"]]$we.eBH<0.05),]
test6 <- KOx.all[["624S"]][which(KOx.all[["624S"]]$we.eBH<0.05),]

test7 <- KOx.all[["012F"]][which(KOx.all[["012F"]]$we.eBH<0.05),]
test8 <- KOx.all[["024F"]][which(KOx.all[["024F"]]$we.eBH<0.05),]
test9 <- KOx.all[["624F"]][which(KOx.all[["624F"]]$we.eBH<0.05),]


sigKO <- list()
for (el in names(KOx.all)){
  sigKO[[el]] <- KOx.all[[el]][which(KOx.all[[el]]$we.eBH<0.01),]
}

sigEC <- list()
for (el in names(x.all)){
  sigEC[[el]] <- x.all[[el]][which(x.all[[el]]$we.eBH<0.01),]
}


for (el in names(sigKO)){
  if (nrow(sigKO[[el]])!=0){
  write.csv2(sigKO[el] , paste(directory,"/KO/","KO_",el,"_BH0.01.csv",sep="")) 
  }
}

for (el in names(sigEC)){
  if (nrow(sigEC[[el]])!=0){
  write.csv2(sigEC[el] , paste(directory,"/EC/","EC_",el,"_BH0.01.csv",sep="")) 
  }
}