#summer project
#essential info to put in report

library(Rcpp)
library("readstata13")
library("Hmisc")
library(plyr)
####function to extract a record####
reveal <- function(x){
  View(liteumk[liteumk$idno_original==x,])
}
####function to identify frequentest value####
tableR <- function(x){
  #function to identify the most frequent 'residence' value
  #taking the earliest one if there are 2 values with equal frequency
  unique_r <- unique(x)
  tr <- rbind(label=unique_r, count=sapply(unique_r,function(y)sum(x==y, na.rm = T)))
  tr[1,which.max(tr[2,])]
}

setwd("~/OneDrive/Summer Project/data/")

#Read lite data file
liteumk <- readRDS("liteumk.RDS")

#InSilicoVA
CODdataset <- readRDS("CODdataset.RDS")

####consistency checks####
####no. of participants####
#114980
length(unique(liteumk$idno_original))

#no. of deaths: 15972
table(liteumk$`_d`, exclude=NULL)
#table(liteumk$`_d`, liteumk$failure, exclude=NULL) #failure var removed
#Because of this resplitting: ‘fail’ variable is now useless as they have been further multiplicated 
#by the split ie. some may have multiple 1 on the ‘fail’ variable

#no. of deaths which had va (main dataset)
table(liteumk$`_d`, liteumk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

#tmp <- liteumk[which(is.na(liteumk$`_d`)),]

#datacleaning
liteumk <- liteumk[which(!is.na(liteumk$`_d`)),] #remove missing from `_d` [already described in report]
length(unique(liteumk[which(is.na(liteumk$`_d`)),]$idno_original)) #2034 cases only instead of 2037
#the discrepancy in # is because some cases have records with NA in `_d` (ie. same entry & exit date)

####Left censoring for 2007 when HIV testing is done for all adults (15+)####
liteumk <- liteumk[which((liteumk$entry>='2007-01-01')),]


#issue with 'residence'
mostFreqRes <- by(liteumk$residence,liteumk$idno_original,function(x){tableR(x)})
freq.id <- table(liteumk$idno_original)

residence2 <- rep(mostFreqRes, freq.id)
liteumk <- cbind(liteumk, residence2)
#tmp5 <- liteumk[is.na(liteumk$residence),] #testing
# reveal(13)
# reveal(14)

#essential strats here
####Merging 2 datasets####
#liteumk & CODdataset
colnames(CODdataset)[1] <- "idno_original" #change to match names
CODdataset <- as.data.frame(CODdataset)
dim(CODdataset) #10171 x2
dim(liteumk) #1271610 x37

dat <- merge(liteumk,CODdataset, by="idno_original", all.x = TRUE)

dim(dat) #1271610 x38
sum(!is.na(dat$COD)) #40315

#further cleaning done here
#the aim is to describe only to describe the final (FINAL) dataset
#the first state of each individual when first enrolled???

#dat <- dat[which(!is.na(dat$COD)),]

#subsetting the first records
dat <- dat[order(dat$entry, decreasing = F),] #first record
#dat <- dat[order(dat$exit, decreasing = T),]
datL <- dat[!duplicated(dat$idno_original),]

####Descriptive statistics of the last records####
varDes <- c('sex','residence','residence2','hadva','failure','agegrp',
            'hivstatus_detail','hivstatus_broad','allinfo_treat_pyramid',
            'hivtreat','hivevertreat','preg_at_art_start',
            'art_available','art_avail_cat', 'COD')
#discarded
#,'age','agelastseen',
des <- list()

for(i in 1:length(varDes)){
  des[[i]] <- table(datL[,which(colnames(datL)==varDes[i])], exclude = NULL)
}
names(des) <- varDes

#check age consistency
#summarize(dat$age,by(dat$agegrp), max)
by(liteumk$age,liteumk$agegrp,summary) #all age groups are consistent
by(liteumk$age,liteumk$agegrp,function(x){sum(is.na(x))})

#residence variable doesn't have data on all episodes
sum(c(by(dat$residence,liteumk$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T) 
#cases with at least 1 missing residence 54855, not 48725
sum(c(by(dat$residence2,dat$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T)
sum(c(by(datL$residence2,datL$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T)
#residence2: 16017, not 14748
#tmp <- datL[which(is.na(datL$residence)),]
#tmp2 <- dat[dat$idno_original %in% tmp$idno_original,]

#do they change their residence during the study period??
# res_umk <- liteumk[!is.na(liteumk$residence),] #removed na!
# sum(c(by(res_umk$residence,res_umk$idno_original,function(x){sum(length(unique(x)))}))>1, na.rm=T) #how many cases have over 1 residency change?
# tmp <- by(res_umk$residence,res_umk$idno_original,function(x){sum(length(unique(x)))}) #is there any residency change? provides count
# 
# tmp2 <- res_umk[which(res_umk$idno_original %in% names(which(tmp>1))),] #subsetting residency change of over 1 time
# by(tmp2$residence,tmp2$idno_original,function(x){sort(table(x), decreasing = T)})
# tmp3 <- by(tmp2$residence,tmp2$idno_original,function(x){tableR(x)})
# freq.id <- table(tmp2$idno_original)
# 
# residence2 <- rep(tmp3, freq.id)
# tmp4 <- cbind(tmp2, residence2)
# #testing
# reveal(173676)
# reveal(173669)
# reveal(170856)
# reveal(167383)
# reveal(165117)
# reveal(15)


