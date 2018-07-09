#summer project
#essential info to put in report

library(Rcpp)
library("readstata13")
library("Hmisc")
####function to extract a record####
reveal <- function(x){
  View(liteumk[liteumk$idno_original==x,])
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

#datacleaning
liteumk <- liteumk[which(!is.na(liteumk$`_d`)),] #remove missing from `_d` [already described in report]

#essential strats here
####Merging 2 datasets####
#liteumk & CODdataset
colnames(CODdataset)[1] <- "idno_original" #change to match names
CODdataset <- as.data.frame(CODdataset)
dim(CODdataset) #10171 x2
dim(liteumk) #2223241 x36

dat <- merge(liteumk,CODdataset, by="idno_original", all.x = TRUE)

dim(dat) #2223241 x37
sum(!is.na(dat$COD)) #147293

#further cleaning done here
#the aim is to describe only to describe the final (FINAL) dataset
#the last state of each individual when last seen???
dat <- dat[which(!is.na(dat$COD)),]

#subsetting the last records

####Descriptive statistics####
varDes <- c('sex','residence','hadva','failure','agegrp',
            'hivstatus_detail','hivstatus_broad','allinfo_treat_pyramid',
            'hivtreat','hivevertreat','preg_at_art_start','agelastseen',
            'art_available','art_avail_cat')
#discarded
#,'age'
des <- list()

for(i in 1:length(varDes)){
  des[[i]] <- table(dat[,which(colnames(dat)==varDes[i])], exclude = NULL)
}

