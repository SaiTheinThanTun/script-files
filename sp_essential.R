#summer project
#essential info to put in report

library(Rcpp)
library("readstata13")
library("Hmisc")

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
table(liteumk$`_d`, liteumk$fail, exclude=NULL)
#table(liteumk$`_d`, liteumk$failure, exclude=NULL) #failure var removed
#Because of this resplitting: ‘fail’ variable is now useless as they have been further multiplicated 
#by the split ie. some may have multiple 1 on the ‘fail’ variable

#no. of deaths which had va (main dataset)
table(liteumk$fail, liteumk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

####checking if there are more than 1 fail in each####
sum(c(by(liteumk$fail,liteumk$idno_original,function(x){sum(x)}))>1, na.rm = T)

####cases that failed but didn't have VA####
head(liteumk[liteumk$failure==1 & !is.na(liteumk$failure) & liteumk$hadva==0,])

####allinfo_treat_pyramid vs hivtreat####
table(liteumk$allinfo_treat_pyramid, liteumk$hivtreat)

####allinfo_treat_pyramid vs hivevertreat####
table(liteumk$allinfo_treat_pyramid, liteumk$hivevertreat)

####battle of 2 sexes####
#table(liteumk$sex, liteumk$male, exclude=NULL) #correct


####Merging 2 datasets####
#liteumk & CODdataset

####function to extract a record####
reveal <- function(x){
  View(liteumk[liteumk$idno_original==x,])
}