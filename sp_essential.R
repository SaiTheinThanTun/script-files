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

####2 FAIL variables####
table(liteumk$fail, liteumk$failure, exclude=NULL)
excludedDeaths <- liteumk[is.na(liteumk$fail) & liteumk$failure==1,]
excludedDeathsID <- excludedDeaths$idno_original
#nrow(liteumk[liteumk$idno_original %in% excludedDeathsID,])
excludedDeathsAllEpisodes <- (liteumk[liteumk$idno_original %in% excludedDeathsID,])
table(excludedDeathsAllEpisodes$idno_original)
table(excludedDeathsAllEpisodes$idno_original, excludedDeathsAllEpisodes$hadva)

CODdataset[CODdataset[,1] %in% excludedDeathsID,] 
#no COD data for excluded cases. so it's ok to proceed

####fail vs HIV####
table(liteumk$hivstatus_detail, liteumk$fail, exclude=NULL)

####no. of participants####
#114980
length(unique(liteumk$idno_original))

#no. of deaths: 15972
table(liteumk$`_d`, liteumk$fail, exclude=NULL)
#table(liteumk$`_d`, liteumk$failure, exclude=NULL) #failure var removed

#no. of deaths which had va (main dataset)
table(liteumk$fail, liteumk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

####checking cases that failed (failure==1)####
head(liteumk[liteumk$failure==1 & !is.na(liteumk$failure),])
case15 <- liteumk[liteumk$idno_original==15,]
case12 <- liteumk[liteumk$idno_original==12,]
case108 <- liteumk[liteumk$idno_original==108,]
case13 <- liteumk[liteumk$idno_original==13,]
case61051 <- liteumk[liteumk$idno_original==61051,]

####checking cases that didn't fail at all####
#sum of failure with each unique case
mini_liteumk <- liteumk[c(1:1000),]
by(mini_liteumk$failure,mini_liteumk$idno_original,function(x){sum(x,na.rm = T)})
by(mini_liteumk$fail,mini_liteumk$idno_original,function(x){sum(x)})

####checking if there are more than 1 fail in each####
sum(c(by(liteumk$fail,liteumk$idno_original,function(x){sum(x)}))>1, na.rm = T)

####cases that failed but didn't have VA####
head(liteumk[liteumk$failure==1 & !is.na(liteumk$failure) & liteumk$hadva==0,])

####idno vs idno_original#####
#they are unique in their own ways :)
length(unique(liteumk$idno))
length(unique(liteumk$idno_original))

####allinfo_treat_pyramid vs hivtreat####
table(liteumk$allinfo_treat_pyramid, liteumk$hivtreat)

####allinfo_treat_pyramid vs hivevertreat####
table(liteumk$allinfo_treat_pyramid, liteumk$hivevertreat)

####battle of 2 sexes####
#table(liteumk$sex, liteumk$male, exclude=NULL) #correct


####Merging 2 datasets####
#liteumk & CODdataset
