#summer project

#library(foreign)
library(Rcpp)
library("readstata13")
library("Hmisc")
createNewFiles <- 0 #switch to turn on/off saving new RDS files. ie. if dta file has been modified in stata

setwd("~/OneDrive/Summer Project/data/")

#Read stata data file
if(createNewFiles){
  umk <- read.dta13("alpha_uMkhanyakude-170601.dta")
  saveRDS(umk,"umk.RDS")
}
if(!createNewFiles){
  umk <- readRDS("umk.RDS") 
}

####CodeBook for main####
if(FALSE){
# #variable name
# names.umk <- names(umk)
# #variable label
# lab.umk <- varlabel(umk)
# #variable type 
# type.umk <- sapply(umk, class)
# #variable label values
# labval.umk <- get.label.tables(umk)
# labval.umk <- sapply(lapply(labval.umk,names),function(x) {paste(x,collapse = ", ")})
# codebk <- cbind(names.umk, lab.umk, type.umk, labval.umk)
# #head(codebk)
# #dim(codebk)
# #write.csv(codebk,"codebook_uMk_R.csv")
}

#subsetting only important variables
selectedVar <- read.csv('selected_var.csv')

####lite version of main dataset####
liteumk <- umk[,names(umk) %in% selectedVar$name]
#will loose all the structure of stata data file
if(createNewFiles){
  saveRDS(liteumk,"liteumk.RDS") 
}

###exploratory####
#InSilicoVA
isv <- read.csv("InSilicoVA_hiv_nophys_output_umkhanyakude.csv") #read.csv("InSilicoVA outputs - uMkhanyakude.csv")
#names(isv)
isvID <- isv[,1]
isv <- isv[,-1]
isvCOD <- names(isv)
COD <- NA
for(i in 1:nrow(isv)){
  COD[i] <- isvCOD[which.max(isv[i,])]
}
CODdataset <- cbind(isvID,COD)
#saveRDS(CODdataset,"CODdataset.RDS")

####another way of getting InSilicoVA data####
if(FALSE){
# isvlimited <- isv[,c(42:52)] #42:52 where injury CoDs are
# describe(isvlimited)
# sum(isvlimited) 
# dim(isvlimited)
# 
# injury <- sapply(isvlimited, sum)
# names(injury) <- names(isvlimited)
# 
# par(las=2)
# barplot(injury)
# par(las=0)
# 
# hist(injury)
# #another way: but not right: hist(isvlimited[,1])
# 
# ####are their id unique?####
# sum(duplicated(isv[,1]))
# 
# ####are the code with injuary most probable?####
# 
# 
# #recombining with ID
# isv2 <- cbind(isv[,1],isvlimited) #10171 rows
# 
# ####are the ID with injury reported unique?####
# isv2.result <- isv2[rowSums(isv2[,c(2:12)])>0,] #953 rows
# sum(duplicated(isv2.result[,1])) #all unique
# 
# ####does any participant have multiple type of injury data? eg. having both traffic accident and assult####
# isv2.result[rowSums(isv2.result[,c(2:12)])>1,] #no
# 
# ####are all IDs from InSilicoVA data in the main data?####
# sum((isv2.result[,1] %in% umk$idno_original)) #yes, all 953 are in 
}

####consistency checks####
####if death occurs, "hadva" variable will be filled####
table(liteumk$hadva, liteumk$failure, exclude = NULL)
sum(liteumk$hadva, na.rm=T)
sum(liteumk$`_d`, na.rm=T)
#can't be checked as hadva has 1 for all episodes of each case 
#whereas fail occurs at most 1 per person
#test <- liteumk[liteumk$hadva==1 & is.na(liteumk$failure),]

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
table(liteumk$hivstatus_detail, liteumk$`_d`, exclude=NULL)

####no. of participants####
#114980
length(unique(liteumk$idno_original))

#no. of deaths: 15972
table(liteumk$`_d`, liteumk$fail, exclude=NULL)
#table(liteumk$`_d`, liteumk$failure, exclude=NULL) #failure var removed
#Because of this resplitting: ‘fail’ variable is now useless as they have been further multiplicated 
#by the split ie. some may have multiple 1 on the ‘fail’ variable


#no. of deaths which had va (main dataset)
table(liteumk$`_d`, liteumk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

####checking cases that failed (failure==1)####
head(liteumk[liteumk$failure==1 & !is.na(liteumk$failure),])
case17 <- liteumk[liteumk$idno_original==17,] #last_neg_date and frst_pos_date present
case12 <- liteumk[liteumk$idno_original==12,] #age check for transtition to above 60
case108 <- liteumk[liteumk$idno_original==108,]
case13 <- liteumk[liteumk$idno_original==13,]
case61051 <- liteumk[liteumk$idno_original==61051,]
case185 <- liteumk[liteumk$idno_original==185,] #age>90

####checking cases that didn't fail at all####
#sum of failure with each unique case
mini_liteumk <- liteumk[c(1:1000),]
by(mini_liteumk$failure,mini_liteumk$idno_original,function(x){sum(x,na.rm = T)})
by(mini_liteumk$fail,mini_liteumk$idno_original,function(x){sum(x)})

####checking if there are more than 1 fail in each####
sum(c(by(liteumk$`_d`,liteumk$idno_original,function(x){sum(x)}))>1, na.rm = T) #none

####cases that failed but didn't have VA####
#head(liteumk[liteumk$failure==1 & !is.na(liteumk$failure) & liteumk$hadva==0,])
dim(liteumk[which(liteumk$`_d`==1 & liteumk$hadva==0),])
#2572 died but didn't have VA 
tmp <- liteumk[which((liteumk$`_d`==1) & (liteumk$hadva==0)),]
#0 are duplicated
sum(duplicated(tmp$idno_original))

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

####checking age####
summary(liteumk$age)
summary(liteumk$`_t0`) #max:122
liteumk[liteumk$age==60,]
over90 <- liteumk[liteumk$`_t0`>90,]
centennials <- liteumk[liteumk$`_t0`>100,]
length(unique(centennials$idno_original))
#write.csv(centennials,'centennials.csv')

####cases with VA but overlapping records####
table(liteumk$`_d`,liteumk$hadva, exclude = NULL)
tmp <- liteumk[liteumk$hadva==1 & !is.na(liteumk$hadva) & is.na(liteumk$`_d`),]
tmp$idno_original %in% CODdataset[,1] #these records have no VA information either

####cases with entry>exit####
tmp <- liteumk[liteumk$entry>liteumk$exit,]
tmp <- liteumk[liteumk$entry==liteumk$exit & !is.na(liteumk$entry) & !is.na(liteumk$exit),]

####what does it mean by HIV status NA####
tmp <- liteumk[which(is.na(liteumk$hivstatus_broad)),]
sum(is.na(tmp$last_neg_date))

tmp2 <- liteumk[which(!is.na(liteumk$hivstatus_broad)),]
sum(is.na(tmp2$last_neg_date))

hivstatus <- table(liteumk$hivstatus_detail,liteumk$hivstatus_broad,  exclude=NULL)
hivstatus2 <- table(liteumk$allinfo_treat_pyramid,liteumk$hivstatus_broad,  exclude=NULL)
#write.csv(hivstatus,'hivstatus.csv')
#write.csv(hivstatus2,'hivstatus2.csv')
posneg <- liteumk[which(liteumk$hivstatus_detai=='Negative' & liteumk$hivstatus_broad=='Positive'),]

table(liteumk$allinfo_treat_pyramid,liteumk$hivstatus_broad,  exclude=NULL)

####Merging 2 datasets####
#liteumk & CODdataset
#look into 'sp_essential.R' file for this and onwards analysis
