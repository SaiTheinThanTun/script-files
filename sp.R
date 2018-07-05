#summer project

#library(foreign)
library(Rcpp)
library("readstata13")
library("Hmisc")

setwd("~/OneDrive/Summer Project/data/")

#Read stata data file
#umk <- read.dta13("alpha_uMkhanyakude-170601.dta")
#saveRDS(umk,"umk.RDS")
umk <- readRDS("umk.RDS")

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


###exploratory####
#InSilicoVA
isv <- read.csv("InSilicoVA outputs - uMkhanyakude.csv")
#names(isv)
isvID <- isv[,1]
isv <- isv[,-1]
isvCOD <- names(isv)
COD <- NA
for(i in 1:nrow(isv)){
  COD[i] <- isvCOD[which.max(isv[i,])]
}
CODdataset <- cbind(isvID,COD)

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
table(umk$hadva, umk$failure, exclude = NULL)
sum(umk$hadva, na.rm=T)
sum(umk$fail, na.rm=T)
#can't be checked as hadva has 1 for all episodes of each case 
#whereas fail occurs at most 1 per person
#test <- umk[umk$hadva==1 & is.na(umk$failure),]

####2 FAIL variables####
table(umk$fail, umk$failure, exclude=NULL)

####fail vs HIV####
table(umk$hivstatus_detail, umk$fail, exclude=NULL)

####no. of participants####
#114980
length(unique(umk$idno_original))

#no. of deaths: 15972
table(umk$`_d`, umk$fail, exclude=NULL)
#table(umk$`_d`, umk$failure, exclude=NULL) #failure var removed

#no. of deaths which had va (main dataset)
table(umk$fail, umk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

####checking cases that failed (failure==1)####
head(umk[umk$failure==1 & !is.na(umk$failure),])
case15 <- umk[umk$idno_original==15,]
case12 <- umk[umk$idno_original==12,]
case108 <- umk[umk$idno_original==108,]
case13 <- umk[umk$idno_original==13,]
case61051 <- umk[umk$idno_original==61051,]

####checking cases that didn't fail at all####
#sum of failure with each unique case
mini_umk <- umk[c(1:1000),]
by(mini_umk$failure,mini_umk$idno_original,function(x){sum(x,na.rm = T)})
by(mini_umk$fail,mini_umk$idno_original,function(x){sum(x)})

####checking if there are more than 1 fail in each####
sum(c(by(umk$fail,umk$idno_original,function(x){sum(x)}))>1, na.rm = T)

####cases that failed but didn't have VA####
head(umk[umk$failure==1 & !is.na(umk$failure) & umk$hadva==0,])

####idno vs idno_original#####
#they are unique in their own ways :)
length(unique(umk$idno))
length(unique(umk$idno_original))

####allinfo_treat_pyramid vs hivtreat####
table(umk$allinfo_treat_pyramid, umk$hivtreat)

####allinfo_treat_pyramid vs hivevertreat####
table(umk$allinfo_treat_pyramid, umk$hivevertreat)

####battle of 2 sexes####
table(umk$sex, umk$male, exclude=NULL)
