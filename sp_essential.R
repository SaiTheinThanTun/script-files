#summer project
#dataflow: sp.do -> sp.R -> sp_essential.R
#essential info to put in report
lfc <- 1 #left censoring at '2007-01-01'
rtc <- 1 #right censoring for age above 100
creation <- FALSE

library(Rcpp)
library("readstata13")
library("Hmisc")
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(ggplot2)
library(data.table)
library(survey) #for svycontrast
source('~/OneDrive/Summer Project/script-files/sp_functions.R')

setwd("~/OneDrive/Summer Project/data/")

#Read lite data file
liteumk <- readRDS("liteumk.RDS")

#InSilicoVA
CODdataset <- readRDS("CODdataset.RDS")

####consistency checks####
####no. of participants####
#114980 -> 114068
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
length(unique(liteumk[which(is.na(liteumk$`_d`)),]$idno_original)) #2034 cases only instead of 2037
liteumk <- liteumk[which(!is.na(liteumk$`_d`)),] #remove missing from `_d` [already described in report]
#additional cases with same entry and exit dates found out later [on date: 20180719] However, see below
#liteumk <- liteumk[which(liteumk$entry!=liteumk$exit),] #this is an artefact of spliting by birthday/positive date
#the discrepancy in # is because some cases have records with NA in `_d` (ie. same entry & exit date)

#cases after removal of missing outcomes
length(unique(liteumk$idno_original)) #114068

####Left censoring (change terms) for 2007 when HIV testing is done for all adults (15+)####
if(lfc){
  liteumk <- liteumk[which((liteumk$entry>='2007-01-01')),] #1271610 remains representing 88693 individuals, 951,631 records discarded
}

#cases after left censoring
length(unique(liteumk$idno_original)) #88693

####Right censoring (change terms) for records with age>100####
if(rtc){
  liteumk <- liteumk[which((liteumk$age<=100)),] #1271372 remains representing 88686 individuals, 238 records discarded
}

#cases after right censoring
length(unique(liteumk$idno_original)) #88686

####issue with 'residence'####
if(FALSE){
  mostFreqRes <- by(liteumk$residence,liteumk$idno_original,function(x){tableR(x)})
  freq.id <- table(liteumk$idno_original)
  
  residence2 <- rep(mostFreqRes, freq.id)
  liteumk <- cbind(liteumk, residence2)
}


#tmp5 <- liteumk[is.na(liteumk$residence),] #testing
# reveal(13)
# reveal(14)

#essential strats here
####Merging 2 datasets####
#liteumk & CODdataset
colnames(CODdataset)[1] <- "idno_original" #change to match names
CODdataset <- as.data.frame(CODdataset)
dim(CODdataset) #12705 x2 #10171 x2
dim(liteumk) #1271372 x43

dat <- merge(liteumk,CODdataset, by="idno_original", all.x = TRUE)

dim(dat) #1271372 x44
sum(!is.na(dat$COD)) #58222 <- #40185

#no. of deaths
dim(dat[dat$`_d`==1,]) #7566 deaths
#sum(c(by(dat$fail0, dat$idno_original, function(x){sum(x)>=1}))) #same as above

#sum(c(by(dat$fail0, dat$idno_original, function(x){sum(x)>1}))) #checking if death in a person is recorded twice, 0

#no. of deaths which had VA
dim(dat[dat$`_d`==1 & !is.na(dat$COD),]) #6222 #4530 #82%

####Adding new fail variable 'fail2'####
#loop, test 2 conditions (liteumk$fail==1 & liteumk$COD %in% accidentNames)
accidentNames <- readRDS('accidentNames.RDS')
dim(dat) #1271372      44
if(FALSE){
  #THIS TAKES TOO LONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  for(i in 1:nrow(dat)){
    if(dat$fail0[i]==1 & dat$COD[i] %in% accidentNames) dat$fail2[i] <- 1
    else dat$fail2[i] <- 0
  }
}
#another faster way of doing the above
fail2 <- (dat$`_d`==1)*(dat$COD %in% accidentNames)
dat <- cbind(dat,fail2)
dim(dat) #1271372      45
#sum(is.na(dat$fail2))
table(dat$fail2, exclude = NULL) #1270599non-case ?giving 773 deaths from accidents only ?817 deaths from prev result


#checking 2 methods of deriving fail2
#tmp1 <- dat[dat$fail2==1,]
#tmp2 <- dat[fail2fast==1,]
#identical(tmp1, tmp2) #TRUE

#missing of main exposure variable
#table(datL$seroconFixed,datL$fail2, exclude = NULL) #5 records with missing hivstatus, seroconFixed
#datL[datL$`_d` ==0 & datL$va==TRUE,] #12 records with no fail but had va, Undeads with COD

#checking hivstatus_broad
#already done in stata, 3rd version of dataset
#dat$neg_yr <- (dat$exit - dat$last_neg_date)/365.25
#dat$negAbove5 <- dat$neg_yr>5 #dat$neg_yr>=5
#labeling
dat$hivstale5y <- as.factor(dat$hivstale5y) 
levels(dat$hivstale5y) <- c("Negative", "Positive" ,"Unknown")
dat$missingFixed <- as.factor(dat$missingFixed) 
levels(dat$missingFixed) <- c("Negative", "Positive" ,"Unknown")
dat$seroconFixed <- as.factor(dat$seroconFixed) 
levels(dat$seroconFixed) <- c("Negative", "Positive" ,"Unknown")
dat$allFixed <- as.factor(dat$allFixed) 
levels(dat$allFixed) <- c("Negative", "Positive" ,"Unknown")

table(dat$negAbove5,dat$hivstatus_broad,exclude = NULL)
table(dat$negAbove5,dat$hivstale5y,exclude = NULL)
table(dat$negAbove5,dat$missingFixed,exclude = NULL)
table(dat$negAbove5,dat$seroconFixed,exclude = NULL)
table(dat$negAbove5,dat$allFixed,exclude = NULL)
#dat$negAbove52 <- dat$neg_yr>5 #dat$neg_yr>=5 #better way to come!!! edited in stata

#checking negative above 5 years after last hiv negative test
#tmp <- dat[which((dat$negAbove5==TRUE) & (dat$hivstale5y=='Negative')),]
#View(tmp)
#length(unique(tmp$idno_original))


#tmp2 <- dat[which(!is.na(dat$last_neg_date) & dat$neg_yr>0),]
#View(tmp2)
#checking HIV status
if(FALSE){
  rnames <- c('idno_original','last_neg_date','frst_pos_date','exit','entry','hivstale5y','hivstatus_detail','neg_yr','negAbove5')
  reveal2(57818, rnames) #borderline
  reveal2(144675, rnames) #case
  reveal2(106727, rnames) #change to unknown after 1 year
  reveal2(16, rnames) #change to unknown after 1 year
  reveal2(17, rnames) #change to unknow on the last_neg_date
  
  
  View(dat[which(dat$idno_original %in% c(240, 272, 358, 492, 1401, 21466, 21490, 153504, 153792, 154144)),]) #change in 2 years
  View(dat[which(dat$idno_original %in% c(56, 287, 288, 359, 739, 789, 1129)),])#change in 3 years
  View(dat[which(dat$idno_original %in% c(17, 84,  238, 278, 14922)),])#seroconverters
}


#further cleaning done here
#the aim is to describe only to describe the final (FINAL) dataset
#the last state of each individual on their last visit???

#dat <- dat[which(!is.na(dat$COD)),]

#subsetting the LAST records
dat <- dat[order(dat$entry, decreasing = T),] #LAST record
#dat <- dat[order(dat$exit, decreasing = T),]
datL <- dat[!duplicated(dat$idno_original),]

#where are missing VA information distributed?####
datL$nova <- 1*(is.na(datL$COD)&(datL$fail0==1))
table(datL$fail0, datL$nova, exclude = NULL) #12 had va but didn't die
#table(datL$sex, datL$nova)/rowSums(table(datL$sex, datL$nova))
dead <- datL[datL$fail0==1,]
#by sex
table(dead$sex, dead$nova)
table(dead$sex, dead$nova)/rowSums(table(dead$sex, dead$nova))
#by agegrp15
table(dead$agegrp15, dead$nova)
table(dead$agegrp15, dead$nova)/rowSums(table(dead$agegrp15, dead$nova))
#by allFixed
table(dead$allFixed, dead$nova)
table(dead$allFixed, dead$nova)/rowSums(table(dead$allFixed, dead$nova))
#by art_status2
table(dead$art_status2, dead$nova)
table(dead$art_status2, dead$nova)/rowSums(table(dead$art_status2, dead$nova))

#no. of external injury death by sex and HIV####
table(dat$fail2,dat$sex)
table(dat$fail2,dat$sex)/sum(dat$fail2)
table(dat$fail2,dat$allFixed)
table(dat$fail2,dat$allFixed)/sum(dat$fail2)
table(dat$fail2,dat$allFixed)/colSums(table(dat$fail2,dat$allFixed))

#% of external injury per calendar year####
table(dat$fail2,dat$ct)
table(dat$fail2,dat$ct)/(table(dat$fail0,dat$ct))[2,]

####Descriptive statistics of the last records####
varDes <- c('sex','residence','hadva','failure','agegrp',
            'hivstatus_detail','hivstale5y','allinfo_treat_pyramid',
            'hivtreat','hivevertreat','preg_at_art_start',
            'art_available','art_avail_cat', 'COD')
#discarded
#,'age','agelastseen',
des <- list()

for(i in 1:length(varDes)){
  des[[i]] <- table(datL[,which(colnames(datL)==varDes[i])], exclude = NULL)
}
names(des) <- varDes


#some checking on hadva and COD
datL$va <- !is.na(datL$COD)
table(datL$hadva, datL$va, exclude = NULL)
table(datL$`_d`, datL$va, exclude = NULL) #12 are not dead but have VA!!!

#get hivstale5y missing to unknown #recoding unknown
#dat$hivstale5y[is.na(dat$hivstale5y)] <- 'Unknown'
names(dat)[names(dat)=="_t0"] <- "time0"
names(dat)[names(dat)=="_t"] <- "timex"
names(dat)[names(dat)=="_d"] <- "fail0" 

#new variable for period #before dim(dat) is 1271372      45
dat$period.2011_15 <- dat$entry>="2011-01-01"
dat$agegrp15 <- cut(dat$age, c(15,30,45,60, 122), include.lowest = T, right = F)
with(dat,table(agegrp15, agegrp))
levels(dat$agegrp15) <- c("15-29","30-44", "45-59","60+")

#new variable for HIV categorization
#1 No more info, 2 HIV negative, 3 offART: Never treated and Interrupted ART, 4 onART: Early ART, Stable ART, 5 HIV+ No more info
dat$art_status <- (1*(dat$allFixed=='Unknown'))+(2*(dat$allFixed=='Negative'))+(3*(((dat$allFixed=='Positive')&(dat$hivtreat=='Never treated' | dat$hivtreat=='Interrupted ART'))|((dat$allFixed=='Positive')&(dat$hivtreat=='HIV negative'))))+(4*((dat$allFixed=='Positive')&(dat$hivtreat=='Early ART' | dat$hivtreat=='Stable ART')))+(5*((dat$allFixed=='Positive')&(dat$hivtreat=='No more info')))
#hiv positive on seroconversion 4777 is added as offART
if('Off ART' %in% dat$hivtreat) print('check recoding for Off ART!!!!') #flag anticipating for more data added which might have Off ART individuals
dat$art_status <- factor(dat$art_status) #, levels = c('No more info','HIV negative', 'offART', 'onART'))
levels(dat$art_status) <- c('HIV Unknown','HIV negative', 'HIV+ offART', 'HIV+ onART', 'HIV+ No more info')
table(dat$hivtreat,dat$art_status)

#1 Unknown, 2 HIV negative, 3 HIV+ never started ART, 4 HIV+ ever started ART, 5 HIV+ No more info
dat$art_status2 <- (1*(dat$allFixed=='Unknown'))+(2*(dat$allFixed=='Negative'))+(3*(((dat$allFixed=='Positive')&(dat$hivtreat=='Never treated'))|((dat$allFixed=='Positive')&(dat$hivtreat=='HIV negative'))|((dat$allFixed=='Positive')&(dat$hivtreat=='No more info'))))+(4*((dat$allFixed=='Positive')&(dat$hivtreat=='Early ART' | dat$hivtreat=='Stable ART'| dat$hivtreat=='Interrupted ART')))
dat$art_status2 <- factor(dat$art_status2) 
levels(dat$art_status2) <- c('HIV Unknown','HIV Negative', 'HIVposUnknownART', 'HIVposART')
table(dat$hivtreat,dat$art_status2)

#new variable for calander time
dat$ct <- year(dat$exit)

#new variable for residence, with missing as a category
dat$residence2 <- dat$residence
dat$residence2[which(is.na(dat$residence2))] <- 'missing'
dat$residence2 <- factor(dat$residence2)
levels(dat$residence2) <- c('Urban', 'Semi-urban', 'Rural', 'missing')
class(dat$residence2)
table(dat$residence,dat$residence2, exclude = NULL)

#summary(dat[dat$period.2011_15==TRUE,]$entry)
#summary(dat[dat$period.2011_15==FALSE,]$entry)

#fail definitions
#STATA defined death, 'failure' is the variable STATA used for this. there's also 'fail' which is not used
#fail2 is external injury related death

#saveRDS(dat, 'dat_lifetable.RDS') #processed data file to be used for lifetable
#library(foreign)
#write.dta(dat,'processed2.dta')

####sex, age, py distributions####
desTableNames <- c('Individuals','Person-years','All cause: deaths', 'Rate', 'Low 95%CI', 'High 95%CI','Ext. Inj: deaths', 'Rate', 'Low 95%CI', 'High 95%CI')
if(FALSE){
per <- 10000

no.individuals <- table(datL$agegrp,datL$sex,  exclude = NULL)
no.individuals.hiv <- table(datL$allFixed,  exclude = NULL)
if(creation) write.csv(no.individuals.hiv,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_numberPerHIV.csv",sep = ""))
datLf <- datL[datL$sex=="Women",]
datLm <- datL[datL$sex=="Men",]
# A <- <- table(datL$agegrp,datL$sex,  exclude = NULL)
# datLallcause <- datL[datL$`_d`==1,]
# datLextInj <- datL[datL$fail2==1,]
# B <- table(datLallcause$agegrp,datLallcause$sex,  exclude = NULL)
# Cee <- table(datLextInj$agegrp,datLextInj$sex,  exclude = NULL)
# agebrkdwn <- cbind(A,B,Cee)


#3 places to change .f , "Women", no.individuals[,2]
#pyaa: person-years by agegroup all causes
pyaa.f <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat[dat$sex=="Women",], scale = 1)
# pyaa.f.summ <- summary(pyaa.f, rate = T, ci.r = T)
pyaa.f.py <- pyaa.f$pyears
pyaa.f.event <- pyaa.f$event
pyaa.f.rate <- pyaa.f.event/pyaa.f.py
pyaa.f.lo <- exp(log(pyaa.f.rate)-(1.96/sqrt(pyaa.f.event)))
pyaa.f.hi <- exp(log(pyaa.f.rate)+(1.96/sqrt(pyaa.f.event)))

#pyae: person-years by agegroup, external injuries
pyae.f <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$sex=="Women",], scale = 1)
# pyae.f.summ <- summary(pyae.f, rate = T, ci.r = T)
pyae.f.py <- pyae.f$pyears
pyae.f.event <- pyae.f$event
pyae.f.rate <- pyae.f.event/pyae.f.py
pyae.f.lo <- exp(log(pyae.f.rate)-(1.96/sqrt(pyae.f.event)))
pyae.f.hi <- exp(log(pyae.f.rate)+(1.96/sqrt(pyae.f.event)))

f.dist <- cbind(no.individuals[,2],round(pyaa.f.py),pyaa.f.event,round(pyaa.f.rate*per,1),round(pyaa.f.lo*per,1),round(pyaa.f.hi*per,1),pyae.f.event,round(pyae.f.rate*per,1),round(pyae.f.lo*per,1),round(pyae.f.hi*per,1))

tot.indi <- sum(no.individuals[,2])
tot.pyaa.f.py <- sum(pyaa.f.py)
tot.pyaa.f.event <- sum(pyaa.f.event)
tot.pyaa.f.rate <- tot.pyaa.f.event/tot.pyaa.f.py
tot.pyaa.f.lo <- exp(log(tot.pyaa.f.rate)-(1.96/sqrt(tot.pyaa.f.event)))
tot.pyaa.f.hi <- exp(log(tot.pyaa.f.rate)+(1.96/sqrt(tot.pyaa.f.event)))
tot.pyae.f.py <- sum(pyae.f.py)
tot.pyae.f.event <- sum(pyae.f.event)
tot.pyae.f.rate <- tot.pyae.f.event/tot.pyae.f.py
tot.pyae.f.lo <- exp(log(tot.pyae.f.rate)-(1.96/sqrt(tot.pyae.f.event)))
tot.pyae.f.hi <- exp(log(tot.pyae.f.rate)+(1.96/sqrt(tot.pyae.f.event)))
Total <- c(tot.indi,round(tot.pyaa.f.py),tot.pyaa.f.event,round(tot.pyaa.f.rate*per,1),round(tot.pyaa.f.lo*per,1),round(tot.pyaa.f.hi*per,1),tot.pyae.f.event,round(tot.pyae.f.rate*per,1),round(tot.pyae.f.lo*per,1),round(tot.pyae.f.hi*per,1))

f.dist <- rbind(f.dist, Total)
colnames(f.dist) <- desTableNames
f.dist
#check result
py.f <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ sex, data=dat, scale = 1)
summary(py.f,rate = T, ci.r = T)
if(creation) write.csv(f.dist,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_Women.csv",sep = "") )

####Men####
#3 places to change .m , "Men", no.individuals[,1]
#pyaa: person-years by agegroup all causes
pyaa.m <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat[dat$sex=="Men",], scale = 1)
# pyaa.m.summ <- summary(pyaa.m, rate = T, ci.r = T)
pyaa.m.py <- pyaa.m$pyears
pyaa.m.event <- pyaa.m$event
pyaa.m.rate <- pyaa.m.event/pyaa.m.py
pyaa.m.lo <- exp(log(pyaa.m.rate)-(1.96/sqrt(pyaa.m.event)))
pyaa.m.hi <- exp(log(pyaa.m.rate)+(1.96/sqrt(pyaa.m.event)))

#pyae: person-years by agegroup, external injuries
pyae.m <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$sex=="Men",], scale = 1)
# pyae.m.summ <- summary(pyae.m, rate = T, ci.r = T)
pyae.m.py <- pyae.m$pyears
pyae.m.event <- pyae.m$event
pyae.m.rate <- pyae.m.event/pyae.m.py
pyae.m.lo <- exp(log(pyae.m.rate)-(1.96/sqrt(pyae.m.event)))
pyae.m.hi <- exp(log(pyae.m.rate)+(1.96/sqrt(pyae.m.event)))

m.dist <- cbind(no.individuals[,1],round(pyaa.m.py),pyaa.m.event,round(pyaa.m.rate*per,1),round(pyaa.m.lo*per,1),round(pyaa.m.hi*per,1),pyae.m.event,round(pyae.m.rate*per,1),round(pyae.m.lo*per,1),round(pyae.m.hi*per,1))

tot.indi <- sum(no.individuals[,1])
tot.pyaa.m.py <- sum(pyaa.m.py)
tot.pyaa.m.event <- sum(pyaa.m.event)
tot.pyaa.m.rate <- tot.pyaa.m.event/tot.pyaa.m.py
tot.pyaa.m.lo <- exp(log(tot.pyaa.m.rate)-(1.96/sqrt(tot.pyaa.m.event)))
tot.pyaa.m.hi <- exp(log(tot.pyaa.m.rate)+(1.96/sqrt(tot.pyaa.m.event)))
tot.pyae.m.py <- sum(pyae.m.py)
tot.pyae.m.event <- sum(pyae.m.event)
tot.pyae.m.rate <- tot.pyae.m.event/tot.pyae.m.py
tot.pyae.m.lo <- exp(log(tot.pyae.m.rate)-(1.96/sqrt(tot.pyae.m.event)))
tot.pyae.m.hi <- exp(log(tot.pyae.m.rate)+(1.96/sqrt(tot.pyae.m.event)))
Total <- c(tot.indi,round(tot.pyaa.m.py),tot.pyaa.m.event,round(tot.pyaa.m.rate*per,1),round(tot.pyaa.m.lo*per,1),round(tot.pyaa.m.hi*per,1),tot.pyae.m.event,round(tot.pyae.m.rate*per,1),round(tot.pyae.m.lo*per,1),round(tot.pyae.m.hi*per,1))

m.dist <- rbind(m.dist, Total)
colnames(m.dist) <- desTableNames
m.dist
#check result
py.m <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ sex, data=dat, scale = 1)
summary(py.m,rate = T, ci.r = T)
if(creation) write.csv(m.dist,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_Men.csv",sep = "") )


#descriptions on HIV status
#HIV_asmr, to calculate standardized rates####
hiv_Negative_rates <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Negative',], scale = 1)
hiv_Negative_rates <- pyears2(hiv_Negative_rates, per=10000)

hiv_Positive_rates <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Positive',], scale = 1)
hiv_Positive_rates <- pyears2(hiv_Positive_rates, per=10000)

hiv_Unknown_rates <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Unknown',], scale = 1)
hiv_Unknown_rates <- pyears2(hiv_Unknown_rates, per=10000)

hiv_asmr <- rbind(hiv_Negative_rates,hiv_Positive_rates,hiv_Unknown_rates)
if(creation) write.csv(hiv_asmr,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_hiv_asmr.csv",sep = "") )

#HIV_asmr, agegrp####
#5year agegroup
#HIV_asmr, Men####
hiv_Negative_rates_Men <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Negative'&dat$sex=='Men',], scale = 1)
hiv_Negative_rates_Men <- pyears2(hiv_Negative_rates_Men, per=10000)
hiv_Negative_rates_Men_Rate <- hiv_Negative_rates_Men$Rate

hiv_Positive_rates_Men <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Positive'&dat$sex=='Men',], scale = 1)
hiv_Positive_rates_Men <- pyears2(hiv_Positive_rates_Men, per=10000)
hiv_Positive_rates_Men_Rate <- hiv_Positive_rates_Men$Rate

hiv_Unknown_rates_Men <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Unknown'&dat$sex=='Men',], scale = 1)
hiv_Unknown_rates_Men <- pyears2(hiv_Unknown_rates_Men, per=10000)
hiv_Unknown_rates_Men_Rate <- hiv_Unknown_rates_Men$Rate

Men_Rate <- as.data.frame(cbind(hiv_Negative_rates_Men_Rate,hiv_Positive_rates_Men_Rate,hiv_Unknown_rates_Men_Rate))
Men_Rate <- as.data.frame(cbind(Men_Rate,row.names(hiv_Negative_rates_Men), rep("Men",nrow(Men_Rate))))
colnames(Men_Rate) <- c('Negative','Positive','Unknown','Age','Sex')

#HIV_asmr, Women####
hiv_Negative_rates_Women <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Negative'&dat$sex=='Women',], scale = 1)
hiv_Negative_rates_Women <- pyears2(hiv_Negative_rates_Women, per=10000)
hiv_Negative_rates_Women_Rate <- hiv_Negative_rates_Women$Rate

hiv_Positive_rates_Women <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Positive'&dat$sex=='Women',], scale = 1)
hiv_Positive_rates_Women <- pyears2(hiv_Positive_rates_Women, per=10000)
hiv_Positive_rates_Women_Rate <- hiv_Positive_rates_Women$Rate

hiv_Unknown_rates_Women <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$allFixed=='Unknown'&dat$sex=='Women',], scale = 1)
hiv_Unknown_rates_Women <- pyears2(hiv_Unknown_rates_Women, per=10000)
hiv_Unknown_rates_Women_Rate <- hiv_Unknown_rates_Women$Rate

Women_Rate <- as.data.frame(cbind(hiv_Negative_rates_Women_Rate,hiv_Positive_rates_Women_Rate,hiv_Unknown_rates_Women_Rate))
Women_Rate <- as.data.frame(cbind(Women_Rate,row.names(hiv_Negative_rates_Women), rep("Women",nrow(Women_Rate))))
colnames(Women_Rate) <- c('Negative','Positive','Unknown','Age','Sex')

#plotting ASMR from external injury by sex and HIV
ASMR_plot <- rbind(Men_Rate,Women_Rate)
ASMR_plot_m <- melt(ASMR_plot, id.vars = c('Age','Sex'))
colnames(ASMR_plot_m) <- c('Age','Sex','HIV_status','Rate')

png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asmr_HIV_sex.png",sep = ""), width = 1100, height = 800)
ggplot(ASMR_plot_m)+
  geom_line(aes(x=Age,y=Rate, col=HIV_status, group=HIV_status))+
  facet_wrap(~Sex)+
  ylab("Age-specific external injury mortality rate per 10,000 PY")+
  xlab("Age group")+
  theme(text = element_text(size=16), legend.position = "bottom")
dev.off()

#HIV_asmr, single agegroup doesn't work

#neg for 5 years + all fixed, final description####
#inj cause
allFixed.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat, scale = 1)
allFixed.rate <- pyears2(allFixed.rate, per = 10000)
if(creation) write.csv(round(allFixed.rate,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_rate.csv",sep = "") )
#all cause
allFixed.allcause.rate <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data=dat, scale = 1)
allFixed.allcause.rate <- pyears2(allFixed.allcause.rate, per = 10000)
if(creation) write.csv(round(allFixed.allcause.rate,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_allcause_rate.csv",sep = "") )
#rates among HIV status, by sex####
#subsetting the LAST records
dat <- dat[order(dat$entry, decreasing = T),] #LAST record
#dat <- dat[order(dat$exit, decreasing = T),]
datL <- dat[!duplicated(dat$idno_original),]
individuals.hiv.sex <- table(datL$allFixed, datL$sex,  exclude = NULL)
#inj cause, Men
allFixed.rate.Men <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Men',], scale = 1)
allFixed.rate.Men <- pyears2(allFixed.rate.Men, per = 10000)
#if(creation) write.csv(round(allFixed.rate.Men,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_rate.csv",sep = "") )
#all cause, Men
allFixed.allcause.rate.Men <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data=dat[dat$sex=='Men',], scale = 1)
allFixed.allcause.rate.Men <- pyears2(allFixed.allcause.rate.Men, per = 10000)
#if(creation) write.csv(round(allFixed.allcause.rate.Men,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_allcause_rate.csv",sep = "") )
HIVtable_Men <- cbind(individuals.hiv.sex[,1],allFixed.allcause.rate.Men,allFixed.rate.Men)
if(creation) write.csv(round(HIVtable_Men,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIVtable_Men.csv",sep = "") )
#inj cause, Women
allFixed.rate.Women <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Women',], scale = 1)
allFixed.rate.Women <- pyears2(allFixed.rate.Women, per = 10000)
#if(creation) write.csv(round(allFixed.rate.Women,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_rate.csv",sep = "") )
#all cause, Women
allFixed.allcause.rate.Women <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data=dat[dat$sex=='Women',], scale = 1)
allFixed.allcause.rate.Women <- pyears2(allFixed.allcause.rate.Women, per = 10000)
#if(creation) write.csv(round(allFixed.allcause.rate.Women,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_allcause_rate.csv",sep = "") )
HIVtable_Women <- cbind(individuals.hiv.sex[,2],allFixed.allcause.rate.Women,allFixed.rate.Women)
if(creation) write.csv(round(HIVtable_Women,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIVtable_Women.csv",sep = "") )


#inj cause by art_status2
art_status2.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ art_status2, data=dat[dat$sex=='Men',], scale = 1)
art_status2.rate <- pyears2(art_status2.rate, per = 10000)
if(creation) write.csv(round(art_status2.rate,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_art_status2_rate_Men.csv",sep = "") )
art_status2.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ art_status2, data=dat[dat$sex=='Women',], scale = 1)
art_status2.rate <- pyears2(art_status2.rate, per = 10000)
if(creation) write.csv(round(art_status2.rate,1),paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_art_status2_rate_Women.csv",sep = "") )

#overall (the whole population) rates####
#all cause
allcause.rate <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ 1, data=dat, scale = 1)
allcause.rate <- pyears2(allcause.rate, per = 10000)
if(creation) write.csv(allcause.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allcause_rate.csv",sep = "") )
#injury cause
inj.cause.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ 1, data=dat, scale = 1)
inj.cause.rate <- pyears2(inj.cause.rate, per = 10000)
if(creation) write.csv(inj.cause.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_cause_rate.csv",sep = "") )
#injury cause_overall
inj.cause.rate2 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat, scale = 1)
inj.cause.rate2 <- pyears2(inj.cause.rate2, per = 10000)
if(creation) write.csv(inj.cause.rate2,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_cause_overall_asmr.csv",sep = "") )
#injury cause per calendar year
inj.cause.rate3 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ ct, data=dat, scale = 1)
inj.cause.rate3 <- pyears2(inj.cause.rate3, per = 10000)
if(creation) write.csv(inj.cause.rate3,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_cause_calendaryear.csv",sep = "") )
#age-standardized injury mortality per calendar year####
WSP <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/WSP.csv")
#overall
std_inj_ct <- NA
for(i in 1:length(names(table(dat$ct)))){
    tmp <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat[dat$ct==names(table(dat$ct))[i],], scale = 1)
    tmp <- pyears2(tmp, per=10000)
    std_inj_ct[[i]] <- sum(tmp$Rate*WSP)
}
#Men
std_inj_ct_Men <- NA
dat_Men <- dat[dat$sex=='Men',]
for(i in 1:length(names(table(dat$ct)))){
  tmp <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat_Men[dat_Men$ct==names(table(dat_Men$ct))[i],], scale = 1)
  tmp <- pyears2(tmp, per=10000)
  std_inj_ct_Men[[i]] <- sum(tmp$Rate*WSP)
}
#Women
std_inj_ct_Women <- NA
dat_Women <- dat[dat$sex=='Women',]
for(i in 1:length(names(table(dat$ct)))){
  tmp <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat_Women[dat_Women$ct==names(table(dat_Women$ct))[i],], scale = 1)
  tmp <- pyears2(tmp, per=10000)
  std_inj_ct_Women[[i]] <- sum(tmp$Rate*WSP)
}
uMkhanyakude <- as.data.frame(rbind(std_inj_ct,std_inj_ct_Men,std_inj_ct_Women))
uMkhanyakude <- cbind(cbind(c('Both', 'Male', 'Female'),rep('uMkhanyakude',3)), uMkhanyakude)
colnames(uMkhanyakude) <- c('sex_name','location_name', names(table(dat$ct)))
uMkhanyakudeMelt <- reshape::melt(uMkhanyakude, id.vars=c('sex_name','location_name'))
colnames(uMkhanyakudeMelt) <- c('sex_name','location_name','year','val')
if(creation) write.csv(uMkhanyakudeMelt,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_std_rate_peryear.csv",sep = ""))
#no. of injury deaths against period for each sex####
#Men
dat.Men.inj <- dat[dat$sex=='Men' & dat$fail2==1,] #618 58
dat.Men.inj <- droplevels(dat.Men.inj)
dat.Men.inj.Tab <- table(dat.Men.inj$COD, dat.Men.inj$period.2011_15)
colnames(dat.Men.inj.Tab) <- c("2007-2010", "2011-2015")

#Women
dat.Women.inj <- dat[dat$sex=='Women' & dat$fail2==1,] #155 58
dat.Women.inj <- droplevels(dat.Women.inj)
dat.Women.inj.Tab <- table(dat.Women.inj$COD, dat.Women.inj$period.2011_15)
colnames(dat.Women.inj.Tab) <- c("2007-2010", "2011-2015")
dat.inj.tab <- cbind(dat.Women.inj.Tab, dat.Men.inj.Tab)
if(creation) write.csv(dat.inj.tab,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_x_period.csv",sep = ""))

#no. of injury deaths against hiv status, allFixed####
dat.inj <- dat[dat$fail2==1,]
dat.inj <- droplevels(dat.inj)
dat.inj$period.2011_15 <- as.factor(dat.inj$period.2011_15)
levels(dat.inj$period.2011_15) <- c("2007-10", "2011-15")
inj_x_hiv <- table(dat.inj$COD, dat.inj$allFixed, dat.inj$period.2011_15)
inj_x_hiv.m <- melt(inj_x_hiv, c("COD","allFixed","period.2011_15"))

# Stacked Percent
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_x_hiv_period.png",sep = ""), width = 800, height = 900)
ggplot(inj_x_hiv.m, aes(fill=COD, y=value, x=allFixed)) +
  facet_grid(. ~ period.2011_15)+
  geom_bar( stat="identity", position="fill") +
  xlab("HIV status")+ylab("Proportion among injury-related deaths in each group")+
  scale_fill_discrete(name = "Deaths from external injury")+
  theme(text = element_text(size=16))
dev.off()

#no. of injury deaths against hiv status by sex, allFixed####
dat.inj <- dat[dat$fail2==1,]
dat.inj <- droplevels(dat.inj)
dat.inj$period.2011_15 <- as.factor(dat.inj$period.2011_15)
levels(dat.inj$period.2011_15) <- c("2007-10", "2011-15")
inj_x_hiv <- table(dat.inj$COD, dat.inj$allFixed, dat.inj$period.2011_15, dat.inj$sex)
inj_x_hiv.m <- melt(inj_x_hiv, c("COD","allFixed","period.2011_15", "sex"))

#inj_x_hiv.abs <- dcast(inj_x_hiv.m, COD ~ allFixed+ period.2011_15+ sex, sum) #dataframe instead of graph
inj_x_hiv.abs <- dcast(inj_x_hiv.m, COD ~ allFixed+ period.2011_15+sex, sum)

lab <- colSums(inj_x_hiv.abs[,-1])
allFixed <- c(rep('Negative',4),rep('Positive',4), rep('Unknown',4))
period.2011_15 <- rep(c(rep('2007-10',2),rep('2011-15',2)),3)
sex <- rep(c('Men','Women'),6)
datlab.hiv.period.sex <- data.frame(x=allFixed,y=1.02, lab, sex,period.2011_15)

# Stacked Percent
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_x_hiv_period_sex.png",sep = ""), width = 1300, height = 900)
ggplot(inj_x_hiv.m) +
  facet_grid(. ~ sex+period.2011_15)+
  geom_bar(aes(fill=COD, y=value, x=allFixed), stat="identity", position="fill") +
  xlab("HIV status")+ylab("Proportion among injury-related deaths in each group")+
  scale_fill_discrete(name = "Deaths from external injury")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=datlab.hiv.period.sex, size = 6)
dev.off()

#facet_grid(. ~ period.2011_15+sex)+

#no. of injury deaths against hiv/art status, art_status2####
dat.inj <- dat[dat$fail2==1,]
dat.inj <- droplevels(dat.inj)
dat.inj$period.2011_15 <- as.factor(dat.inj$period.2011_15)
levels(dat.inj$period.2011_15) <- c("2007-10", "2011-15")
inj_x_hiv <- table(dat.inj$COD, dat.inj$art_status2, dat.inj$period.2011_15)
inj_x_hiv.m <- melt(inj_x_hiv, c("COD","art_status2","period.2011_15"))

# Stacked Percent
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_inj_x_hivart_period.png",sep = ""), width = 1100, height = 900)
ggplot(inj_x_hiv.m, aes(fill=COD, y=value, x=art_status2)) +
  facet_grid(. ~ period.2011_15)+
  geom_bar( stat="identity", position="fill") +
  xlab("HIV/ART status")+ylab("Proportion among injury-related deaths in each group")+
  scale_fill_discrete(name = "Deaths from external injury")+
  theme(text = element_text(size=14))
dev.off()
}

#py_agegrp_extInj <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat, scale = 1)
#pyae <- summary(py_agegrp_extInj, rate = T, ci.r = T)

#write.csv(agebrkdwn,)

#write.csv(varDes,'descript_names.csv')
#lapply(des, function(x){write.table(data.frame(x),'descriptives.csv',append = T,sep = ',')})


####Survival Analysis####

testOClist <- c("hivstatus_broad","hivstale5y", "missingFixed", "seroconFixed", "allFixed")
#testOC <- testOClist[1]

###SURVIVAL CURVE####
#by hivstatus exposure
for(i in 1:length(testOClist)){
  testOC <- testOClist[i]
  
  png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_",testOC,"_allCause.png",sep = ""), width = 800, height = 600)
  sv <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ dat[,which(names(dat) %in% testOC)], data=dat)
  #sv <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y, data=dat) #ALL CAUSES OF DEATH
  plot(sv, conf.int=FALSE, mark.time=FALSE, col=1:3, main=paste(testOC, ", all causes", sep = "")) #by HIV status
  legend('topright',legend = c('neg','pos','unk'), col=1:3, lty=1)
  dev.off()
  
  png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_",testOC,"_injuryCause.png",sep = ""), width = 800, height = 600)
  sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ dat[,which(names(dat) %in% testOC)], data=dat) #Deaths due to injuries alone: fail2
  #summary(sv2)
  #summary(sv2, times = c(20, 30, 40))
  plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main=paste(testOC, ", external injury", sep = "")) #by HIV status
  legend('bottomleft',legend = c('neg','pos','unk'), col=1:3, lty=1)
  #ggsurvplot(sv2) #takes too much time and memory
  dev.off()
}

#by hivstatus exposure, stratified by sex, #reported
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_injuryCause_Women.png",sep = ""), width = 800, height = 600)
sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Women',]) #Deaths due to injuries alone: fail2
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main=paste("allFixed", ", external injury, Women", sep = "")) #by HIV status
legend('bottomleft',legend = c('neg','pos','unk'), col=1:3, lty=1)
dev.off()

png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_injuryCause_Men.png",sep = ""), width = 800, height = 600)
sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Men',]) #Deaths due to injuries alone: fail2
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main=paste("allFixed", ", external injury, Men", sep = "")) #by HIV status
legend('bottomleft',legend = c('neg','pos','unk'), col=1:3, lty=1)
dev.off()

#4x survival curve, #reported####

png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_survivalcurve_report.png",sep = ""), width = 1000, height = 1100)
par(mfrow=c(2,2))
sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data=dat[dat$sex=='Women',], type="kaplan-meier") #Deaths due to all
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main="Survival from any cause of death, \nWomen", xlab = "Age", ylab = "Estimates of survival function") #by HIV status
legend(title='HIV status','bottomleft',legend = c('Negative','Positive','Unknown'), col=1:3, lty=1)

sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data=dat[dat$sex=='Men',], type="kaplan-meier") #Deaths due to all
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main="Survival from any cause of death, \nMen", xlab = "Age", ylab = "Estimates of survival function") #by HIV status
legend(title='HIV status','bottomleft',legend = c('Negative','Positive','Unknown'), col=1:3, lty=1)

sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Women',], type="kaplan-meier") #Deaths due to injuries alone: fail2
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .75, main="Survival from deaths related to \nexternal injury, Women", xlab = "Age", ylab = "Estimates of survival function") #by HIV status
legend(title='HIV status','bottomleft',legend = c('Negative','Positive','Unknown'), col=1:3, lty=1)

sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat[dat$sex=='Men',], type="kaplan-meier") #Deaths due to injuries alone: fail2
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .75, main="Survival from deaths related to \nexternal injury, Men", xlab = "Age", ylab = "Estimates of survival function") #by HIV status
legend(title='HIV status','bottomleft',legend = c('Negative','Positive','Unknown'), col=1:3, lty=1)
dev.off()
par(mfrow=c(1,1))

#survival curve for residence, notreported####
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_residence.png",sep = ""), width = 800, height = 600)
sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ residence, data=dat) 
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main=paste("survival from injury death by residence", sep = "")) #by HIV status
legend('bottomleft',legend = c('urban','semiurban','rural'), col=1:3, lty=1)
dev.off()


png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_residence2.png",sep = ""), width = 800, height = 600)
sv2 <- survfit(Surv(time=time0, time2 = timex, event = fail2) ~ residence2, data=dat) 
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8, main=paste("survival from injury death by residence", sep = "")) #by HIV status
legend('bottomleft',legend = c('urban','semiurban','rural', 'missing'), col=1:3, lty=1)
dev.off()


###COX REGRESSION####
#!!need to check proportionality of hazards
#basic
#svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ factor(hivstale5y), data=dat) 
#summary(svcox)
testOClist <- c("hivstatus_broad","hivstale5y", "missingFixed", "seroconFixed", "allFixed")
testOC <- testOClist[5] #c("hivstatus_broad","hivstale5y", "missingFixed", "seroconFixed", "allFixed")
svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)]), data=dat) 
summary(svcox) #Positive     1.022     0.9788    0.8151     1.281
base.ph <- cox.zph(svcox)
base.ph
plot(base.ph) #TEST PROPORTIONALITY
#cumulative hazard
plot(survfit(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data = dat),fun='cloglog',xlab = 'Years', ylab="Cumulative Hazard (log)")
plot(survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data = dat),fun='cloglog',xlab = 'Years', ylab="Cumulative Hazard (log)",lty=1:3)

po_con <- c('sex','age','residence')
#? timeposneg??? > Negative test expiry date in years
#? hivevertreat???
# NA assumed to be Unknown

#basic+res
#res.svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y+residence, data=dat) 
#summary(res.svcox)
res.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(residence), data=dat) 
summary(res.svcox) #Positive     1.112     0.8996    0.8606     1.436
res.ph <- cox.zph(res.svcox)
res.ph
plot(res.ph) #TEST PROPORTIONALITY

#basic+sex *
#sex.svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y+sex, data=dat) 
#summary(sex.svcox)
sex.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(sex), data=dat) 
summary(sex.svcox) #Positive    1.3333      0.750    1.0621    1.6737
sex.ph <- cox.zph(sex.svcox)
sex.ph
plot(sex.ph) #TEST PROPORTIONALITY

#basic+age 
#age.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(hivstale5y)+factor(agegrp), data=dat) #error
age.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age, data=dat)
summary(age.svcox) #Positive     1.048     0.9546    0.8353     1.314
age.ph <- cox.zph(age.svcox)
age.ph
plot(age.ph) #TEST PROPORTIONALITY

#basic+agegrp15
age15.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(agegrp15), data=dat)
summary(age15.svcox) #Positive    1.0615     0.9420    0.8398     1.342

#basic+period
period.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(period.2011_15), data=dat)
summary(period.svcox) #Positive  0.04195   1.04285  0.11542  0.364  0.71623

#basic+age+sex *
as.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(sex), data=dat) 
summary(as.svcox) #Positive    1.3070     0.7651    1.0419    1.6396
as.ph <- cox.zph(as.svcox)
as.ph
plot(as.ph) #TEST PROPORTIONALITY

#basic+age+res+sex **
ars.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(residence)+factor(sex), data=dat) 
summary(ars.svcox) #Positive    1.4302     0.6992    1.1058    1.8499
ars.ph <- cox.zph(ars.svcox)
ars.ph
plot(ars.ph) #TEST PROPORTIONALITY

#basic+agegrp+res+sex *
a15rs.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(agegrp15)+factor(residence)+factor(sex), data=dat) 
summary(a15rs.svcox) #Positive    1.4004     0.7141    1.0689     1.835

#unadjusted, report
svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ allFixed, data=dat) 
summary(svcox)
#basic+agegrp+res+sex+period ** report
a15rsp.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(agegrp15)+factor(residence)+factor(sex)+factor(period.2011_15), data=dat) 
summary(a15rsp.svcox) #Positive    1.4369     0.6959    1.0967    1.8826

#basic+agegrp+res+sex+period+interaction ** period strata specific 
a15rspi.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(agegrp15)+factor(residence)+factor(sex)+factor(period.2011_15)+allFixed:factor(period.2011_15), data=dat) 
summary(a15rspi.svcox) #Positive   1.6619     0.6017    1.1501     2.401

#stratam specific hazard ratios####
unk2007 <- svycontrast(a15rspi.svcox,c("allFixedUnknown"=1))
unk2011 <- svycontrast(a15rspi.svcox,c("allFixedUnknown"=1, "allFixedUnknown:factor(period.2011_15)TRUE"=1))
unk2007ci <- c(exp(unk2007[1]),exp(unk2007-(1.96*.1878))[1],exp(unk2007+(1.96*.1878))[1], NA) 
unk2011ci <- c(exp(unk2011[1]),exp(unk2011-(1.96*.1878))[1],exp(unk2011+(1.96*.1878))[1], NA) 
stra_spec_unk <- rbind(unk2007ci,unk2011ci)
colnames(stra_spec_unk) <- c('HR','lower','upper','p-value')
write.csv(stra_spec_unk,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_periodStratum_spec_Unk.csv",sep = ""))

pos2007 <- svycontrast(a15rspi.svcox,c("allFixedPositive"=1))
pos2011 <- svycontrast(a15rspi.svcox,c("allFixedPositive"=1, "allFixedPositive:factor(period.2011_15)TRUE"=1))

#pvalue for interaction: 0.272413
pos2007ci <- c(exp(pos2007[1]),exp(pos2007-(1.96*.1878))[1],exp(pos2007+(1.96*.1878))[1], 0.272413) #1.661889 1.150121 2.401379 
pos2011ci <- c(exp(pos2011[1]),exp(pos2011-(1.96*.1878))[1],exp(pos2011+(1.96*.1878))[1], NA) #1.2495098 0.8647311 1.8055031

stra_spec <- rbind(pos2007ci,pos2011ci)
colnames(stra_spec) <- c('HR','lower','upper','p-value')
write.csv(stra_spec,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_periodStratum_spec.csv",sep = ""))


#basic+agegrp+res+sex+period+interaction
a15rsip.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(agegrp15)+factor(residence)+factor(sex)+factor(period.2011_15)+allFixed:factor(sex), data=dat) 
summary(a15rsip.svcox) #Positive  1.5647     0.6391    1.1299    2.1668

a15irsp.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(agegrp15)+factor(residence)+factor(sex)+factor(period.2011_15)+allFixed:factor(agegrp15), data=dat) 
summary(a15irsp.svcox) #Positive     1.8999     0.5263    1.1869    3.0414

a15risp.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(agegrp15)+factor(residence)+factor(sex)+factor(period.2011_15)+allFixed:factor(residence), data=dat) 
summary(a15risp.svcox) #Positive   3.6109     0.2769   0.99050   13.1633

#anova
anova(a15rsp.svcox, a15rspi.svcox) #p .5109
anova(a15rspi.svcox, a15rs.svcox)

#producing RR (crude, partial, fully adj) and p values
results <- list()
for(i in 1:length(testOClist)){
  testOC <- testOClist[i]
  
  #crude
  svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)]), data=dat) 
  
  cph_summ <- summary(svcox) #summary(coxph object here)
  co_pos <- round(cph_summ$conf.int[1,1],4)
  lo_pos <- round(cph_summ$conf.int[1,3],4)
  hi_pos <- round(cph_summ$conf.int[1,4],4)
  co_un <- round(cph_summ$conf.int[2,1],4)
  lo_un <- round(cph_summ$conf.int[2,3],4)
  hi_un <- round(cph_summ$conf.int[2,4],4)
  waldp <- round(cph_summ$waldtest[3],6)
  rr_crude <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
  
  #partial
  #basic+age+sex
  as.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(sex), data=dat) 
  
  cph_summ <- summary(as.svcox) #summary(coxph object here)
  co_pos <- round(cph_summ$conf.int[1,1],4)
  lo_pos <- round(cph_summ$conf.int[1,3],4)
  hi_pos <- round(cph_summ$conf.int[1,4],4)
  co_un <- round(cph_summ$conf.int[2,1],4)
  lo_un <- round(cph_summ$conf.int[2,3],4)
  hi_un <- round(cph_summ$conf.int[2,4],4)
  waldp <- round(cph_summ$waldtest[3],6)
  rr_partial <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
  
  #fully adj
  #basic+age+res+sex
  ars.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(residence)+factor(sex), data=dat) 
  
  cph_summ <- summary(ars.svcox) #summary(coxph object here)
  co_pos <- round(cph_summ$conf.int[1,1],4)
  lo_pos <- round(cph_summ$conf.int[1,3],4)
  hi_pos <- round(cph_summ$conf.int[1,4],4)
  co_un <- round(cph_summ$conf.int[2,1],4)
  lo_un <- round(cph_summ$conf.int[2,3],4)
  hi_un <- round(cph_summ$conf.int[2,4],4)
  waldp <- round(cph_summ$waldtest[3],6)
  rr_full <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
  
  results[[i]] <- c(rr_crude, rr_partial, rr_full)
}
rr_results <- matrix(unlist(results),length(testOClist),3*3, byrow=T)
write.csv(rr_results,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_rateratios.csv",sep = ""))

#EXTRACTING COEFFICIENTS
if(FALSE){
  results <- list() #base structure, sample
  for(i in 1:length(testOClist)){
    testOC <- testOClist[i]
    
    cph_summ <- tmp #summary(coxph object here)
    #extracting exp(coef) and 95% CI # example
    co_pos <- round(cph_summ$conf.int[1,1],4)
    lo_pos <- round(cph_summ$conf.int[1,3],4)
    hi_pos <- round(cph_summ$conf.int[1,4],4)
    co_un <- round(cph_summ$conf.int[2,1],4)
    lo_un <- round(cph_summ$conf.int[2,3],4)
    hi_un <- round(cph_summ$conf.int[2,4],4)
    waldp <- round(cph_summ$waldtest[3],6)
    rr_crude <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
    rr_partial <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
    rr_full <- c(paste(co_pos," (", lo_pos," to ", hi_pos,")", sep = ""), paste(co_un," (", lo_un," to ", hi_un,")", sep = ""), waldp)
    results[[i]] <- c(rr_crude, rr_partial, rr_full)
  }
  final_results <- matrix(unlist(results),length(testOClist),3*3, byrow=T)
  
}

#cox regression revisited####
ars.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ allFixed+age+factor(residence)+factor(sex), data=dat) 
summary(ars.svcox)
ars.ph <- cox.zph(ars.svcox)
ars.ph
plot(ars.ph)

arsp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ allFixed+age+factor(residence)+factor(sex)+factor(period.2011_15), data=dat) 
summary(arsp.cox)
arsp.cox.ph <- cox.zph(arsp.cox)
arsp.cox.ph

#conceptual framework####
#univariate analyses, reported also in crude hazard ratios
#H0: Injury related mortality increases in male
summary(coxph(Surv(time=timex-time0, event = fail2) ~ factor(sex), data=dat)) #factor(sex)Women -1.52823   0.21692  0.08983 -17.01   <2e-16 ***

#H0: Injury related mortality increases in older age
summary(coxph(Surv(time=timex-time0, event = fail2) ~ factor(agegrp15), data=dat)) #factor(agegrp15)60+    0.49491   1.64035  0.11403  4.340 1.42e-05 ***

#H0: Injury related mortality is lower in 2011-15 because of improved care and living standard
summary(coxph(Surv(time=timex-time0, event = fail2) ~ factor(period.2011_15), data=dat)) #factor(period.2011_15)TRUE -0.24661   0.78145  0.07196 -3.427 0.000611 ***

#H0: Injury related mortality decreases in people living in rural area
summary(coxph(Surv(time=timex-time0, event = fail2) ~ factor(residence), data=dat)) #factor(residence)2 0.02143   1.02166  0.17342 0.124    0.902

#H0: HIV prevalence differs between urban, semi-urban and rural area
table(datL$allFixed, datL$residence)
prop.table(table(datL$allFixed, datL$residence),2)
chisq.test(datL$allFixed,datL$residence) #X-squared = 1230.4, df = 4, p-value < 2.2e-16
#chisq.test(table(datL$allFixed, datL$residence)) #same as above

#H0: HIV prevalence differs between sex
table(datL$allFixed, datL$sex)
prop.table(table(datL$allFixed, datL$sex),2)
chisq.test(datL$allFixed,datL$sex) #X-squared = 3443.4, df = 2, p-value < 2.2e-16

#H0: HIV prevalence differs by age
table(datL$allFixed, datL$agegrp15)
prop.table(table(datL$allFixed, datL$agegrp15),2)
chisq.test(datL$allFixed,datL$agegrp15) #X-squared = 10038, df = 6, p-value < 2.2e-16

#H0: HIV prevalence differs by follow up period
table(datL$allFixed, datL$period.2011_15)
prop.table(table(datL$allFixed, datL$period.2011_15),2)
chisq.test(datL$allFixed,datL$period.2011_15) #X-squared = 727.4, df = 2, p-value < 2.2e-16

#cox regression, reported####
crude.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed), data=dat) 
summary(crude.cox)

arsp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(residence), data=dat) 
summary(arsp.cox)
arsp.cox.ph <- cox.zph(arsp.cox)
arsp.cox.ph

arspi.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(residence)+factor(allFixed):factor(period.2011_15), data=dat) 
anova(arsp.cox, arspi.cox)


#cox regression, testing out residence2, not going well####
# crude.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ allFixed, data=dat) 
# summary(crude.cox)
a15rsp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(residence2), data=dat) 
summary(a15rsp.cox)
a15rsp.cox.ph <- cox.zph(a15rsp.cox)
a15rsp.cox.ph
#interaction between HIV status and period
a15rspi.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(residence2)+factor(allFixed):factor(period.2011_15), data=dat) 
summary(a15rspi.cox)
#anova, using residence2
anova(a15rsp.cox, a15rspi.cox) #p .1809
anova(a15sp.cox,a15rsp.cox) #p 0.1276, adding residence to the model makes no difference


#without residence in the model
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat)
a15sp.cox <- coxph(Surv(time=exit-entry, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat) 
a15sp.cox <- coxph(Surv(time2=as.numeric(exit), time=as.numeric(entry), event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox) #, transform = 'identity')
a15sp.cox.ph
survminer::ggcoxzph(a15sp.cox.ph)

for(i in 1:(nrow(a15sp.cox.ph$table)-1)){
  png(paste("~/OneDrive/Summer Project/output/fit/",gsub("\\:","",Sys.time()),"_",rownames(a15sp.cox.ph$table)[i],"_fit.png",sep = ""), width = 1100, height = 800)
  plot(a15sp.cox.ph[i], main=(paste(rownames(a15sp.cox.ph$table)[i])))
  abline(0,0, col=2)
  abline(h= a15sp.cox$coef[i], col=3, lwd=2, lty=2)
  
  dev.off()
}

#interaction, without residence in the model
a15spi.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(allFixed):factor(period.2011_15), data=dat) 
summary(a15spi.cox)
anova(a15sp.cox,a15spi.cox) #p 0.1812

anova(asp.cox,arsp.cox) #arsp is with missing value, thus it will not run!

#without residence in the model+calendar time
a15spt.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+ct, data=dat) 
summary(a15spt.cox)
a15spt.cox.ph <- cox.zph(a15spt.cox)
a15spt.cox.ph
#interaction, without residence in the model+calendar time
a15spti.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(allFixed):factor(period.2011_15), data=dat) 
summary(a15spti.cox)
anova(a15spt.cox,a15spti.cox)

#testing out new art_status category####
#without residence and with new ART categorization
dat$art_status2 <- relevel(dat$art_status2, ref = 'HIV negative')
a15spART2.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(art_status2)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat) 
summary(a15spART2.cox)
a15spART2.cox.ph <- cox.zph(a15spART2.cox)
a15spART2.cox.ph
a15spART2i.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(art_status2)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(art_status2):factor(period.2011_15), data=dat) 
summary(a15spART2i.cox)
anova(a15spART2.cox, a15spART2i.cox)

#without residence and with new ART categorization+calendar time
dat$art_status2 <- relevel(dat$art_status2, ref = 'HIV negative')
a15spART2t.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(art_status2)+factor(sex)+factor(agegrp15)+ct, data=dat) 
summary(a15spART2t.cox)
a15spART2t.cox.ph <- cox.zph(a15spART2t.cox)
a15spART2t.cox.ph
#no period, no interaction to test
a15spART2ti.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(art_status2)+factor(sex)+factor(agegrp15)+ct+factor(art_status2):factor(period.2011_15), data=dat) 
summary(a15spART2ti.cox)
anova(a15spART2t.cox, a15spART2ti.cox)

#lexis
if(FALSE){
  lx <- Lexis(entry = list(per=cal.yr(entry)), exit = list(per=cal.yr(exit)+.000001), id=idno_original, exit.status=fail0, data=dat)
  #saveRDS(lx,'lexisData20180719.RDS')
  #lx <- readRDS('lexisData20180719.RDS')
  #+.000001 was added in exit year in order to overcome artefact of having 
  #same exit and entry date from splitting at first positive dates, and birthdates
  summary(lx)
  summary(lx, by='hivstale5y', Rates = T, scale = 1000)
  #print(summary(lx, by='hivstale5y', Rates = T), digits=4)
  
  
  if(FALSE){
    View(lx[lx$lex.id==67180,])
    st <- survtab(Surv(time = lex.dur, event = lex.Xst) ~ hivstale5y, data = lx, 
                  surv.type = "surv.obs",
                  breaks = list(lex.dur = seq(15, 100, 1)))
  }
  
  
  #du <- cal.yr(dat$exit)-cal.yr(dat$entry)
  #sum(du==0) #4486
  #by(du, dat$hivstale5y, sum)
  #c(97645.92,79372.85,362437.1)/365
  #c(952,2329,4285)/(c(97645.92,79372.85,362437.1)/365)
  
  #records with same entry and exit dates
  #samedates <- dat[dat$exit==dat$entry,]
  #samedates <- dat[which(du==0),]
  #View(samedates)
  #length(unique(samedates$idno_original)) #1525
  #sum(!is.na(samedates$COD)) #of which 372 has CoD
  
  #using only accident CoD
  lx2 <- Lexis(entry = list(per=cal.yr(entry)), exit = list(per=cal.yr(exit)+.000001), id=idno_original, exit.status=fail2, data=dat)
  #saveRDS(lx,'lexisData20180719.RDS')
  #lx <- readRDS('lexisData20180719.RDS')
  #+.000001 was added in exit year in order to overcome artefact of having 
  #same exit and entry date from splitting at first positive dates, and birthdates
  summary(lx2)
  summary(lx2, by='hivstale5y', Rates = T, scale = 1000)
}

#others
if(FALSE){
  #the rest# delete after a while
  #plot(survfit(Surv(time0, timex, fail0) ~ 1, data=dat), conf.int=FALSE, mark.time=FALSE) #no comparison
  #sv <- with(dat,Surv(time0, timex, fail0))
  sv <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y, data=dat) #errors
  sv
  summary(sv)
  summary(sv, times = c(20, 30, 40))
  
  svRes <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ residence, data=dat) #errors
  svRes
  plot(svRes, conf.int=FALSE, mark.time=FALSE, col=1:3)
  legend('bottomleft',legend = c("urban","peri-urban", "rural"), col=1:3, lty=1)
  
  
  svSex <- survfit(Surv(time=time0, time2 = timex, event = fail0) ~ sex, data=dat) #errors
  svSex
  plot(svSex, conf.int=FALSE, mark.time=FALSE, col=1:2)
  legend('bottomleft',legend = c("men","women"), col=1:2, lty=1)
  
  
  if(FALSE){
    #this method uses each line as a single case, which is not true
    #person years
    py <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ 1, data=dat, scale = 1)
    py_hiv <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y, data=dat, scale = 1)
    summary(py_hiv, rate = T, ci.r = T)
    
    py_hiv2 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y, data=dat, scale = 1)
    summary(py_hiv2, rate = T, ci.r = T)
    
    py_hiv2 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat, scale = 1)
    summary(py_hiv2, rate = T, ci.r = T)
    
    py_hiv2 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ seroconFixed, data=dat, scale = 1)
    summary(py_hiv2, rate = T, ci.r = T)
    
  }
  
  plot(svcox, conf.int=FALSE, mark.time=FALSE, col=1:3, ymin = .8) #by HIV status, DOESN'T plot anything
  legend('bottomleft',legend = c('neg','pos','unk'), col=1:3, lty=1)
  
  svcoxid <- coxph(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y+cluster(idno_original), data=dat) #errors
  svcoxx <- coxph(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y+sex, data=dat) #errors
  ggsurvplot(svcox)
  survdiff(Surv(time=time0, time2 = timex, event = fail0) ~ hivstale5y, data=dat)
  
}

#Further checks before the analyses are here####
if(FALSE){
  #check age consistency
  #summarize(dat$age,by(dat$agegrp), max)
  by(liteumk$age,liteumk$agegrp,summary) #all age groups are consistent
  by(liteumk$age,liteumk$agegrp,function(x){sum(is.na(x))})
  
  #residence variable doesn't have data on all episodes
  sum(c(by(dat$residence,liteumk$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T) 
  #cases with at least 1 missing residence 54855, not 48725
  
  #sum(c(by(dat$residence2,dat$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T)
  #sum(c(by(datL$residence2,datL$idno_original,function(x){sum(is.na(x))}))>=1, na.rm=T)
  
  #residence2: 16017, not 14748
  
  #checking residence and residence2
  # > sum(c(6317, 18997, 35766, 27613))
  # [1] 88693
  # > sum(c(7621, 22863, 42192, 16017 ))
  # [1] 88693
  
  #removing missing values #obsolete now
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
  
}

