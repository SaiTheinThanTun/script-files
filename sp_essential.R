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
  liteumk <- liteumk[which((liteumk$entry>='2007-01-01')),] 
}

#cases after left censoring
length(unique(liteumk$idno_original)) #88693

####Right censoring (change terms) for records with age>100####
if(rtc){
  liteumk <- liteumk[which((liteumk$age<=100)),] 
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
#fail definitions
#STATA defined death, 'failure' is the variable STATA used for this. there's also 'fail' which is not used
#fail2 is external injury related death

#saveRDS(dat, 'dat_lifetable.RDS') #processed data file to be used for lifetable


####sex, age, py distributions####
desTableNames <- c('Individuals','Person-years','All cause: deaths', 'Rate', 'Low 95%CI', 'High 95%CI','Ext. Inj: deaths', 'Rate', 'Low 95%CI', 'High 95%CI')
if(FALSE){
per <- 10000

no.individuals <- table(datL$agegrp,datL$sex,  exclude = NULL)
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

f.dist <- cbind(no.individuals[,2],round(pyaa.f.py),pyaa.f.event,round(pyaa.f.rate*per,2),round(pyaa.f.lo*per,2),round(pyaa.f.hi*per,2),pyae.f.event,round(pyae.f.rate*per,2),round(pyae.f.lo*per,2),round(pyae.f.hi*per,2))

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
Total <- c(tot.indi,round(tot.pyaa.f.py),tot.pyaa.f.event,round(tot.pyaa.f.rate*per,2),round(tot.pyaa.f.lo*per,2),round(tot.pyaa.f.hi*per,2),tot.pyae.f.event,round(tot.pyae.f.rate*per,2),round(tot.pyae.f.lo*per,2),round(tot.pyae.f.hi*per,2))

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

m.dist <- cbind(no.individuals[,1],round(pyaa.m.py),pyaa.m.event,round(pyaa.m.rate*per,2),round(pyaa.m.lo*per,2),round(pyaa.m.hi*per,2),pyae.m.event,round(pyae.m.rate*per,2),round(pyae.m.lo*per,2),round(pyae.m.hi*per,2))

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
Total <- c(tot.indi,round(tot.pyaa.m.py),tot.pyaa.m.event,round(tot.pyaa.m.rate*per,2),round(tot.pyaa.m.lo*per,2),round(tot.pyaa.m.hi*per,2),tot.pyae.m.event,round(tot.pyae.m.rate*per,2),round(tot.pyae.m.lo*per,2),round(tot.pyae.m.hi*per,2))

m.dist <- rbind(m.dist, Total)
colnames(m.dist) <- desTableNames
m.dist
#check result
py.m <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ sex, data=dat, scale = 1)
summary(py.m,rate = T, ci.r = T)
if(creation) write.csv(m.dist,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_Men.csv",sep = "") )


#descriptions on HIV status
#original rate
hivstatus_broad.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ hivstatus_broad, data=dat, scale = 1)
hivstatus_broad.rate <- pyears2(hivstatus_broad.rate, per = 10000)
if(creation) write.csv(hivstatus_broad.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_hivstatus_broad_rate.csv",sep = "") )

#neg for 5 years
hivstale5y.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y, data=dat, scale = 1)
hivstale5y.rate <- pyears2(hivstale5y.rate, per = 10000)
if(creation) write.csv(hivstale5y.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_hivstale5y_rate.csv",sep = "") )

#neg for 5 years + missing data fixed
missingFixed.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ missingFixed, data=dat, scale = 1)
missingFixed.rate <- pyears2(missingFixed.rate, per = 10000)
if(creation) write.csv(missingFixed.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_missingFixed_rate.csv",sep = "") )

#neg for 5 years + seroconverter data fixed
seroconFixed.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ seroconFixed, data=dat, scale = 1)
seroconFixed.rate <- pyears2(seroconFixed.rate, per = 10000)
if(creation) write.csv(seroconFixed.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_seroconFixed_rate.csv",sep = "") )

#neg for 5 years + all fixed
allFixed.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat, scale = 1)
allFixed.rate <- pyears2(allFixed.rate, per = 10000)
if(creation) write.csv(allFixed.rate,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_allFixed_rate.csv",sep = "") )

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

#by hivstatus exposure, stratified by sex
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

###COX REGRESSION####
#!!need to check proportionality of hazards
#basic
#svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ factor(hivstale5y), data=dat) 
#summary(svcox)
testOC <- testOClist[5] #c("hivstatus_broad","hivstale5y", "missingFixed", "seroconFixed", "allFixed")
svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)]), data=dat) 
summary(svcox)
base.ph <- cox.zph(svcox)
plot(base.ph) #TEST PROPORTIONALITY
#cumulative hazard
plot(survfit(Surv(time=time0, time2 = timex, event = fail0) ~ allFixed, data = dat),fun='cloglog',xlab = 'Years', ylab="Cumulative Hazard (log)")
plot(survfit(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data = dat),fun='cloglog',xlab = 'Years', ylab="Cumulative Hazard (log)",lty=1:3)

po_con <- c('sex','age','residence')
#? timeposneg??? > Negative test expiry date in years
#? hivevertreat???
# NA assumed to be Unknown

#basic+res *
#res.svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y+residence, data=dat) 
#summary(res.svcox)
res.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(residence), data=dat) 
summary(res.svcox)
res.ph <- cox.zph(res.svcox)
plot(res.ph) #TEST PROPORTIONALITY

#basic+sex *
#sex.svcox <- coxph(Surv(time=time0, time2 = timex, event = fail2) ~ hivstale5y+sex, data=dat) 
#summary(sex.svcox)
sex.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+factor(sex), data=dat) 
summary(sex.svcox)
sex.ph <- cox.zph(sex.svcox)
plot(sex.ph) #TEST PROPORTIONALITY

#basic+age *?
#age.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(hivstale5y)+factor(agegrp), data=dat) #error
age.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age, data=dat)
summary(age.svcox)
age.ph <- cox.zph(age.svcox)
plot(age.ph) #TEST PROPORTIONALITY


#basic+age+sex
as.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(sex), data=dat) 
summary(as.svcox)
as.ph <- cox.zph(as.svcox)
as.ph
plot(as.ph) #TEST PROPORTIONALITY

#basic+age+res+sex
ars.svcox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(dat[,which(names(dat) %in% testOC)])+age+factor(residence)+factor(sex), data=dat) 
summary(ars.svcox)
ars.ph <- cox.zph(ars.svcox)
ars.ph
plot(ars.ph) #TEST PROPORTIONALITY


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

