#summer project
#essential info to put in report
lfc <- 1 #left censoring at '2007-01-01'
rtc <- 1 #right censoring for age above 100

library(Rcpp)
library("readstata13")
library("Hmisc")
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(Epi)
library(popEpi)
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
table(liteumk$fail0, exclude=NULL)
#table(liteumk$fail0, liteumk$failure, exclude=NULL) #failure var removed
#Because of this resplitting: ‘fail’ variable is now useless as they have been further multiplicated 
#by the split ie. some may have multiple 1 on the ‘fail’ variable

#no. of deaths which had va (main dataset)
table(liteumk$fail0, liteumk$hadva, exclude=NULL)
#2572+13400=15972
#2572 didn't have va
#13400 had va

#tmp <- liteumk[which(is.na(liteumk$fail0)),]

#datacleaning
liteumk <- liteumk[which(!is.na(liteumk$fail0)),] #remove missing from fail0 [already described in report]
#additional cases with same entry and exit dates found out later [on date: 20180719] However, see below
#liteumk <- liteumk[which(liteumk$entry!=liteumk$exit),] #this is an artefact of spliting by birthday/positive date
length(unique(liteumk[which(is.na(liteumk$fail0)),]$idno_original)) #2034 cases only instead of 2037
#the discrepancy in # is because some cases have records with NA in fail0 (ie. same entry & exit date)

#cases after removal of missing outcomes
length(unique(liteumk$idno_original))

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
dim(liteumk) #1271372 x36

dat <- merge(liteumk,CODdataset, by="idno_original", all.x = TRUE)

dim(dat) #1271372 x37
sum(!is.na(dat$COD)) #58222 #40185

#no. of deaths
dim(dat[dat$fail0==1,]) #7566 deaths
#sum(c(by(dat$fail0, dat$idno_original, function(x){sum(x)>=1}))) #same as above

#sum(c(by(dat$fail0, dat$idno_original, function(x){sum(x)>1}))) #checking if death in a person is recorded twice, 0

#no. of deaths which had VA
dim(dat[dat$fail0==1 & !is.na(dat$COD),]) #6222 #4530

####Adding new fail variable 'fail2'####
#loop, test 2 conditions (liteumk$fail==1 & liteumk$COD %in% accidentNames)
accidentNames <- readRDS('accidentNames.RDS')
dim(dat) #1271372      37
if(FALSE){
  #THIS TAKES TOO LONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  for(i in 1:nrow(dat)){
    if(dat$fail[i]==1 & dat$COD[i] %in% accidentNames) dat$fail2[i] <- 1
    else dat$fail2[i] <- 0
  }
}
#another faster way of doing the above
fail2 <- (dat$fail==1)*(dat$COD %in% accidentNames)
dat <- cbind(dat,fail2)
dim(dat) #1271372      38
#sum(is.na(dat$fail2))
table(dat$fail2, exclude = NULL) #only 817 deaths from accidents

#tmp1 <- dat[dat$fail2==1,]
#tmp2 <- dat[fail2fast==1,]
#identical(tmp1, tmp2) #TRUE

#further cleaning done here
#the aim is to describe only to describe the final (FINAL) dataset
#the first state of each individual when first enrolled???

#dat <- dat[which(!is.na(dat$COD)),]

#subsetting the LAST records
dat <- dat[order(dat$entry, decreasing = T),] #LAST record
#dat <- dat[order(dat$exit, decreasing = T),]
datL <- dat[!duplicated(dat$idno_original),]

####Descriptive statistics of the last records####
varDes <- c('sex','residence','hadva','failure','agegrp',
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

#descriptives <- do.call(rbind.data.frame, des)

#get hivstatus_broad missing to unknown
dat$hivstatus_broad[is.na(dat$hivstatus_broad)] <- 'Unknown'
names(dat)[names(dat)=="_t0"] <- "entry"
names(dat)[names(dat)=="_t"] <- "exit"
names(dat)[names(dat)==""] <- "fail0"

####Survival Analysis####
lx <- Lexis(entry = list(per=cal.yr(entry)), exit = list(per=cal.yr(exit)+.000001), id=idno_original, exit.status=fail0, data=dat)
#saveRDS(lx,'lexisData20180719.RDS')
#lx <- readRDS('lexisData20180719.RDS')
#+.000001 was added in exit year in order to overcome artefact of having 
#same exit and entry date from splitting at first positive dates, and birthdates
summary(lx)
summary(lx, by='hivstatus_broad', Rates = T, scale = 1000)
#print(summary(lx, by='hivstatus_broad', Rates = T), digits=4)

View(lx[lx$lex.id==67180,])
st <- survtab(Surv(time = lex.dur, event = lex.Xst) ~ hivstatus_broad, data = lx, 
              surv.type = "surv.obs",
              breaks = list(lex.dur = seq(15, 100, 1)))

#du <- cal.yr(dat$exit)-cal.yr(dat$entry)
#sum(du==0) #4486
by(du, dat$hivstatus_broad, sum)
c(97645.92,79372.85,362437.1)/365
c(952,2329,4285)/(c(97645.92,79372.85,362437.1)/365)

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
summary(lx2, by='hivstatus_broad', Rates = T, scale = 1000)

sv2 <- survfit(Surv(time=entry, time2 = exit, event = fail2) ~ hivstatus_broad, data=dat) #errors
summary(sv2)
summary(sv2, times = c(20, 30, 40))
plot(sv2, conf.int=FALSE, mark.time=FALSE, col=1:3) #by HIV status
legend('bottomleft',legend = c('neg','pos','unk'), col=1:3, lty=1)

#the rest# delete after a while
#plot(survfit(Surv(entry, exit, fail0) ~ 1, data=dat), conf.int=FALSE, mark.time=FALSE) #no comparison
#sv <- with(dat,Surv(entry, exit, fail0))
sv <- survfit(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad, data=dat) #errors
sv
summary(sv)
summary(sv, times = c(20, 30, 40))

sv3 <- survfit(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad+cluster(idno_original), data=dat) #errors
sv3

if(FALSE){
  #this method uses each line as a single case, which is not true
  py <- pyears(Surv(time=entry, time2 = exit, event = fail0) ~ 1, data=dat, scale = 1)
  py_hiv <- pyears(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad, data=dat, scale = 1)
  
}

svcox <- coxph(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad, data=dat) #errors
svcoxid <- coxph(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad+cluster(idno_original), data=dat) #errors
svcoxx <- coxph(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad+sex, data=dat) #errors
ggsurvplot(svcox)
survdiff(Surv(time=entry, time2 = exit, event = fail0) ~ hivstatus_broad, data=dat)

plot(sv, conf.int=FALSE, mark.time=FALSE, col=1:3) #by HIV status
legend('topright',legend = c('neg','pos','unk'), col=1:3, lty=1)

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

