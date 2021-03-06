#sp_lifetable, started 20180729
#all in agegrp = 15
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- FALSE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')
injNames <- readRDS('accidentNames.RDS')

#putting new variables for each injury
for(i in 1: length(injNames)){
  x <- NA
  x <- (dat$fail2==1)*(dat$COD %in% injNames[i])
  dat <- cbind(dat,x)
  colnames(dat)[ncol(dat)] <- injNames[i]
}
#1271372      57

#calander year 2007 until "2015-11-26", it would be a period LT with a long period####
#last exit on "2015-12-01"
#unabridged lifetables####
#all cause
asmr.all <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat, scale = 1)
#summary(asmr.all, rate = T, ci.r = T)
#lt(asmr.all, per=10000, ci95=T)
lt.all <- lt(asmr.all, ageint = 1)
#if(creation) write.csv(lt.all,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_all.csv",sep = ""))

#ext.inj cause
asmr.inj <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat, scale = 1)
lt.inj <- lt(asmr.inj, ageint = 1)
#if(creation) write.csv(lt.inj,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_inj.csv",sep = ""))

if(FALSE){
#ext. injury deleted lifetable####
#overall
asdt.inj.del <- asdt(allcause = lt.all, i_cause = lt.inj, ageint = 1)
if(creation) write.csv(asdt.inj.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_del.csv",sep = ""))
asdt.inj.del$s_ex[1]-asdt.inj.del$ex[1] #1.84 years

#all cause by 15 year age group
asmr.all.15 <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat, scale = 5)
lt.all.15 <- lt(asmr.all.15, ageint = 15)
asmr.inj.15 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat, scale = 5)
lt.inj.15 <- lt(asmr.inj.15, ageint = 15)
asdt.inj.del <- asdt(allcause = lt.all.15, i_cause = lt.inj.15, ageint = 15)
asdt.inj.del$s_ex[1]-asdt.inj.del$ex[1] #1.868421 years

}
#e15 gain, comparison within itself in HIV group, 15 year agegrp15####
#hivstatus allFixed has beeen used
#1 year age group doesn't work as data is inadequate: producing NaN
#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Women.Negative.del <- asdt(allcause = asmr.all.Women.Negative, i_cause = asmr.inj.Women.Negative, ageint = 15)
if(creation) write.csv(asdt.inj.Women.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Neg_del.csv",sep = ""))
Women.Negative <- asdt.inj.Women.Negative.del$s_ex[1]-asdt.inj.Women.Negative.del$ex[1] 
Women.Negative # 1.860363, previously NaN in 1 year age group
#Positive
asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Women.Positive.del <- asdt(allcause = asmr.all.Women.Positive, i_cause = asmr.inj.Women.Positive, ageint = 15)
if(creation) write.csv(asdt.inj.Women.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Pos_del.csv",sep = ""))
Women.Positive <- asdt.inj.Women.Positive.del$s_ex[1]-asdt.inj.Women.Positive.del$ex[1] 
Women.Positive #0.6764863, previously NaN in 1 year age group
#Unknown
asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Women.Unknown.del <- asdt(allcause = asmr.all.Women.Unknown, i_cause = asmr.inj.Women.Unknown, ageint = 15)
if(creation) write.csv(asdt.inj.Women.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Unk_del.csv",sep = ""))
Women.Unknown <- asdt.inj.Women.Unknown.del$s_ex[1]-asdt.inj.Women.Unknown.del$ex[1] 
Women.Unknown # 0.917153 for 5 year age group #0.684899 for 1 year age group

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Men.Negative.del <- asdt(allcause = asmr.all.Men.Negative, i_cause = asmr.inj.Men.Negative, ageint = 15)
if(creation) write.csv(asdt.inj.Men.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Neg_del.csv",sep = ""))
Men.Negative <- asdt.inj.Men.Negative.del$s_ex[1]-asdt.inj.Men.Negative.del$ex[1] 
Men.Negative # 3.535453
#Positive
asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Men.Positive.del <- asdt(allcause = asmr.all.Men.Positive, i_cause = asmr.inj.Men.Positive, ageint = 15)
if(creation) write.csv(asdt.inj.Men.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Pos_del.csv",sep = ""))
Men.Positive <- asdt.inj.Men.Positive.del$s_ex[1]-asdt.inj.Men.Positive.del$ex[1] 
Men.Positive #1.537419
#Unknown
asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Men.Unknown.del <- asdt(allcause = asmr.all.Men.Unknown, i_cause = asmr.inj.Men.Unknown, ageint = 15)
if(creation) write.csv(asdt.inj.Men.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Unk_del.csv",sep = ""))
Men.Unknown <- asdt.inj.Men.Unknown.del$s_ex[1]-asdt.inj.Men.Unknown.del$ex[1] 
Men.Unknown #3.309036

#Decomposition by age, sex and CoD from injury####
#Women
decom.Women.Pos.Neg <- decom(allcause.A=asmr.all.Women.Positive,allcause.B =  asmr.all.Women.Negative, i_cause.A = asmr.inj.Women.Positive,i_cause.B =  asmr.inj.Women.Negative, ageint=1)
if(creation) write.csv(decom.Women.Pos.Neg, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decom_Women_Pos_Neg.csv",sep = ""))
#sum(decom.Women.Pos.Neg$ndeltax)
#sum(decom.Women.Pos.Neg$ndeltax.i)
decom.Women.Unk.Neg <- decom(allcause.A=asmr.all.Women.Unknown,allcause.B =  asmr.all.Women.Negative, i_cause.A = asmr.inj.Women.Unknown,i_cause.B =  asmr.inj.Women.Negative, ageint=1)
if(creation) write.csv(decom.Women.Unk.Neg, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decom_Women_Unk_Neg.csv",sep = ""))


#Men
decom.Men.Pos.Neg <- decom(allcause.A=asmr.all.Men.Positive,allcause.B =  asmr.all.Men.Negative, i_cause.A = asmr.inj.Men.Positive,i_cause.B =  asmr.inj.Men.Negative, ageint=1)
if(creation) write.csv(decom.Men.Pos.Neg, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decom_Men_Pos_Neg.csv",sep = ""))
#sum(decom.Men.Pos.Neg$ndeltax)
#sum(decom.Men.Pos.Neg$ndeltax.i)
decom.Men.Unk.Neg <- decom(allcause.A=asmr.all.Men.Unknown,allcause.B =  asmr.all.Men.Negative, i_cause.A = asmr.inj.Men.Unknown,i_cause.B =  asmr.inj.Men.Negative, ageint=1)
if(creation) write.csv(decom.Men.Unk.Neg, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decom_Men_Unk_Neg.csv",sep = ""))


#lapply this
#pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat, scale = 1)$event
#lapply(injNames,function(x){x})
# injEventList <- lapply(injNames,function(x){
#   pyears(Surv(time=time0, time2 = timex, event = dat[,which(colnames(dat)==x)]) ~ agegrp15, data=dat, scale = 1)$event
# })
#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns####
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Women
decomList.Women.Positive.Negative <- decomList(allcause.A = asmr.all.Women.Positive, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Positive, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Positive.Negative) <- c(colnames(decomList.Women.Positive.Negative)[1:7],injNames)
#decomList.Women.Positive.Negative
if(creation) write.csv(decomList.Women.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_Negative.csv",sep = ""))

decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative.csv",sep = ""))

#Men
decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
#decomList.Men.Positive.Negative
if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative.csv",sep = ""))

decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative.csv",sep = ""))

#decomposing within periods####
dat.original <- dat

#2007-2010####
dat <- dat.original[dat.original$period.2011_15==FALSE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2007 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Women.Negative.del <- asdt(allcause = asmr.all.Women.Negative, i_cause = asmr.inj.Women.Negative, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Neg_del.csv",sep = ""))
Women.Negative <- asdt.inj.Women.Negative.del$s_ex[1]-asdt.inj.Women.Negative.del$ex[1] 
Women.Negative # 1.860363, previously NaN in 1 year age group
#Positive
asmr.all.Women.Positive.2007 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Women.Positive.del <- asdt(allcause = asmr.all.Women.Positive, i_cause = asmr.inj.Women.Positive, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Pos_del.csv",sep = ""))
Women.Positive <- asdt.inj.Women.Positive.del$s_ex[1]-asdt.inj.Women.Positive.del$ex[1] 
Women.Positive #0.6764863, previously NaN in 1 year age group
#Unknown
asmr.all.Women.Unknown.2007 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Women.Unknown.del <- asdt(allcause = asmr.all.Women.Unknown, i_cause = asmr.inj.Women.Unknown, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Unk_del.csv",sep = ""))
Women.Unknown <- asdt.inj.Women.Unknown.del$s_ex[1]-asdt.inj.Women.Unknown.del$ex[1] 
Women.Unknown # 0.917153 for 5 year age group #0.684899 for 1 year age group

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2007 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Men.Negative.del <- asdt(allcause = asmr.all.Men.Negative, i_cause = asmr.inj.Men.Negative, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Neg_del.csv",sep = ""))
Men.Negative <- asdt.inj.Men.Negative.del$s_ex[1]-asdt.inj.Men.Negative.del$ex[1] 
Men.Negative # 3.535453
#Positive
asmr.all.Men.Positive.2007 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Men.Positive.del <- asdt(allcause = asmr.all.Men.Positive, i_cause = asmr.inj.Men.Positive, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Pos_del.csv",sep = ""))
Men.Positive <- asdt.inj.Men.Positive.del$s_ex[1]-asdt.inj.Men.Positive.del$ex[1] 
Men.Positive #1.537419
#Unknown
asmr.all.Men.Unknown.2007 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Men.Unknown.del <- asdt(allcause = asmr.all.Men.Unknown, i_cause = asmr.inj.Men.Unknown, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Unk_del.csv",sep = ""))
Men.Unknown <- asdt.inj.Men.Unknown.del$s_ex[1]-asdt.inj.Men.Unknown.del$ex[1] 
Men.Unknown #3.309036

#2007-2010#
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive.2007 <- injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative.2007 <- injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown.2007 <- injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive.2007 <- injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative.2007 <- injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown.2007 <- injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Women
decomList.Women.Positive.Negative <- decomList(allcause.A = asmr.all.Women.Positive, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Positive, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Positive.Negative) <- c(colnames(decomList.Women.Positive.Negative)[1:7],injNames)
#decomList.Women.Positive.Negative
if(creation) write.csv(decomList.Women.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_Negative_2007.csv",sep = ""))

decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative_2007.csv",sep = ""))

#Men
decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
#decomList.Men.Positive.Negative
if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative_2007.csv",sep = ""))

decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative_2007.csv",sep = ""))

#2011-2015####
dat <- dat.original[dat.original$period.2011_15==TRUE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2011 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Women.Negative.del <- asdt(allcause = asmr.all.Women.Negative, i_cause = asmr.inj.Women.Negative, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Neg_del.csv",sep = ""))
Women.Negative <- asdt.inj.Women.Negative.del$s_ex[1]-asdt.inj.Women.Negative.del$ex[1] 
Women.Negative # 1.860363, previously NaN in 1 year age group
#Positive
asmr.all.Women.Positive.2011 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Women.Positive.del <- asdt(allcause = asmr.all.Women.Positive, i_cause = asmr.inj.Women.Positive, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Pos_del.csv",sep = ""))
Women.Positive <- asdt.inj.Women.Positive.del$s_ex[1]-asdt.inj.Women.Positive.del$ex[1] 
Women.Positive #0.6764863, previously NaN in 1 year age group
#Unknown
asmr.all.Women.Unknown.2011 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Women.Unknown.del <- asdt(allcause = asmr.all.Women.Unknown, i_cause = asmr.inj.Women.Unknown, ageint = 15)
#if(creation) write.csv(asdt.inj.Women.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Unk_del.csv",sep = ""))
Women.Unknown <- asdt.inj.Women.Unknown.del$s_ex[1]-asdt.inj.Women.Unknown.del$ex[1] 
Women.Unknown # 0.917153 for 5 year age group #0.684899 for 1 year age group

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2011 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asdt.inj.Men.Negative.del <- asdt(allcause = asmr.all.Men.Negative, i_cause = asmr.inj.Men.Negative, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Neg_del.csv",sep = ""))
Men.Negative <- asdt.inj.Men.Negative.del$s_ex[1]-asdt.inj.Men.Negative.del$ex[1] 
Men.Negative # 3.535453
#Positive
asmr.all.Men.Positive.2011 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asdt.inj.Men.Positive.del <- asdt(allcause = asmr.all.Men.Positive, i_cause = asmr.inj.Men.Positive, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Pos_del.csv",sep = ""))
Men.Positive <- asdt.inj.Men.Positive.del$s_ex[1]-asdt.inj.Men.Positive.del$ex[1] 
Men.Positive #1.537419
#Unknown
asmr.all.Men.Unknown.2011 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asdt.inj.Men.Unknown.del <- asdt(allcause = asmr.all.Men.Unknown, i_cause = asmr.inj.Men.Unknown, ageint = 15)
#if(creation) write.csv(asdt.inj.Men.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Unk_del.csv",sep = ""))
Men.Unknown <- asdt.inj.Men.Unknown.del$s_ex[1]-asdt.inj.Men.Unknown.del$ex[1] 
Men.Unknown #3.309036

#2011-2015#
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive.2011 <- injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative.2011 <- injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown.2011 <- injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive.2011 <- injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative.2011 <- injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown.2011 <- injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#Women
decomList.Women.Positive.Negative <- decomList(allcause.A = asmr.all.Women.Positive, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Positive, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Positive.Negative) <- c(colnames(decomList.Women.Positive.Negative)[1:7],injNames)
#decomList.Women.Positive.Negative
if(creation) write.csv(decomList.Women.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_Negative_2011.csv",sep = ""))

decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 15)
colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative_2011.csv",sep = ""))

#Men
decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
#decomList.Men.Positive.Negative
if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative_2011.csv",sep = ""))

decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative_2011.csv",sep = ""))

dat <- dat.original #restoring#

#decomposing between periods####
#POSITIVES
decomList.Women.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Women.Positive.2007, allcause.B = asmr.all.Women.Positive.2011, i_cause.A=injEventList.Women.Positive.2007, i_cause.B = injEventList.Women.Positive.2011, ageint = 15)
colnames(decomList.Women.Positive.2007.2011) <- c(colnames(decomList.Women.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_2011_2007.csv",sep = ""))

decomList.Men.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Men.Positive.2007, allcause.B = asmr.all.Men.Positive.2011, i_cause.A=injEventList.Men.Positive.2007, i_cause.B = injEventList.Men.Positive.2011, ageint = 15)
colnames(decomList.Men.Positive.2007.2011) <- c(colnames(decomList.Men.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_2011_2007.csv",sep = ""))

#NEGATIVES
decomList.Women.Negative.2007.2011 <- decomList(allcause.A = asmr.all.Women.Negative.2007, allcause.B = asmr.all.Women.Negative.2011, i_cause.A=injEventList.Women.Negative.2007, i_cause.B = injEventList.Women.Negative.2011, ageint = 15)
colnames(decomList.Women.Negative.2007.2011) <- c(colnames(decomList.Women.Negative.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.Negative.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Negative_2011_2007.csv",sep = ""))

decomList.Men.Negative.2007.2011 <- decomList(allcause.A = asmr.all.Men.Negative.2007, allcause.B = asmr.all.Men.Negative.2011, i_cause.A=injEventList.Men.Negative.2007, i_cause.B = injEventList.Men.Negative.2011, ageint = 15)
colnames(decomList.Men.Negative.2007.2011) <- c(colnames(decomList.Men.Negative.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.Negative.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Negative_2011_2007.csv",sep = ""))

#UNKNOWNS
decomList.Women.Unknown.2007.2011 <- decomList(allcause.A = asmr.all.Women.Unknown.2007, allcause.B = asmr.all.Women.Unknown.2011, i_cause.A=injEventList.Women.Unknown.2007, i_cause.B = injEventList.Women.Unknown.2011, ageint = 15)
colnames(decomList.Women.Unknown.2007.2011) <- c(colnames(decomList.Women.Unknown.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.Unknown.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_2011_2007.csv",sep = ""))

decomList.Men.Unknown.2007.2011 <- decomList(allcause.A = asmr.all.Men.Unknown.2007, allcause.B = asmr.all.Men.Unknown.2011, i_cause.A=injEventList.Men.Unknown.2007, i_cause.B = injEventList.Men.Unknown.2011, ageint = 15)
colnames(decomList.Men.Unknown.2007.2011) <- c(colnames(decomList.Men.Unknown.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.Unknown.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_2011_2007.csv",sep = ""))
