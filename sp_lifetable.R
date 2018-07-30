#sp_lifetable, started 20180729
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- FALSE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')

#calander year 2007 until "2015-11-26", it would be a period LT with a long period####
#all cause
asmr.all <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat, scale = 1)
#summary(asmr.all, rate = T, ci.r = T)
#lt(asmr.all, per=10000, ci95=T)
lt.all <- lt(asmr.all, ageint = 1)
if(creation) write.csv(lt.all,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_all.csv",sep = ""))

#ext.inj cause
asmr.inj <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat, scale = 1)
lt.inj <- lt(asmr.inj, ageint = 1)
if(creation) write.csv(lt.inj,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_inj.csv",sep = ""))

#ext. injury deleted lifetable####
#overall
asdt.inj.del <- asdt(allcause = lt.all, i_cause = lt.inj, ageint = 1)
if(creation) write.csv(asdt.inj.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_del.csv",sep = ""))
asdt.inj.del$s_ex[1]-asdt.inj.del$ex[1] #1.84 years

#all cause by 5 year age group
asmr.all.5 <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat, scale = 5)
lt.all.5 <- lt(asmr.all.5, ageint = 5)
asmr.inj.5 <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat, scale = 5)
lt.inj.5 <- lt(asmr.inj.5, ageint = 5)
asdt.inj.del <- asdt(allcause = lt.all.5, i_cause = lt.inj.5, ageint = 5)
asdt.inj.del$s_ex[1]-asdt.inj.del$ex[1] #1.868421 years

#comparison within itself in HIV group####
#hivstatus allFixed has beeen used
#1 year age group doesn't work as data is inadequate: producing NaN
#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
asdt.inj.Women.Negative.del <- asdt(allcause = asmr.all.Women.Negative, i_cause = asmr.inj.Women.Negative, ageint = 5)
if(creation) write.csv(asdt.inj.Women.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Neg_del.csv",sep = ""))
Women.Negative <- asdt.inj.Women.Negative.del$s_ex[1]-asdt.inj.Women.Negative.del$ex[1] 
Women.Negative # 1.860363, previously NaN in 1 year age group
#Positive
asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
asdt.inj.Women.Positive.del <- asdt(allcause = asmr.all.Women.Positive, i_cause = asmr.inj.Women.Positive, ageint = 5)
if(creation) write.csv(asdt.inj.Women.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Pos_del.csv",sep = ""))
Women.Positive <- asdt.inj.Women.Positive.del$s_ex[1]-asdt.inj.Women.Positive.del$ex[1] 
Women.Positive #0.6764863, previously NaN in 1 year age group
#Unknown
asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)
asdt.inj.Women.Unknown.del <- asdt(allcause = asmr.all.Women.Unknown, i_cause = asmr.inj.Women.Unknown, ageint = 5)
if(creation) write.csv(asdt.inj.Women.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Women_Unk_del.csv",sep = ""))
Women.Unknown <- asdt.inj.Women.Unknown.del$s_ex[1]-asdt.inj.Women.Unknown.del$ex[1] 
Women.Unknown # 0.917153 for 5 year age group #0.684899 for 1 year age group

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
asdt.inj.Men.Negative.del <- asdt(allcause = asmr.all.Men.Negative, i_cause = asmr.inj.Men.Negative, ageint = 5)
if(creation) write.csv(asdt.inj.Men.Negative.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Neg_del.csv",sep = ""))
Men.Negative <- asdt.inj.Men.Negative.del$s_ex[1]-asdt.inj.Men.Negative.del$ex[1] 
Men.Negative # 3.535453
#Positive
asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
asdt.inj.Men.Positive.del <- asdt(allcause = asmr.all.Men.Positive, i_cause = asmr.inj.Men.Positive, ageint = 5)
if(creation) write.csv(asdt.inj.Men.Positive.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Pos_del.csv",sep = ""))
Men.Positive <- asdt.inj.Men.Positive.del$s_ex[1]-asdt.inj.Men.Positive.del$ex[1] 
Men.Positive #1.537419
#Unknown
asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)
asdt.inj.Men.Unknown.del <- asdt(allcause = asmr.all.Men.Unknown, i_cause = asmr.inj.Men.Unknown, ageint = 5)
if(creation) write.csv(asdt.inj.Men.Unknown.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_Men_Unk_del.csv",sep = ""))
Men.Unknown <- asdt.inj.Men.Unknown.del$s_ex[1]-asdt.inj.Men.Unknown.del$ex[1] 
Men.Unknown #3.309036

Women.e15.increase <- c(Women.Negative, Women.Positive, Women.Unknown)
Men.e15.increase <- c(Men.Negative, Men.Positive, Men.Unknown)
e15.increase <- cbind(Women.e15.increase, Men.e15.increase)
if(creation) write.csv(e15.increase, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_e15_increase.csv",sep = ""))


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
