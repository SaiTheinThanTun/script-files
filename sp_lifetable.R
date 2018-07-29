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
asmr.ext <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat, scale = 1)
lt.ext <- lt(asmr.ext, ageint = 1)
if(creation) write.csv(lt.ext,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_ext.csv",sep = ""))

#all cause by 5 year age group
asmr.all.5 <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat, scale = 1)
lt(asmr.all.5, ageint = 5)

#ext. injury deleted lifetable####
#overall
asdt.inj.del <- asdt(allcause = lt.all, i_cause = lt.ext, ageint = 1)
if(creation) write.csv(asdt.inj.del,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_asdt_inj_del.csv",sep = ""))
asdt.inj.del$s_ex[1]-asdt.inj.del$ex[1] #1.84 years

#comparison between HIV groups####
#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 1)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 1)
asdt.inj.Women.Negative.del <- asdt(allcause = asmr.all.Women.Negative, i_cause = asmr.inj.Women.Negative, ageint = 1)
Women.Negative <- asdt.inj.Women.Negative.del$s_ex[1]-asdt.inj.Women.Negative.del$ex[1] 
Women.Negative
#Positive
asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 1)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 1)
asdt.inj.Women.Positive.del <- asdt(allcause = asmr.all.Women.Positive, i_cause = asmr.inj.Women.Positive, ageint = 1)
Women.Positive <- asdt.inj.Women.Positive.del$s_ex[1]-asdt.inj.Women.Positive.del$ex[1] 
Women.Positive
#Unknown
asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 1)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 1)
asdt.inj.Women.Unknown.del <- asdt(allcause = asmr.all.Women.Unknown, i_cause = asmr.inj.Women.Unknown, ageint = 1)
Women.Unknown <- asdt.inj.Women.Unknown.del$s_ex[1]-asdt.inj.Women.Unknown.del$ex[1] 
Women.Unknown #0.684899



