#sp_lifetable, started 20180729
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')

#calander year 2007 until "2015-11-26", it would be a period LT with a long period####
#all cause
asmr.all <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ age, data=dat, scale = 1)
#summary(asmr.all, rate = T, ci.r = T)
#lt(asmr.all, per=10000, ci95=T)
lt.all <- lt(asmr.all)
write.csv(lt.all,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_all.csv",sep = ""))

#ext.inj cause
asmr.ext <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat, scale = 1)
lt.ext <- lt(asmr.ext)
write.csv(lt.ext,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_ext.csv",sep = ""))
