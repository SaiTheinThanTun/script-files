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
lt.all <- lt(asmr.all, ageint = 1)
write.csv(lt.all,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_all.csv",sep = ""))

#ext.inj cause
asmr.ext <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ age, data=dat, scale = 1)
lt.ext <- lt(asmr.ext, ageint = 1)
write.csv(lt.ext,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_lt_ext.csv",sep = ""))

#all cause by 5 year age group
asmr.all.5 <- pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat, scale = 1)
lt(asmr.all.5, ageint = 5)

asdt <- function(allcause, i_cause, ageint, deletion=TRUE){
  #Associated single decrement life table
  #allcause and i_cause are dataframe resulted from `lt` function
  #ageint defines age interval
  #deletion specifies if othercause (allcause `minus` i_cause) is used
  tableNames <- c('Person-years','Event', 'Rate', 'nqx', 'npx', 'lx', 'ndx', 'nLx', 'Tx', 'ex')
  
  #check if the dateset input are correct
  if(sum(colnames(allcause) %in% tableNames)!=length(tableNames)) stop("allcause dataset incorrect")
  if(sum(colnames(i_cause) %in% tableNames)!=length(tableNames)) stop("i_cause dataset incorrect")
  
  #subsetting only necessary sections
  allcause <- allcause[,which(colnames(allcause %in% tableNames))]
  i_cause <- i_cause[,which(colnames(i_cause %in% tableNames))]
  
  if(deletion){
    new_nmx <- allcause$rate - i_cause$rate
  }
  else {new_nmx <- i_cause$rate}
  
  s_npx <- exp(-ageint*new_nmx) #s denotes star or asterisk*
  s_npx[length(allcause$event)] <- 0
  
  s_lx <- NA
  for(i in 1:length(allcause$event)){
    if(i==1){s_lx[i] <- 1} #radix of 1 so that ex can be directly got from Tx
    else {s_lx[i] <- s_lx[i-1]*s_npx[i-1]}
  }
  
  s_ndx <- s_lx*(1-s_npx)
  
  s_nLx <- s_ndx/new_nmx
  
  
}
