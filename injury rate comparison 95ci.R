#age-standardized injury mortality per calendar year, version 2####
#this version include 95%CI
library(survival)
library(ggplot2)
library(data.table)
library(reshape)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- FALSE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')


WSP <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/WSP.csv")
#Men
std_inj_ct_Men <- std_inj_ct_Men_low <- std_inj_ct_Men_hi <- NA
dat_Men <- dat[dat$sex=='Men',]
for(i in 1:length(names(table(dat$ct)))){
  tmp <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat_Men[dat_Men$ct==names(table(dat_Men$ct))[i],], scale = 1)
  tmp <- pyears2(tmp, per=10000)
  std_inj_ct_Men[[i]] <- sum(tmp$Rate*WSP, na.rm = T)
  std_inj_ct_Men_low[[i]] <- sum(tmp[,4]*WSP, na.rm = T)
  std_inj_ct_Men_hi[[i]] <- sum(tmp[,5]*WSP, na.rm = T)
}
#Women
std_inj_ct_Women <- std_inj_ct_Women_low <- std_inj_ct_Women_hi <- NA
dat_Women <- dat[dat$sex=='Women',]
for(i in 1:length(names(table(dat$ct)))){
  tmp <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat_Women[dat_Women$ct==names(table(dat_Women$ct))[i],], scale = 1)
  tmp <- pyears2(tmp, per=10000)
  std_inj_ct_Women[[i]] <- sum(tmp$Rate*WSP, na.rm = T)
  std_inj_ct_Women_low[[i]] <- sum(tmp[,4]*WSP, na.rm = T)
  std_inj_ct_Women_hi[[i]] <- sum(tmp[,5]*WSP, na.rm = T)
}

uMkhanyakude <- as.data.frame(rbind(std_inj_ct_Men,std_inj_ct_Men_low, std_inj_ct_Men_hi, std_inj_ct_Women, std_inj_ct_Women_low, std_inj_ct_Women_hi))

uMkhanyakude <- cbind(cbind(c(rep('Male',3), rep('Female',3)),rep(c('estimate','low','high'),2),rep('uMkhanyakude',2*3)), uMkhanyakude)
colnames(uMkhanyakude) <- c('sex_name','metric','location_name', as.numeric(names(table(dat$ct))))
uMkhanyakudeMelt <- reshape::melt(uMkhanyakude, id.vars=c('sex_name','metric','location_name'))
colnames(uMkhanyakudeMelt) <- c('sex_name','metric','location_name','year','val')
umk_rate <- dcast(uMkhanyakudeMelt, sex_name+location_name+year ~ metric)
colnames(umk_rate) <- c('sex_name','location_name','year','val', 'upper','lower')
umk_rate$year <- as.numeric(levels(umk_rate$year))[umk_rate$year]
#if(creation) write.csv(uMkhanyakudeMelt,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_std_rate_peryear.csv",sep = ""))
# 
# umk_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/_std_rate_peryear.csv")
# umk_rate <- umk_rate[,colnames(umk_rate) %in% c('sex_name','location_name','year','val')]

SA_rate <- read.csv("/Users/sai/OneDrive/Summer Project/data/GBD/rate_injury_SA_2007-15.csv")
SA_rate <- SA_rate[,colnames(SA_rate) %in% c('sex_name','location_name','year','val', 'upper','lower')]
SA_rate <- SA_rate[SA_rate$sex_name!='Both',]
SA_rate$val <- SA_rate$val/10
SA_rate$upper <- SA_rate$upper/10
SA_rate$lower <- SA_rate$lower/10

both_rate <- rbind(umk_rate,SA_rate)
if(TRUE) both_rate <- both_rate[both_rate$year!=2015,]

#sort location name, year and sex
both_rate <- both_rate[order(both_rate$location_name, both_rate$sex_name,both_rate$year),]
row.names(both_rate) <- 1:nrow(both_rate)
names(both_rate)[1] <- 'sex'
names(both_rate)[2] <- 'location'

png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_std_inj_rate_peryear_compare.png",sep = ""), width = 800, height = 760)
ggplot(both_rate, aes(x=year, y=val, linetype=location, col=sex))+
  geom_line(size=1)+
  ggtitle("Standardized injury-related mortality rate in \n uMkhanyakude vs national average")+
  ylab("Age-standardized injury-related mortality rate per 10,000 PY")+
  theme(text = element_text(size=16), legend.position = 'bottom')+
  scale_linetype_manual(values=c("solid","dotdash"))+
  scale_size_manual(values=c(3, 2))+
  geom_ribbon(data=both_rate[both_rate$location=='uMkhanyakude',],aes(ymin=lower, ymax=upper), linetype='dotted', alpha=0.15)
dev.off()