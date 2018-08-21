#compare HIV+Known ART vs HIV+Unknown ART####
#HIV+Known ART - HIV+Unknown ART
#sp_lifetable, started 20180729
#all in agegrp = 15
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
library(ggplot2)
library(data.table)
library(reshape)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- FALSE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')
#write.csv(dat,'processed.csv')
injNames <- readRDS('accidentNames.RDS')

#putting new variables for each injury
for(i in 1: length(injNames)){
  x <- NA
  x <- (dat$fail2==1)*(dat$COD %in% injNames[i])
  dat <- cbind(dat,x)
  colnames(dat)[ncol(dat)] <- injNames[i]
}

#preparing lifetables for decomposing
#Women
dat.Women <- dat[dat$sex=='Women',]
#HIVposART
asmr.all.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 15)
asmr.inj.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 15)
#HIVposUnknownART
asmr.all.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#HIVposART
asmr.all.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 15)
asmr.inj.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 15)

#HIVposUnknownART
asmr.all.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)


#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns
#Women
datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Women.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposART"),]
injEventList.Women.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Men
datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Men.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposART"),]
injEventList.Men.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Women
decomList.Women.HIVposUnknownART.HIVposART <- decomList(allcause.A = asmr.all.Women.HIVposUnknownART, allcause.B = asmr.all.Women.HIVposART, i_cause.A=injEventList.Women.HIVposUnknownART, i_cause.B = injEventList.Women.HIVposART, ageint = 15)
colnames(decomList.Women.HIVposUnknownART.HIVposART) <- c(colnames(decomList.Women.HIVposUnknownART.HIVposART)[1:7],injNames)
#decomList.Women.HIVposUnknownART.HIVposART
if(creation) write.csv(decomList.Women.HIVposUnknownART.HIVposART, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_HIVposUnknownART_HIVposART.csv",sep = ""))

#Men
decomList.Men.HIVposUnknownART.HIVposART <- decomList(allcause.A = asmr.all.Men.HIVposUnknownART, allcause.B = asmr.all.Men.HIVposART, i_cause.A=injEventList.Men.HIVposUnknownART, i_cause.B = injEventList.Men.HIVposART, ageint = 15)
colnames(decomList.Men.HIVposUnknownART.HIVposART) <- c(colnames(decomList.Men.HIVposUnknownART.HIVposART)[1:7],injNames)
#decomList.Men.HIVposUnknownART.HIVposART
if(creation) write.csv(decomList.Men.HIVposUnknownART.HIVposART, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_HIVposUnknownART_HIVposART.csv",sep = ""))


#plot 
#decomList.Men.HIVposUnknownART.HIVposART #data to use
max(decomList.Men.HIVposUnknownART.HIVposART[,c(8:18)])
max(decomList.Women.HIVposUnknownART.HIVposART[,c(8:18)])

Men.HIV.decom <- decomList.Men.HIVposUnknownART.HIVposART[,c(8:18)]*12

Women.HIV.decom <- decomList.Women.HIVposUnknownART.HIVposART[,c(8:18)]*12

HIV.decom <- rbind(Women.HIV.decom,Men.HIV.decom)
HIV.decom <- HIV.decom[,which(colSums(HIV.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("Women","Men"), each=4)
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))

labdat_hiv <- data.frame(x=4, y=2.5, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.HIVposUnknownART.HIVposART[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.HIVposUnknownART.HIVposART[,8:18]),2),"months")), sex=c("Men","Women"))
labdat_hiv2 <- data.frame(x=4, y=2.5, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.HIVposUnknownART.HIVposART[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.HIVposUnknownART.HIVposART[,7]),2),"years")), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_ART_decom.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_hiv2, size = 5, vjust=4)
dev.off()


#counter plot of Known vs Unknown ART####
#sp_lifetable, started 20180729
#all in agegrp = 15
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
library(ggplot2)
library(data.table)
library(reshape)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- FALSE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')
#write.csv(dat,'processed.csv')
injNames <- readRDS('accidentNames.RDS')

#putting new variables for each injury
for(i in 1: length(injNames)){
  x <- NA
  x <- (dat$fail2==1)*(dat$COD %in% injNames[i])
  dat <- cbind(dat,x)
  colnames(dat)[ncol(dat)] <- injNames[i]
}

#preparing lifetables for decomposing
#Women
dat.Women <- dat[dat$sex=='Women',]
#HIVposUnknownART
asmr.all.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
#HIVposART
asmr.all.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 15)
asmr.inj.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#HIVposUnknownART
asmr.all.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)

#HIVposART
asmr.all.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 15)
asmr.inj.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 15)


#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns
#Women
datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposART"),]
injEventList.Women.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Women.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Men
datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposART"),]
injEventList.Men.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Men.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Women
decomList.Women.HIVposART.HIVposUnknownART <- decomList(allcause.A = asmr.all.Women.HIVposART, allcause.B = asmr.all.Women.HIVposUnknownART, i_cause.A=injEventList.Women.HIVposART, i_cause.B = injEventList.Women.HIVposUnknownART, ageint = 15)
colnames(decomList.Women.HIVposART.HIVposUnknownART) <- c(colnames(decomList.Women.HIVposART.HIVposUnknownART)[1:7],injNames)
#decomList.Women.HIVposART.HIVposUnknownART
if(creation) write.csv(decomList.Women.HIVposART.HIVposUnknownART, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_HIVposART_HIVposUnknownART.csv",sep = ""))

#Men
decomList.Men.HIVposART.HIVposUnknownART <- decomList(allcause.A = asmr.all.Men.HIVposART, allcause.B = asmr.all.Men.HIVposUnknownART, i_cause.A=injEventList.Men.HIVposART, i_cause.B = injEventList.Men.HIVposUnknownART, ageint = 15)
colnames(decomList.Men.HIVposART.HIVposUnknownART) <- c(colnames(decomList.Men.HIVposART.HIVposUnknownART)[1:7],injNames)
#decomList.Men.HIVposART.HIVposUnknownART
if(creation) write.csv(decomList.Men.HIVposART.HIVposUnknownART, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_HIVposART_HIVposUnknownART.csv",sep = ""))


#plot
#decomList.Men.HIVposART.HIVposUnknownART #data to use
max(decomList.Men.HIVposART.HIVposUnknownART[,c(8:18)])
max(decomList.Women.HIVposART.HIVposUnknownART[,c(8:18)])

Men.HIV.decom <- decomList.Men.HIVposART.HIVposUnknownART[,c(8:18)]*12

Women.HIV.decom <- decomList.Women.HIVposART.HIVposUnknownART[,c(8:18)]*12

HIV.decom <- rbind(Women.HIV.decom,Men.HIV.decom)
HIV.decom <- HIV.decom[,which(colSums(HIV.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("Women","Men"), each=4)
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))

labdat_hiv <- data.frame(x=4, y=5, lab=c("Total difference: 9.3 months", "Total difference: 1.3 months"), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_ART_decom_opposite.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 6, vjust=1)
dev.off()


#HIV negative vs HIV+ ART unknown####
#all in agegrp = 15
library(survival)
library(survminer)
library(Epi)
library(popEpi)
library(Hmisc)
library(ggplot2)
library(data.table)
library(reshape)
source('~/OneDrive/Summer Project/script-files/sp_functions.R')
creation <- TRUE #WRITE NEW FILES

setwd('~/OneDrive/Summer Project/data')
dat <- readRDS('dat_lifetable.RDS')
#write.csv(dat,'processed.csv')
injNames <- readRDS('accidentNames.RDS')

#putting new variables for each injury
for(i in 1: length(injNames)){
  x <- NA
  x <- (dat$fail2==1)*(dat$COD %in% injNames[i])
  dat <- cbind(dat,x)
  colnames(dat)[ncol(dat)] <- injNames[i]
}

#preparing lifetables for decomposing
#Women
dat.Women <- dat[dat$sex=='Women',]
#HIV_Negative
asmr.all.Women.HIV_Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIV Negative',], scale = 1), ageint = 15)
asmr.inj.Women.HIV_Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIV Negative',], scale = 1), ageint = 15)
#HIVposUnknownART
asmr.all.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#HIV_Negative
asmr.all.Men.HIV_Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIV Negative',], scale = 1), ageint = 15)
asmr.inj.Men.HIV_Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIV Negative',], scale = 1), ageint = 15)

#HIVposUnknownART
asmr.all.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)
asmr.inj.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 15)


#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns
#Women
datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Women.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$art_status2=='HIV Negative'),]
injEventList.Women.HIV_Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Men
datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Men.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$art_status2=='HIV Negative'),]
injEventList.Men.HIV_Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Women
decomList.Women.HIVposUnknownART.HIV_Negative <- decomList(allcause.A = asmr.all.Women.HIVposUnknownART, allcause.B = asmr.all.Women.HIV_Negative, i_cause.A=injEventList.Women.HIVposUnknownART, i_cause.B = injEventList.Women.HIV_Negative, ageint = 15)
colnames(decomList.Women.HIVposUnknownART.HIV_Negative) <- c(colnames(decomList.Women.HIVposUnknownART.HIV_Negative)[1:7],injNames)
#decomList.Women.HIVposUnknownART.HIV_Negative
if(creation) write.csv(decomList.Women.HIVposUnknownART.HIV_Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_HIVposUnknownART_HIV_Negative.csv",sep = ""))

#Men
decomList.Men.HIVposUnknownART.HIV_Negative <- decomList(allcause.A = asmr.all.Men.HIVposUnknownART, allcause.B = asmr.all.Men.HIV_Negative, i_cause.A=injEventList.Men.HIVposUnknownART, i_cause.B = injEventList.Men.HIV_Negative, ageint = 15)
colnames(decomList.Men.HIVposUnknownART.HIV_Negative) <- c(colnames(decomList.Men.HIVposUnknownART.HIV_Negative)[1:7],injNames)
#decomList.Men.HIVposUnknownART.HIV_Negative
if(creation) write.csv(decomList.Men.HIVposUnknownART.HIV_Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_HIVposUnknownART_HIV_Negative.csv",sep = ""))


#plot 
#decomList.Men.HIVposUnknownART.HIV_Negative #data to use
max(decomList.Men.HIVposUnknownART.HIV_Negative[,c(8:18)])
max(decomList.Women.HIVposUnknownART.HIV_Negative[,c(8:18)])

Men.HIV.decom <- decomList.Men.HIVposUnknownART.HIV_Negative[,c(8:18)]*12

Women.HIV.decom <- decomList.Women.HIVposUnknownART.HIV_Negative[,c(8:18)]*12

HIV.decom <- rbind(Women.HIV.decom,Men.HIV.decom)
HIV.decom <- HIV.decom[,which(colSums(HIV.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("Women","Men"), each=4)
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))

labdat_hiv <- data.frame(x=4, y=2.5, lab=c(paste("Total difference:",round(12*sum(decomList.Men.HIVposUnknownART.HIV_Negative[,8:18]),2),"months"), paste("Total difference:",round(12*sum(decomList.Women.HIVposUnknownART.HIV_Negative[,8:18]),2),"months")), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_NegVsposUnkART_decom.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 6, vjust=1)
dev.off()
