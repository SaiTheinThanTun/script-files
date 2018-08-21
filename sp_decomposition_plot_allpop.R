#decomposing between periods all population####
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


#2007-2010####
dat.original <- dat
dat <- dat.original[dat.original$period.2011_15==FALSE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2007 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
#Positive
asmr.all.Women.Positive.2007 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
#Unknown
asmr.all.Women.Unknown.2007 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2007 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
#Positive
asmr.all.Men.Positive.2007 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
#Unknown
asmr.all.Men.Unknown.2007 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)

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

#between HIV status within 2007
# #Women
# decomList.Women.Positive.Negative <- decomList(allcause.A = asmr.all.Women.Positive, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Positive, i_cause.B = injEventList.Women.Negative, ageint = 15)
# colnames(decomList.Women.Positive.Negative) <- c(colnames(decomList.Women.Positive.Negative)[1:7],injNames)
# #decomList.Women.Positive.Negative
# if(creation) write.csv(decomList.Women.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_Negative_2007.csv",sep = ""))
# 
# decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 15)
# colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative_2007.csv",sep = ""))
# 
# #Men
# decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 15)
# colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
# #decomList.Men.Positive.Negative
# if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative_2007.csv",sep = ""))
# 
# decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 15)
# colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative_2007.csv",sep = ""))

#2011-2015####
dat <- dat.original[dat.original$period.2011_15==TRUE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2011 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
#Positive
asmr.all.Women.Positive.2011 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
#Unknown
asmr.all.Women.Unknown.2011 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2011 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
#Positive
asmr.all.Men.Positive.2011 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)

#Unknown
asmr.all.Men.Unknown.2011 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)


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
#POSITIVES
decomList.Women.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Women.Positive.2007, allcause.B = asmr.all.Women.Positive.2011, i_cause.A=injEventList.Women.Positive.2007, i_cause.B = injEventList.Women.Positive.2011, ageint = 15)
colnames(decomList.Women.Positive.2007.2011) <- c(colnames(decomList.Women.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_2011_2007.csv",sep = ""))

decomList.Men.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Men.Positive.2007, allcause.B = asmr.all.Men.Positive.2011, i_cause.A=injEventList.Men.Positive.2007, i_cause.B = injEventList.Men.Positive.2011, ageint = 15)
colnames(decomList.Men.Positive.2007.2011) <- c(colnames(decomList.Men.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_2011_2007.csv",sep = ""))

#plot for comparison of positives between period###
#max(decomList.Women.Positive.2007.2011[,c(8:18)])
#max(decomList.Men.Positive.2007.2011[,c(8:18)])

Men.period.decom <- decomList.Men.Positive.2007.2011[,c(8:18)]*12

Women.period.decom <- decomList.Women.Positive.2007.2011[,c(8:18)]*12

period.decom <- rbind(Women.period.decom,Men.period.decom)
period.decom <- period.decom[,which(colSums(period.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("HIV+ Women","HIV+ Men"), each=4)
period.decom <- cbind(period.decom, age, sex)
#colnames(period.decom)[1] <- 'age'
period.decom.L <- reshape::melt(period.decom, id.vars=c("age", "sex"))

#labdat_period <- data.frame(x=4, y=4, lab=c("Total difference: 8.9 months", "Total difference: 3.2 months"), sex=c("HIV+ Men","HIV+ Women"))
labdat_period <- data.frame(x=4, y=4, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.Positive.2007.2011[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.Positive.2007.2011[,8:18]),2),"months")), sex=c("HIV+ Men","HIV+ Women"))
labdat_period2 <- data.frame(x=4, y=4, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.Positive.2007.2011[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.Positive.2007.2011[,7]),2),"years")), sex=c("HIV+ Men","HIV+ Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_positive_period_decom.png",sep = ""), width = 1100, height = 800)
ggplot(period.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = 'bottom')+
  geom_text(aes(x,y,label=lab), data=labdat_period, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_period2, size = 5, vjust=4)
dev.off()