#decomposition plot all 4, 20180829
#decomposition in 5year age group
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
#putting new variables for each injury
for(i in 1: length(injNames)){
  x <- NA
  x <- (dat$fail2==1)*(dat$COD %in% injNames[i])
  dat <- cbind(dat,x)
  colnames(dat)[ncol(dat)] <- injNames[i]
}

#2007-2010####
dat.original <- dat
dat <- dat.original[dat.original$period.2011_15==FALSE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#all
asmr.all.Women.all.2007 <- asmr.all.Women.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women, scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#all
asmr.all.Men.all.2007 <- asmr.all.Men.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men, scale = 1), ageint = 5)

#2007-2010#
#Women

datL <- dat[which(dat$sex=="Women"),]
injEventList.Women.all.2007 <- injEventList.Women.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})


#Men

datL <- dat[which(dat$sex=="Men"),]
injEventList.Men.all.2007 <- injEventList.Men.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#2011-2015####
dat <- dat.original[dat.original$period.2011_15==TRUE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#all
asmr.all.Women.all.2011 <- asmr.all.Women.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women, scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#all
asmr.all.Men.all.2011 <- asmr.all.Men.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men, scale = 1), ageint = 5)


#2011-2015#
#Women

datL <- dat[which(dat$sex=="Women"),]
injEventList.Women.all.2011 <- injEventList.Women.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})


#Men

datL <- dat[which(dat$sex=="Men"),]
injEventList.Men.all.2011 <- injEventList.Men.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#all,
decomList.Women.all.2007.2011 <- decomList(allcause.A = asmr.all.Women.all.2007, allcause.B = asmr.all.Women.all.2011, i_cause.A=injEventList.Women.all.2007, i_cause.B = injEventList.Women.all.2011, ageint = 5)
colnames(decomList.Women.all.2007.2011) <- c(colnames(decomList.Women.all.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.all.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_all_2011_2007.csv",sep = ""))

decomList.Men.all.2007.2011 <- decomList(allcause.A = asmr.all.Men.all.2007, allcause.B = asmr.all.Men.all.2011, i_cause.A=injEventList.Men.all.2007, i_cause.B = injEventList.Men.all.2011, ageint = 5)
colnames(decomList.Men.all.2007.2011) <- c(colnames(decomList.Men.all.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.all.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_all_2011_2007.csv",sep = ""))

#plot for comparison of alls between period
#max(decomList.Women.all.2007.2011[,c(8:18)])
#max(decomList.Men.all.2007.2011[,c(8:18)])

Men.period.decom <- decomList.Men.all.2007.2011[,c(8:18)]*12



Women.period.decom <- decomList.Women.all.2007.2011[,c(8:18)]*12



period.decom <- rbind(Women.period.decom,Men.period.decom)
period.decom <- period.decom[,which(colSums(period.decom)!=0)]
age <- rep(names(table(dat$agegrp)),2)
sex <- rep(c("All Women","All Men"), each=length(table(dat$agegrp)))
period.decom <- cbind(period.decom, age, sex)
#colnames(period.decom)[1] <- 'age'
period.decom.L <- reshape::melt(period.decom, id.vars=c("age", "sex"))
summ.period.decom.L <- dcast(period.decom.L, variable ~ sex, sum)
write.csv(period.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_period_decom_all.csv",sep = "") )
write.csv(summ.period.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_period_decom_all_summary.csv",sep = "") )

#labdat_period <- data.frame(x=4, y=4, lab=c("Total difference: 8.9 months", "Total difference: 3.2 months"), sex=c("HIV+ Men","HIV+ Women"))
labdat_period <- data.frame(x=3, y=2, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.all.2007.2011[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.all.2007.2011[,8:18]),2),"months")), sex=c("All Men","All Women"))
labdat_period2 <- data.frame(x=3, y=2, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.all.2007.2011[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.all.2007.2011[,7]),2),"years")), sex=c("All Men","All Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_all_period_decom.png",sep = ""), width = 1100, height = 800)
ggplot(period.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = 'bottom')+
  geom_text(aes(x,y,label=lab), data=labdat_period, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_period2, size = 5, vjust=4)
dev.off()

dat <- dat.original #restroing original data for further analysis


#decomposing between HIV status####
#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
#Positive
asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
#Unknown
asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)

#Positive
asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)

#Unknown
asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)


#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#Women
decomList.Women.Positive.Negative <- decomList(allcause.A = asmr.all.Women.Positive, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Positive, i_cause.B = injEventList.Women.Negative, ageint = 5)
colnames(decomList.Women.Positive.Negative) <- c(colnames(decomList.Women.Positive.Negative)[1:7],injNames)
#decomList.Women.Positive.Negative
if(creation) write.csv(decomList.Women.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_Negative.csv",sep = ""))

# decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 5)
# colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative.csv",sep = ""))

#Men
decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 5)
colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
#decomList.Men.Positive.Negative
if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative.csv",sep = ""))

# decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 5)
# colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative.csv",sep = ""))

#plot for between HIV status comparison
#decomList.Men.Positive.Negative #data to use
max(decomList.Men.Positive.Negative[,c(8:18)])
max(decomList.Women.Positive.Negative[,c(8:18)])

Men.HIV.decom <- decomList.Men.Positive.Negative[,c(8:18)]*12



Women.HIV.decom <- decomList.Women.Positive.Negative[,c(8:18)]*12


HIV.decom <- rbind(Women.HIV.decom,Men.HIV.decom)
HIV.decom <- HIV.decom[,which(colSums(HIV.decom)!=0)]
age <- rep(names(table(dat$agegrp)),2)
sex <- rep(c("Women","Men"), each=length(table(dat$agegrp)))
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))
write.csv(HIV.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_hiv_decom.csv",sep = "") )
summ.HIV.decom.L <- dcast(HIV.decom.L, variable ~ sex, sum)
write.csv(summ.HIV.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIV_decom_summary.csv",sep = "") )

#labdat_hiv <- data.frame(x=4, y=5, lab=c("Total difference: 9.3 months", "Total difference: 1.3 months"), sex=c("Men","Women"))
labdat_hiv <- data.frame(x=6, y=2, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.Positive.Negative[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.Positive.Negative[,8:18]),2),"months")), sex=c("Men","Women"))
labdat_hiv2 <- data.frame(x=6, y=2, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.Positive.Negative[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.Positive.Negative[,7]),2),"years")), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIV_decom.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_hiv2, size = 5, vjust=4)
dev.off()

#compare HIV+Known ART vs HIV+Unknown ART####
#preparing lifetables for decomposing
#Women
dat.Women <- dat[dat$sex=='Women',]
#HIVposART
asmr.all.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 5)
asmr.inj.Women.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$art_status2=='HIVposART',], scale = 1), ageint = 5)
#HIVposUnknownART
asmr.all.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 5)
asmr.inj.Women.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$art_status2=='HIVposUnknownART',], scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#HIVposART
asmr.all.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 5)
asmr.inj.Men.HIVposART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$art_status2=='HIVposART',], scale = 1), ageint = 5)

#HIVposUnknownART
asmr.all.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 5)
asmr.inj.Men.HIVposUnknownART <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$art_status2=='HIVposUnknownART',], scale = 1), ageint = 5)


#DECOMPOSITION BY AGE, SEX, HIV STATUS AND EACH INJURY COD, ie many columns
#Women
datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Women.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$art_status2=="HIVposART"),]
injEventList.Women.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})


#Men
datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposUnknownART"),]
injEventList.Men.HIVposUnknownART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$art_status2=="HIVposART"),]
injEventList.Men.HIVposART <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})


#Women
decomList.Women.HIVposUnknownART.HIVposART <- decomList(allcause.A = asmr.all.Women.HIVposUnknownART, allcause.B = asmr.all.Women.HIVposART, i_cause.A=injEventList.Women.HIVposUnknownART, i_cause.B = injEventList.Women.HIVposART, ageint = 5)
colnames(decomList.Women.HIVposUnknownART.HIVposART) <- c(colnames(decomList.Women.HIVposUnknownART.HIVposART)[1:7],injNames)
#decomList.Women.HIVposUnknownART.HIVposART
if(creation) write.csv(decomList.Women.HIVposUnknownART.HIVposART, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_HIVposUnknownART_HIVposART.csv",sep = ""))

#Men
decomList.Men.HIVposUnknownART.HIVposART <- decomList(allcause.A = asmr.all.Men.HIVposUnknownART, allcause.B = asmr.all.Men.HIVposART, i_cause.A=injEventList.Men.HIVposUnknownART, i_cause.B = injEventList.Men.HIVposART, ageint = 5)
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
age <- rep(names(table(dat$agegrp)),2)
sex <- rep(c("Women","Men"), each=length(table(dat$agegrp)))
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))
write.csv(HIV.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_art_decom.csv",sep = "") )

summ.HIV.decom.L <- dcast(HIV.decom.L, variable ~ sex, sum)
write.csv(summ.HIV.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_ART_decom_summary.csv",sep = "") )

labdat_hiv <- data.frame(x=8, y=2, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.HIVposUnknownART.HIVposART[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.HIVposUnknownART.HIVposART[,8:18]),2),"months")), sex=c("Men","Women"))
labdat_hiv2 <- data.frame(x=8, y=2, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.HIVposUnknownART.HIVposART[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.HIVposUnknownART.HIVposART[,7]),2),"years")), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_ART_decom.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_hiv2, size = 5, vjust=4)
dev.off()



#decomposing between periods for HIV positive####
#2007-2010####
dat <- dat.original[dat.original$period.2011_15==FALSE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2007 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
#Positive
asmr.all.Women.Positive.2007 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
#Unknown
asmr.all.Women.Unknown.2007 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2007 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
#Positive
asmr.all.Men.Positive.2007 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
#Unknown
asmr.all.Men.Unknown.2007 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)

#2007-2010#
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive.2007 <- injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative.2007 <- injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown.2007 <- injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive.2007 <- injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative.2007 <- injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown.2007 <- injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#2011-2015####
dat <- dat.original[dat.original$period.2011_15==TRUE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative.2011 <- asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 5)
#Positive
asmr.all.Women.Positive.2011 <- asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 5)
#Unknown
asmr.all.Women.Unknown.2011 <- asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 5)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative.2011 <- asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 5)
#Positive
asmr.all.Men.Positive.2011 <- asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 5)

#Unknown
asmr.all.Men.Unknown.2011 <- asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 5)


#2011-2015#
#Women
datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Positive"),]
injEventList.Women.Positive.2011 <- injEventList.Women.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Negative"),]
injEventList.Women.Negative.2011 <- injEventList.Women.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Women" & dat$allFixed=="Unknown"),]
injEventList.Women.Unknown.2011 <- injEventList.Women.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

#Men
datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Positive"),]
injEventList.Men.Positive.2011 <- injEventList.Men.Positive <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Negative"),]
injEventList.Men.Negative.2011 <- injEventList.Men.Negative <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})

datL <- dat[which(dat$sex=="Men" & dat$allFixed=="Unknown"),]
injEventList.Men.Unknown.2011 <- injEventList.Men.Unknown <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp, data=datL, scale = 1)$event
})
#POSITIVES
decomList.Women.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Women.Positive.2007, allcause.B = asmr.all.Women.Positive.2011, i_cause.A=injEventList.Women.Positive.2007, i_cause.B = injEventList.Women.Positive.2011, ageint = 5)
colnames(decomList.Women.Positive.2007.2011) <- c(colnames(decomList.Women.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Positive_2011_2007.csv",sep = ""))

decomList.Men.Positive.2007.2011 <- decomList(allcause.A = asmr.all.Men.Positive.2007, allcause.B = asmr.all.Men.Positive.2011, i_cause.A=injEventList.Men.Positive.2007, i_cause.B = injEventList.Men.Positive.2011, ageint = 5)
colnames(decomList.Men.Positive.2007.2011) <- c(colnames(decomList.Men.Positive.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.Positive.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_2011_2007.csv",sep = ""))


Men.period.decom <- decomList.Men.Positive.2007.2011[,c(8:18)]*12


Women.period.decom <- decomList.Women.Positive.2007.2011[,c(8:18)]*12


period.decom <- rbind(Women.period.decom,Men.period.decom)
period.decom <- period.decom[,which(colSums(period.decom)!=0)]
age <- rep(names(table(dat$agegrp)),2)
sex <- rep(c("HIV+ Women","HIV+ Men"), each=length(table(dat$agegrp)))
period.decom <- cbind(period.decom, age, sex)
#colnames(period.decom)[1] <- 'age'
period.decom.L <- reshape::melt(period.decom, id.vars=c("age", "sex"))
write.csv(period.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIVpos_period_decom.csv",sep = "") )

summ.period.decom.L <- dcast(period.decom.L, variable ~ sex, sum)
write.csv(summ.period.decom.L,paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_positive_period_decom_summary.csv",sep = "") )

#labdat_period <- data.frame(x=4, y=4, lab=c("Total difference: 8.9 months", "Total difference: 3.2 months"), sex=c("HIV+ Men","HIV+ Women"))
labdat_period <- data.frame(x=8, y=-2, lab=c(paste("Total difference, injury:",round(12*sum(decomList.Men.Positive.2007.2011[,8:18]),2),"months"), paste("Total difference, injury:",round(12*sum(decomList.Women.Positive.2007.2011[,8:18]),2),"months")), sex=c("HIV+ Men","HIV+ Women"))
labdat_period2 <- data.frame(x=8, y=-2, lab=c(paste("Total difference, all-cause:",round(sum(decomList.Men.Positive.2007.2011[,7]),2),"years"), paste("Total difference, all-cause:",round(sum(decomList.Women.Positive.2007.2011[,7]),2),"years")), sex=c("HIV+ Men","HIV+ Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_positive_period_decom.png",sep = ""), width = 1200, height = 800)
ggplot(period.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = 'bottom')+
  geom_text(aes(x,y,label=lab), data=labdat_period, size = 5, vjust=1)+
  geom_text(aes(x,y,label=lab), data=labdat_period2, size = 5, vjust=4)
dev.off()

dat <- dat.original #restroing original data for further analysis