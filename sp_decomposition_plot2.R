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

#preparing lifetables for decomposing####
#Women
dat.Women <- dat[dat$sex=='Women',]
#Negative
asmr.all.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Women.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Negative',], scale = 1), ageint = 15)
#Positive
asmr.all.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Women.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Positive',], scale = 1), ageint = 15)
#Unknown
asmr.all.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Women.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Women[dat.Women$allFixed=='Unknown',], scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#Negative
asmr.all.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)
asmr.inj.Men.Negative <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Negative',], scale = 1), ageint = 15)

#Positive
asmr.all.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)
asmr.inj.Men.Positive <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Positive',], scale = 1), ageint = 15)

#Unknown
asmr.all.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)
asmr.inj.Men.Unknown <- lt(pyears(Surv(time=time0, time2 = timex, event = fail2) ~ agegrp15, data=dat.Men[dat.Men$allFixed=='Unknown',], scale = 1), ageint = 15)


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

# decomList.Women.Unknown.Negative <- decomList(allcause.A = asmr.all.Women.Unknown, allcause.B = asmr.all.Women.Negative, i_cause.A=injEventList.Women.Unknown, i_cause.B = injEventList.Women.Negative, ageint = 15)
# colnames(decomList.Women.Unknown.Negative) <- c(colnames(decomList.Women.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Women.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_Unknown_Negative.csv",sep = ""))

#Men
decomList.Men.Positive.Negative <- decomList(allcause.A = asmr.all.Men.Positive, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Positive, i_cause.B = injEventList.Men.Negative, ageint = 15)
colnames(decomList.Men.Positive.Negative) <- c(colnames(decomList.Men.Positive.Negative)[1:7],injNames)
#decomList.Men.Positive.Negative
if(creation) write.csv(decomList.Men.Positive.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Positive_Negative.csv",sep = ""))

# decomList.Men.Unknown.Negative <- decomList(allcause.A = asmr.all.Men.Unknown, allcause.B = asmr.all.Men.Negative, i_cause.A=injEventList.Men.Unknown, i_cause.B = injEventList.Men.Negative, ageint = 15)
# colnames(decomList.Men.Unknown.Negative) <- c(colnames(decomList.Men.Unknown.Negative)[1:7],injNames)
# if(creation) write.csv(decomList.Men.Unknown.Negative, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_Unknown_Negative.csv",sep = ""))

#plot for between HIV status comparison####
#decomList.Men.Positive.Negative #data to use
max(decomList.Men.Positive.Negative[,c(8:18)])
max(decomList.Women.Positive.Negative[,c(8:18)])

Men.HIV.decom <- decomList.Men.Positive.Negative[,c(8:18)]*12

Women.HIV.decom <- decomList.Women.Positive.Negative[,c(8:18)]*12

HIV.decom <- rbind(Women.HIV.decom,Men.HIV.decom)
HIV.decom <- HIV.decom[,which(colSums(HIV.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("Women","Men"), each=4)
HIV.decom <- cbind(HIV.decom, age, sex)
#colnames(HIV.decom)[1] <- 'age'
HIV.decom.L <- reshape::melt(HIV.decom, id.vars=c("age", "sex"))

labdat_hiv <- data.frame(x=4, y=5, lab=c("Total difference: 9.3 months", "Total difference: 1.3 months"), sex=c("Men","Women"))
png(paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_HIV_decom.png",sep = ""), width = 1100, height = 800)
ggplot(HIV.decom.L) +
  geom_bar(aes(x=age, y=value, fill=variable), stat="identity") + facet_grid(. ~ sex)+ coord_flip()+
  ylab("Life-year difference in months")+
  scale_fill_discrete(name = "External injury type")+
  theme(text = element_text(size=16), legend.position = "bottom")+
  geom_text(aes(x,y,label=lab), data=labdat_hiv, size = 6, vjust=1)
dev.off()
