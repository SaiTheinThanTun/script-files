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
asmr.all.Women.all.2007 <- asmr.all.Women.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women, scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#all
asmr.all.Men.all.2007 <- asmr.all.Men.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men, scale = 1), ageint = 15)

#2007-2010#
#Women

datL <- dat[which(dat$sex=="Women"),]
injEventList.Women.all.2007 <- injEventList.Women.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Men

datL <- dat[which(dat$sex=="Men"),]
injEventList.Men.all.2007 <- injEventList.Men.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#2011-2015####
dat <- dat.original[dat.original$period.2011_15==TRUE,]

#Women
dat.Women <- dat[dat$sex=='Women',]
#all
asmr.all.Women.all.2011 <- asmr.all.Women.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Women, scale = 1), ageint = 15)

#Men
dat.Men <- dat[dat$sex=='Men',]
#all
asmr.all.Men.all.2011 <- asmr.all.Men.all <- lt(pyears(Surv(time=time0, time2 = timex, event = fail0) ~ agegrp15, data=dat.Men, scale = 1), ageint = 15)


#2011-2015#
#Women

datL <- dat[which(dat$sex=="Women"),]
injEventList.Women.all.2011 <- injEventList.Women.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})


#Men

datL <- dat[which(dat$sex=="Men"),]
injEventList.Men.all.2011 <- injEventList.Men.all <- lapply(injNames,function(x){
  pyears(Surv(time=time0, time2 = timex, event = datL[,which(colnames(datL)==x)]) ~ agegrp15, data=datL, scale = 1)$event
})

#all,
decomList.Women.all.2007.2011 <- decomList(allcause.A = asmr.all.Women.all.2007, allcause.B = asmr.all.Women.all.2011, i_cause.A=injEventList.Women.all.2007, i_cause.B = injEventList.Women.all.2011, ageint = 15)
colnames(decomList.Women.all.2007.2011) <- c(colnames(decomList.Women.all.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Women.all.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Women_all_2011_2007.csv",sep = ""))

decomList.Men.all.2007.2011 <- decomList(allcause.A = asmr.all.Men.all.2007, allcause.B = asmr.all.Men.all.2011, i_cause.A=injEventList.Men.all.2007, i_cause.B = injEventList.Men.all.2011, ageint = 15)
colnames(decomList.Men.all.2007.2011) <- c(colnames(decomList.Men.all.2007.2011)[1:7],injNames)
if(creation) write.csv(decomList.Men.all.2007.2011, paste("~/OneDrive/Summer Project/output/",gsub("\\:","",Sys.time()),"_decomList_Men_all_2011_2007.csv",sep = ""))

#plot for comparison of alls between period###
#max(decomList.Women.all.2007.2011[,c(8:18)])
#max(decomList.Men.all.2007.2011[,c(8:18)])

Men.period.decom <- decomList.Men.all.2007.2011[,c(8:18)]*12

Women.period.decom <- decomList.Women.all.2007.2011[,c(8:18)]*12

period.decom <- rbind(Women.period.decom,Men.period.decom)
period.decom <- period.decom[,which(colSums(period.decom)!=0)]
age <- rep(c("15-29","30-44","45-59","60+"),2)
sex <- rep(c("All Women","All Men"), each=4)
period.decom <- cbind(period.decom, age, sex)
#colnames(period.decom)[1] <- 'age'
period.decom.L <- reshape::melt(period.decom, id.vars=c("age", "sex"))

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
