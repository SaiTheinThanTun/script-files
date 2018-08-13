#without residence in the model+interaction with age
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed):strata(agegrp15)+factor(sex)+factor(period.2011_15), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph
#interaction, without residence in the model
a15spi.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(allFixed):factor(period.2011_15), data=dat) 
summary(a15spi.cox)
anova(a15sp.cox,a15spi.cox)

#without residence in the model+ period strata
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+strata(period.2011_15), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph
#interaction, without residence in the model
a15spi.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+factor(allFixed):factor(period.2011_15), data=dat) 
summary(a15spi.cox)
anova(a15sp.cox,a15spi.cox)

#without residence in the model+agegrp strata *
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+strata(agegrp15)+factor(period.2011_15), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+strata(agegrp15)+ct, data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#time0
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+strata(agegrp15)+factor(period.2011_15)+factor(allFixed):time0, data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#timex
a15sp.cox <- coxph(Surv(time=timex, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat) #+factor(allFixed):timex
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#time
dat$time <- dat$timex-dat$time0
a15sp.cox <- coxph(Surv(time=time, event = fail2) ~ factor(allFixed)+factor(sex)+strata(agegrp15)+factor(period.2011_15)+factor(allFixed):time, data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph
survminer::ggcoxzph(a15sp.cox.ph)


#without residence in the model+residence strata ** better than before
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+strata(residence), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#without residence in the model+residence strata *** better than before
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+strata(agegrp15)+factor(period.2011_15)+strata(residence), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph


#without residence in the model+ residence2 strata
a15sp.cox <- coxph(Surv(time=timex-time0, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15)+strata(residence2), data=dat) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph


#testing tt
#without residence in the model+interaction with age
dat$time <- dat$timex-dat$time0
dat2 <- dat[which(colnames(dat) %in% c('time','fail2','allFixed','agegrp','sex','period.2011_15'))]
a15sp.cox <- coxph(Surv(time=time, event = fail2) ~ tt(allFixed)+factor(agegrp)+factor(sex)+factor(period.2011_15), data=dat2,tt=function(x,t,...) (x=='Positive')*log(t)) 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#coxphw
library(coxphw)

a15sp.cox <- coxphw(Surv(time=time, event = fail2) ~ factor(allFixed)+factor(sex)+factor(agegrp15)+factor(period.2011_15), data=dat, template="AHR") 
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#timesplit, category
vet2 <- survSplit(Surv(time=time, event = fail2) ~ ., data= dat, cut=c(.5), episode= "tgroup", id="id2")
a15sp.cox <- coxph(Surv(time=time, event = fail2) ~ factor(allFixed):strata(tgroup)+factor(sex)+factor(period.2011_15), data=vet2) #factor(agegrp15)
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#timesplit, continuous
a15sp.cox <- coxph(Surv(time=time, event = fail2) ~ factor(allFixed)+factor(allFixed):age+factor(sex)+factor(period.2011_15), data=dat) #factor(agegrp15)
summary(a15sp.cox)
a15sp.cox.ph <- cox.zph(a15sp.cox)
a15sp.cox.ph

#
tmp <- coxph(Surv(time = time0,time2 = timex,event = fail2) ~ factor(allFixed)+factor(sex)+factor(period.2011_15)+factor(agegrp15), data=dat)
summary(tmp)
tmp.ph <- cox.zph(tmp)
tmp.ph

(xtabs( ~ fail2+agegrp15, data=dat))
