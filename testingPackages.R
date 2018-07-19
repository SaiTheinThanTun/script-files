#survival exercise
#https://bioconnector.org/r-survival.html
library(dplyr)
library(survival)
library(survminer)

lung <- as_tibble(lung)

s <- Surv(lung$time, lung$status)
class(s)
s
head(lung)

sfit <- survfit(Surv(time, status)~1, data=lung)
sfit
summary(sfit)

sfit <- survfit(Surv(time, status)~sex, data=lung)
sfit
summary(sfit)

summary(sfit, times=seq(0, 1000, 100))

ggsurvplot(sfit)


#another package called Epi and popEpi####