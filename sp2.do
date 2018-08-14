clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
use "processed2.dta"
log using "sp2_log_cox.log", replace
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)
//stset exit, id(idno_original) failure(failure) time0(entry) origin(time dob) scale(365.25)

sts graph, by(allFixed) ylabel(minmax) //this is why sex has to be a strata, since it differentiates the trends

stcox i.allFixed

stcox i.allFixed i.sex //1.312

stcox i.allFixed i.sex i.agegrp15 //age is already there!

stcox i.residence //overall p=0.0566

stcox i.allFixed i.sex i.period_2011_15 //1.34, p=0.018 #best1#

stcox i.allFixed i.sex i.period_2011_15 i.residence //1.34, p=0.038, disregard residence

stcox i.allFixed i.sex i.period_2011_15 ct //1.36, p=0.012

stcox i.allFixed i.sex ct //1.38, p=0.009 #best2#
estat phtest, detail //global p=0

stcox i.allFixed ct, strata(sex) //1.46, p=0.002 #best3#
estat phtest, detail //global p=0.0048. HIV unknown, p=0.0034, is bringing down. it has time dependency

stcox i.allFixed ct, strata(sex) tvc(allFixed) texp(_t) nolog //positive: 1.887565   .2918111     4.11   0.000     1.394149    2.555612
//estat phtest, detail // can't do phtest after tvc

stcox i.allFixed i.sex ct, tvc(allFixed) texp(_t) nolog


//testing hazards proportionality
stphplot, by(allFixed) //all crossing each other


clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_cox_strata_interaction.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)

stcox i.allFixed i.sex ct, strata(agegrp15)
estat phtest, detail

stcox i.allFixed i.sex ct, strata(residence2)
estat phtest, detail

stcox i.allFixed ct, strata(sex residence2) //1.43, p=0.004, best4, because it doesn't violate PH assumption
estimates store A
estat phtest, detail

//stcox i.allFixed ct, strata(sex) //1.46, p=0.002 #best3#
//estimates store B
//lrtest A B //can't do lrtest

stcox i.allFixed##ct, strata(sex residence2) //best5##
estimates store B
lrtest A B 
estat phtest, detail

//tvc as category/factor
//stcox i.allFixed i.sex ct, tvc(i.allFixed) texp(_t) nolog //doesn't work with factor (i)

stcox i.allFixed i.sex ct i.residence2, tvc(allFixed) texp(_t) nolog

capture log close


clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_best.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)
stcox i.allFixed

stcox i.sex

stcox i.agegrp15

stcox i.period_2011_15

stcox ct

stcox i.residence //overall p=0.0566

stcox i.residence2

stcox i.allFixed i.sex //1.312

stcox i.allFixed i.sex i.period_2011_15 //1.34, p=0.018 #best1#

stcox i.allFixed i.sex ct //1.38, p=0.009 #best2#
estat phtest, detail //global p=0

stcox i.allFixed ct, strata(sex) //1.46, p=0.002 #best3#
estimates store A
estat phtest, detail //global p=0.0048. HIV unknown, p=0.0034, is bringing down. it has time dependency

//stcox i.allFixed ct, strata(sex residence2) //1.43, p=0.004, best4, because it doesn't violate PH assumption
//estimates store A
//estat phtest, detail


//stcox i.allFixed##ct, strata(sex residence2) //best5##
stcox i.allFixed##ct, strata(sex) //best5##
estimates store B
lrtest A B 
estat phtest, detail

log close

//onoffART
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_onoffART.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)
stcox ib2.art_status

stcox ib2.art_status i.sex //1.312

stcox ib2.art_status i.sex i.period_2011_15 //1.34, p=0.018 #best1#

stcox ib2.art_status i.sex ct //1.38, p=0.009 #best2#
estat phtest, detail //global p=0

stcox ib2.art_status ct, strata(sex) //1.46, p=0.002 #best3#
estat phtest, detail //global p=0.0048. HIV unknown, p=0.0034, is bringing down. it has time dependency

stcox ib2.art_status ct, strata(sex residence2) //1.43, p=0.004, best4, because it doesn't violate PH assumption
estimates store A
estat phtest, detail


stcox ib2.art_status##ct, strata(sex residence2) //best5##
estimates store B
lrtest A B 
estat phtest, detail

log close


//everART
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_everART.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)
stcox ib2.art_status2

stcox ib2.art_status2 i.sex //1.312

stcox ib2.art_status2 i.sex i.period_2011_15 //1.34, p=0.018 #best1#

stcox ib2.art_status2 i.sex ct //1.38, p=0.009 #best2#
estat phtest, detail //global p=0

stcox ib2.art_status2 ct, strata(sex) //1.46, p=0.002 #best3#
estat phtest, detail //global p=0.0048. HIV unknown, p=0.0034, is bringing down. it has time dependency

stcox ib2.art_status2 ct, strata(sex residence2) //1.43, p=0.004, best4, because it doesn't violate PH assumption
estimates store A
estat phtest, detail


stcox ib2.art_status2##ct, strata(sex residence2) //best5##
estimates store B
lrtest A B 
estat phtest, detail

log close


//tvc without hiv positive
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_tvc.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)

ta allFixed, gen(allFixed)
des allFixed*

stcox i.allFixed i.sex ct i.residence2
estat phtest, detail 
stcox i.allFixed2 i.allFixed3 i.sex ct i.residence2
estat phtest, detail 

stcox i.allFixed2 i.allFixed3 i.sex ct i.residence2, tvc(allFixed3) texp(_t) nolog

stcox i.allFixed2 i.allFixed3 i.sex ct i.residence2, tvc(allFixed3 sex) texp(_t) nolog

log close

//final best model
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "sp2_log_finalbest.log", replace
use "processed2.dta"
set more off 

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)

ta allFixed, gen(allFixed)

stcox i.allFixed ct, strata(sex) //1.46, p=0.002 #best3#
estimates store A
estat phtest, detail //global p=0.0048. HIV unknown, p=0.0034, is bringing down. it has time dependency

stcox i.allFixed##ct, strata(sex) //best5##
estimates store B
lrtest A B 
estat phtest, detail

stcox i.allFixed2 i.allFixed3 i.sex ct, tvc(allFixed3 sex) texp(_t) nolog

log close
