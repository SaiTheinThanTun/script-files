//alternative approach
//using allFixed

//third version
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
log using "_allFixed.log", replace
use "processed2.dta"
set more off 

drop clinic_last_visit_date last_seen_date last_seen_where first_seen_date 
drop years_one fiveyear preg_at_art_start agelastseen age_at_exit2 age_at_exit3
drop hivstale5y neg_yr negAbove5 seroDur missingFixed seroconFixed

ta art_status2, gen(art_status2)

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)


//new ART status variable
stcox i.sex

stcox ct

stcox residence2

stcox ib2.art_status2 

stcox ib2.art_status2 i.sex 

stcox ib2.art_status2 i.sex ct 
estimates store A
estat phtest, detail

//stcox ib2.art_status2 i.sex ct i.residence2// not significant

//separate analysis by sex
stcox ib2.art_status2 if sex==1 //crude male
estimates store M_cr
stcox ct if sex==1 //crude male
stcox i.residence2 if sex==1 //crude male
//stcox ib2.art_status2 ct i.residence2 if sex==1 //male, test residence2
stcox ib2.art_status2 ct if sex==1 //male
estimates store MA
estat phtest, detail //HIV Unknown violates PH assumption

lrtest MA M_cr //test improvement of ct

stcox ib2.art_status2##ct if sex==1 //male
estimates store MB
lrtest MA MB //interaction between HIV and calendar time, in male

//tvc HIV Unknown for male
stcox i.art_status21 i.art_status23 i.art_status24 ct if sex==1, tvc(art_status21) texp(_t) nolog

stcox ib2.art_status2 if sex==2 //crude female
estimates store F_cr
stcox ct if sex==2 //crude female
stcox i.residence2 if sex==2 //crude female
//stcox ib2.art_status2 ct i.residence2 if sex==2 //female, test residence2
stcox ib2.art_status2 ct if sex==2 //female
estimates store FA
estat phtest, detail //HIV Unknown doesn't violate PH assumption

lrtest FA F_cr //test improvement of ct

stcox ib2.art_status2##ct if sex==2 //female
estimates store FB
lrtest FA FB //interaction between HIV and calendar time, in female


//sex as an interaction covariate

stcox ib2.art_status2##i.sex ct, nofvlabel
lincom 1.art_status2 + 1.art_status2#1.sex, hr //maleUnknownHIV, 1.05 (.86-1.27), p=0.657
lincom 1.art_status2 + 1.art_status2#2.sex, hr //femaleUnknownHIV, .71 (.49-1.02), p=.065
lincom 3.art_status2 + 3.art_status2#1.sex, hr //male+unknownART, 1.65 (1.21-2.24), p=0.001
lincom 3.art_status2 + 3.art_status2#2.sex, hr //female+unknownART, .94 (.58-1.51), p=.785
lincom 4.art_status2 + 4.art_status2#1.sex, hr //male+ART, 1.4 (.88-2.25), p=.154
lincom 4.art_status2 + 4.art_status2#2.sex, hr //female+ART, 1.1 (.59-2.07), p=.768
estimates store B


lrtest A B //p=.1623, not a strong interaction

//stratification, no longer reported in version 2
//stcox ib2.art_status2 ct, strata(sex)
//estat phtest, detail

//tvc as a whole is obsolete

//stcox i.art_status21 i.art_status23 i.art_status24 i.sex ct, tvc(art_status21 sex) texp(_t) nolog

stcox ib2.art_status2
estimates store artbase

stcox ib2.art_status2 i.sex
estimates store artsex
lrtest artbase artsex

log close

