//alternative approach
//using allFixed

//third version
clear
cd "/Users/sai/OneDrive/Summer Project/stata/"
capture log close
//log using "_allFixed.log", replace
log using "_.log", replace
use "processed2.dta"
set more off 

drop clinic_last_visit_date last_seen_date last_seen_where first_seen_date 
drop years_one fiveyear preg_at_art_start agelastseen age_at_exit2 age_at_exit3
drop hivstale5y neg_yr negAbove5 seroDur missingFixed seroconFixed

ta allFixed, gen(allFixed)

stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)

//ssc install stpm2
//stpm2 fail2, scale(hazard) tvc(allFixed) df(3) dftvc(2) // not working
/*
stpm2 allFixed, df(4) scale(hazard)
predict h, hazard ci
predict s, survival ci
gen h10000=h*10000
twoway (line h10000 _t if allFixed==2 & sex==1, sort clpattern(solid)) ///
		(line h10000 _t if allFixed==1 & sex==1, sort clpattern(longdash)) ///
		(line h10000 _t if allFixed==3 & sex==1, sort clpattern(thin))

twoway (line h10000 _t if allFixed==2 & sex==2, sort clpattern(solid)) ///
		(line h10000 _t if allFixed==1 & sex==2, sort clpattern(longdash)) */
		
stpm2 allFixed2-allFixed3 if sex==1, df(3) scale(hazard)
predict h, hazard ci
//predict s, survival ci
gen h10000=h*10000
twoway (line h10000 _t if allFixed==1, sort clpattern(solid)) ///
		(line h10000 _t if allFixed==2, sort clpattern(longdash)) ///
		(line h10000 _t if allFixed==3, sort clpattern(shortdash)) ///
		, scheme(sj) ytitle("Injury Mortality Rate per 10,000 PY") ///
		xtitle("Age (years)") ///
		title("Injury Mortality Rate in Men") ///
		legend(order(1 "Negative" 2 "Positive" 3 "Unknown") ///
				ring(0) pos(4) cols(1) subtitle("HIV status")) name(h,replace)

stpm2 allFixed2-allFixed3 if sex==2, df(3) scale(hazard)
predict h2, hazard ci
//predict s2, survival ci
gen h100002=h2*10000
twoway (line h100002 _t if allFixed==1, sort clpattern(solid)) ///
		(line h100002 _t if allFixed==2, sort clpattern(longdash)) ///
		(line h100002 _t if allFixed==3, sort clpattern(shortdash)) ///
		, scheme(sj) ylabel(0(10)50) ytitle("Injury Mortality Rate per 10,000 PY") ///
		yscale(range(0(10)50)) ///
		xtitle("Age (years)") ///
		title("Injury Mortality Rate in Women") ///
		legend(order(1 "Negative" 2 "Positive" 3 "Unknown") ///
				ring(0) pos(4) cols(1) subtitle("HIV status")) name(h2,replace)

		/*
stpm2 allFixed2 if sex==2, df(3) scale(hazard)
predict h3, hazard ci
gen h100003=h3*10000
//twoway (line h100003 _t, sort clpattern(solid)) 
twoway (line h100003 _t if allFixed==1, sort clpattern(solid)) ///
		(line h100003 _t if allFixed==2, sort clpattern(longdash)) ///
		(line h100003 _t if allFixed==3, sort clpattern(thin))

stpm2 allFixed2 if sex==1, df(3) scale(hazard)
predict h4, hazard ci
gen h100004=h4*10000
twoway (line h100003 _t, sort clpattern(solid)) 

*/
		
//smooth hazards
gen HIVstatus = allFixed
label define hivstat 1 "Negative" 2 "Positive" 3 "Unknown"
label value HIVstatus hivstat
sts graph,hazard by(HIVstatus)
sts graph if sex==1,hazard by(HIVstatus) yscale(range(0 .008))
sts graph if sex==2,hazard by(HIVstatus) yscale(range(0 .008))

//new ART status variable
stcox i.sex

stcox ct

stcox i.residence2

stcox i.allFixed 

stcox i.allFixed i.sex 

stcox i.allFixed i.sex ct 
estimates store A
estat phtest, detail

stcox i.allFixed i.sex ct i.residence2 // not significant

//separate analysis by sex
stcox i.allFixed if sex==1 //crude male
estimates store M_cr
stcox ct if sex==1 //crude male
stcox i.residence2 if sex==1 //crude male
stcox i.allFixed ct i.residence2 if sex==1 //male, test residence2
estimates store MAR
stcox i.allFixed ct if sex==1 //male
estimates store MA
estat phtest, detail //HIV Unknown violates PH assumption

lrtest MA M_cr //test improvement of ct
lrtest MA MAR //test improvement by residence

stcox i.allFixed##ct if sex==1 //male
estimates store MB
lrtest MA MB //interaction between HIV and calendar time, in male

//tvc HIV Unknown for male
stcox i.allFixed2 i.allFixed3 ct if sex==1, tvc(allFixed3) texp(_t) nolog

stcox i.allFixed if sex==2 //crude female
estimates store F_cr
stcox ct if sex==2 //crude female
stcox i.residence2 if sex==2 //crude female
stcox allFixed ct i.residence2 if sex==2 //female, test residence2
estimates store FAR
stcox i.allFixed ct if sex==2 //female
estimates store FA
estat phtest, detail //HIV Unknown doesn't violate PH assumption

lrtest FA F_cr //test improvement of ct
lrtest FA FAR //test improvement by residence

stcox i.allFixed##ct if sex==2 //female
estimates store FB
lrtest FA FB //interaction between HIV and calendar time, in female


//sex as an interaction covariate

stcox i.allFixed##i.sex ct, nofvlabel
lincom 2.allFixed + 2.allFixed#1.sex, hr //malePosHIV, 1.05 (.86-1.27), p=0.657
lincom 2.allFixed + 2.allFixed#2.sex, hr //femalePosHIV, .71 (.49-1.02), p=.065
lincom 3.allFixed + 3.allFixed#1.sex, hr //maleUnknownHIv, 1.65 (1.21-2.24), p=0.001
lincom 3.allFixed + 3.allFixed#2.sex, hr //femaleUnknownHIv, .94 (.58-1.51), p=.785
//lincom 4.allFixed + 4.allFixed#1.sex, hr //male+ART, 1.4 (.88-2.25), p=.154
//lincom 4.allFixed + 4.allFixed#2.sex, hr //female+ART, 1.1 (.59-2.07), p=.768
estimates store B


lrtest A B //p=.1623, not a strong interaction

//stratification, no longer reported in version 2
//stcox allFixed ct, strata(sex)
//estat phtest, detail

//tvc as a whole is obsolete

//stcox i.art_status2 i.art_status3 i.sex ct, tvc(art_status3 sex) texp(_t) nolog

stcox i.allFixed
estimates store hivbase

stcox i.allFixed i.sex
estimates store hivsex
lrtest hivbase hivsex

log close

