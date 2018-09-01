/*************************************************************************
Breast cancer PH analysis for Stata Journal
Paul Lambert (paul.lambert@le.ac.uk)
*************************************************************************/
clear all
use "breast_under50"

stset survtime, failure(dead == 1) exit(time 5) id(ident)

/* Cox Model */
stcox dep2-dep5, noshow nolog 

/* stpm2 PH model */
stpm2 dep2-dep5, df(5) scale(hazard) eform nolog

/* predictions */
predict xb, xb 
predict s, survival
predict h, hazard 

/* survival */
twoway 	(line s _t  if dep1==1, sort clpattern(solid) clwidth(thin)) ///
		(line s _t if dep2==1, sort clpattern(longdash) clwidth(thin)) ///
		(line s _t if dep3==1, sort clpattern(dash) clwidth(thin)) ///
		(line s _t if dep4==1, sort clpattern(shortdash) clwidth(thin)) ///
		(line s _t if dep5==1, sort clpattern(dot) clwidth(thin)) ///
		, scheme(sj) ylabel(,angle(h)) ytitle(Survival) ///
		xtitle("Time from Diagnosis (years)") ///
		title("(c)") ///
		legend(order(1 "Least Deprived" 2 "2" 3 "3" 4 "4" 5 "Most Deprived") ///
				ring(0) pos(7) cols(1) subtitle("Deprivation Group") bmargin(large) nobox) ///
		name(s,replace)

/* log cumulative hazard */
twoway (line xb _t if dep1==1, sort clpattern(solid) clwidth(thin)) ///
		(line xb _t if dep2==1, sort clpattern(longdash) clwidth(thin)) ///
		(line xb _t if dep3==1, sort clpattern(dash) clwidth(thin)) ///
		(line xb _t if dep4==1, sort clpattern(shortdash) clwidth(thin)) ///
		(line xb _t if dep5==1, sort clpattern(dot) clwidth(thin)) ///
		, scheme(sj) ylabel(,angle(h)) ytitle(log cumulative hazard) ///
		xtitle("Time from Diagnosis (years)") ///
		title("(a)") ///
		legend(order(1 "Least Deprived" 2 "2" 3 "3" 4 "4" 5 "Most Deprived") ///
				ring(0) pos(5) cols(5) ///
				subtitle("Deprivation Group", size(*0.7)) size(*0.8) nobox ///
				region(lstyle(none))) ///
		name(xb1,replace)

/* log cumulative hazard  vs log(time) */
twoway (line xb _t if dep1==1, sort clwidth(thin) clpattern(solid)) ///
		(line xb _t if dep2==1, sort clwidth(thin) clpattern(longdash)) ///
		(line xb _t if dep3==1, sort clwidth(thin) clpattern(dash)) ///
		(line xb _t if dep4==1, sort clwidth(thin) clpattern(shortdash)) ///
		(line xb _t if dep5==1, sort clwidth(thin) clpattern(dot)) ///
		, scheme(sj) ylabel(,angle(h)) ytitle(log cumulative hazard) ///
		xtitle("Time from Diagnosis (years)") xscale(log) ///
		title("(b)") ///
		legend(order(1 "Least Deprived" 2 "2" 3 "3" 4 "4" 5 "Most Deprived") ///
				ring(0) pos(5) cols(1) subtitle("Deprivation Group")) name(xb2,replace)

/* hazard */
gen h1000 = h * 1000

twoway (line h1000 _t if dep1==1, sort clpattern(solid)) ///
		(line h1000 _t if dep2==1, sort clpattern(longdash)) ///
		(line h1000 _t if dep3==1, sort clpattern(dash)) ///
		(line h1000 _t if dep4==1, sort clpattern(shortdash)) ///
		(line h1000 _t if dep5==1, sort clpattern(dot)) ///
		, scheme(sj) ylabel(0(25)150,angle(h)) ytitle("Mortality Rate (per 1000 py)") ///
		yscale(range(0(25)150)) ///
		xtitle("Time from Diagnosis (years)") ///
		title("(d)") ///
		legend(order(1 "Least Deprived" 2 "2" 3 "3" 4 "4" 5 "Most Deprived") ///
				ring(0) pos(1) cols(1) subtitle("Deprivation Group")) name(h,replace)

/* only keep deprivation groups 1 and 5 */
drop if inrange(dep,2,4)

/* Time-dependent effect (cumulative hazard scale) */
/* define timevar for predictions */
range timevar 0.003 5 500

/* PH model */
stpm2 dep5, scale(hazard) df(5) nolog
predict h_ph, hazard timevar(timevar)

stpm2 dep5, df(5) scale(hazard) tvc(dep5) dftvc(3) nolog
predict h_tvc, hazard timevar(timevar)

twoway	(line h_ph timevar if dep5 == 0 , sort clcolor(gs8) clpattern(solid) clwidth(thin)) ///
		(line h_ph timevar if dep5 == 1 , sort clcolor(gs8) clpattern(dash) clwidth(thin)) ///
		(line h_tvc timevar if dep5 == 0 , sort clcolor(black) clpattern(solid) clwidth(thick)) ///
		(line h_tvc timevar if dep5 == 1 , sort clcolor(black) clpattern(dash) clwidth(thick)) ///
		, scheme(sj) ///
		legend(order(3 "Least deprived" 4 "Most deprived") ring(0) pos(1) cols(1) ///
				nobox region(lstyle(none))) ///
		ylabel(,angle(h)) ytitle(hazard rate) yscale(log) ///
			xtitle("Time from Diagnosis (years)") ///
		note("Thinner lines are predictions from proportional hazards model")

/* predict hazard ratio, hazard difference and survival difference */
predict hr, hrnum(dep5 1) hrdenom(dep5 0) timevar(timevar) ci
predict hdiff, hdiff1(dep5 1) hdiff2(dep5 0) timevar(timevar) ci
predict sdiff, sdiff1(dep5 1) sdiff2(dep5 0) timevar(timevar) ci

/* hazard ratio */
twoway	(rarea hr_lci hr_uci timevar if hr_uci<4, sort pstyle(ci)) ///
		(line hr timevar, sort clpattern(solid) clwidth(thick)) ///
		, scheme(sj) ///
		legend(off) ///
		ylabel(,angle(h)) ytitle(hazard ratio) ///
		xtitle("Time from Diagnosis (years)") ///
		title("(a)") ///
		yline(1, lpattern(dash) lwidth(thin)) name(hr,replace)

/* hazard difference */
replace hdiff = hdiff*1000
replace hdiff_lci = hdiff_lci*1000
replace hdiff_uci = hdiff_uci*1000
twoway	(rarea hdiff_lci hdiff_uci timevar, sort pstyle(ci)) ///
		(line hdiff timevar, sort clpattern(solid) clwidth(thick)) ///
		, scheme(sj) ///
		legend(off) ///
		ylabel(,angle(h)) ytitle("Difference in hazard rate (per 1000 py's)") ///
		xtitle("Time from Diagnosis (years)") ///
		title("(b)") ///
		yline(0, lpattern(dash) lwidth(thin))  name(hd,replace)


/* survival */
predict s0, s zeros timevar(timevar)
predict s1, s at(dep5 1) timevar(timevar)
twoway	(line s0 timevar, sort clpattern(solid) clwidth(thick)) ///
		(line s1 timevar, sort clpattern(dash) clwidth(thick)) ///
		, scheme(sj) ///
		legend(order(1 "Least deprived" 2 "Most deprived") ring(0) pos(7) cols(1) ///
				nobox region(lstyle(none))) ///
		ylabel(0.4(0.2)1,angle(h)) ytitle("Survival") ///
		title("(c)") ///
		xtitle("Time from Diagnosis (years)") name(s,replace)

/* survival difference */
twoway	(rarea sdiff_lci sdiff_uci timevar, sort pstyle(ci)) ///
		(line sdiff timevar, sort clpattern(solid) clwidth(thick)) ///
		, scheme(sj) ///
		legend(off) ///
		ylabel(-0.1(0.02)0.02,angle(h)) ytitle("Difference in Survival Curves") ///
		xtitle("Time from Diagnosis (years)") ///
		title("(d)") ///
		yline(0, lpattern(dash) lwidth(thin)) name(sd,replace)

/* combine graphs */
graph combine hr hd s sd, nocopies scheme(sj)

/* difference df for time-dependent effect */
forvalues i = 1/5 {
	stpm2 dep5, scale(hazard) df(5) tvc(dep5) dftvc(`i')
	estimates store df`i'
	predict hr`i', hrnum(dep5 1) timevar(timevar) ci
	count if _d == 1
}

/* plot time-dependent hazard ratios */
twoway	(line hr1 timevar) ///
		(line hr2 timevar) ///
		(line hr3 timevar) ///
		(line hr4 timevar) ///
		(line hr5 timevar) ///						
		, legend(order(1 "`df1'" 2 "`df2'" 3 "`df3'" 4 "`df4'" 5 "`df5'") ring(0) pos(1) cols(1) ///
			nobox region(lstyle(none))) ///
		scheme(sj) ///
		ylabel(,angle(h)) ytitle(hazard ratio) ///
		xtitle("Time from Diagnosis (years)") 

/* Two time-dependent effects */
stpm2 dep5 yeardiag, scale(hazard) df(5) tvc(dep5 yeardiag) dftvc(3) nolog

predict hr_early, hrnum(dep5 1 yeardiag 1986) hrdenom(dep5 0 yeardiag 1986) ///
		timevar(timevar) ci
predict hr_late, hrnum(dep5 1 yeardiag 1990) hrdenom(dep5 0 yeardiag 1990) ///
		timevar(timevar) ci
								
twoway	(line hr_early hr_early_lci hr_early_uci timevar, clpattern(solid solid solid) clwidth(thick)) ///
		(line hr_late hr_late_lci hr_late_uci timevar, clpattern(dash dash dash) clwidth(thick)) ///
		, legend(order(1 "1986" 4 "1990") ring(0) pos(1) cols(1) ///
			nobox region(lstyle(none))) ///
		scheme(sj) yscale(log) ///
		ylabel(,angle(h)) ytitle(hazard ratio) ///
		xtitle("Time from Diagnosis (years)") ///
		title("Hazard Ratio for Deprivation Group")
