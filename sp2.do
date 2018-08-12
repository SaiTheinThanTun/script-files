clear
cd "/Users/sai/Documents/stata/"
use "processed2.dta"


stset exit, id(idno_original) failure(fail2) time0(entry) origin(time dob) scale(365.25)
//stset exit, id(idno_original) failure(failure) time0(entry) origin(time dob) scale(365.25)

sts graph, by(allFixed) ylabel(minmax)

stcox i.allFixed

stcox i.allFixed i.sex //1.312

stcox i.allFixed i.sex i.agegrp15 //age is already there!

stcox i.residence //overall p=0.0566

stcox i.allFixed i.sex i.period_2011_15 //1.34, p=0.018 #best#

stcox i.allFixed i.sex i.period_2011_15 i.residence //1.34, p=0.038


//testing hazards proportionality
stphplot, by(allFixed)
