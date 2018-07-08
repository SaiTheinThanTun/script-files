//Summer project: MSc Demography and Health
clear
cd "/Users/sai/OneDrive/Summer Project/data"
//log using "~\stata\proj_log.log", replace
//use "~\data\alpha_uMkhanyakude-170601.dta"

//log using "proj_log.log", replace

use "alpha_uMkhanyakude-170601.dta"
format last_neg_date frst_pos_date %td
//drop va_* sr_* study_* clinic_* retro_* c_*
//Exploratory

//codebook //separately saved as codebook.log
//describe


drop age
stsplit age, at(1 2 (1) 122)
browse if idno_original==12 //age check
browse if idno_original==17 //date format checks 
//Data cleaning

//Analysis

//Graphing & Tabulating


//save "alpha_uMkhanyakude-170601.dta", replace
