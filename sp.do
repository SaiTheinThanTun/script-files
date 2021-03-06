//Summer project: MSc Demography and Health
clear
cd "/Users/sai/OneDrive/Summer Project/data"
//log using "~\stata\proj_log.log", replace
//use "~\data\alpha_uMkhanyakude-170601.dta"

//log using "proj_log.log", replace

//use "alpha_uMkhanyakude-170601.dta" //2nd original version with stsplit by age

* new hivstat variable based on 5 years post-negative follow up time.
/*
gen hivstale5y = hivstatus_broad
replace hivstale5y=3 if hivstatus_detail == 4  //Post Negative < year
replace hivstale5y=1 if hivstale5y==3 &(hivstatus_detail==4 | hivstatus_detail ==8) & exit <= last_neg_date + 365.25*5 & last_neg_date !=.

label define hivstale5y 1 "Negative" 2 "Positive" 3 "Unknown"
label value hivstale5y hivstale5y

gen neg_yr= (exit-last_neg_date)/365.25
gen negAbove5=0 if neg_yr<=5 //for checking ! careful b/c . is the largest value in Stata!!!!
replace negAbove5=1 if neg_yr!=. & neg_yr>5

drop if _d==. //this is also done in R
gen seroDur = (frst_pos_date - last_neg_date)/365.25

//positive status
//browse if hivstale5y==. & frst_pos_date!=. & (exit-frst_pos_date>0) //should be zero record

//allFixed

gen allFixed = hivstale5y
ta negAbove5 allFixed, mi
replace allFixed=1 if neg_yr <= 5 & neg_yr >0 & allFixed==. & frst_pos_date==. //23352
replace allFixed=3 if neg_yr!=. & neg_yr > 5 & allFixed==. & frst_pos_date==. //8981 becareful of greater than sign

replace allFixed = 1 if hivstatus_detail==7 & (neg_yr <= seroDur/2)
replace allFixed = 2 if hivstatus_detail==7 & (neg_yr > seroDur/2)

label define allFixed 1 "Negative" 2 "Positive" 3 "Unknown"
label value allFixed allFixed
ta negAbove5 allFixed, mi


//missingFixed
gen missingFixed=hivstale5y
ta negAbove5 missingFixed, mi
replace missingFixed=1 if neg_yr <= 5 & neg_yr >0 & missingFixed==. & frst_pos_date==. //23352
replace missingFixed=3 if neg_yr!=. & neg_yr > 5 & missingFixed==. & frst_pos_date==. //8981 becareful of greater than sign
label value missingFixed hivstatus_broad


label define missingFixed 1 "Negative" 2 "Positive" 3 "Unknown"
label value missingFixed missingFixed
ta negAbove5 missingFixed, mi


//seroconFixed
gen seroconFixed=hivstale5y
replace seroconFixed = 1 if hivstatus_detail==7 & (neg_yr <= seroDur/2)
replace seroconFixed = 2 if hivstatus_detail==7 & (neg_yr > seroDur/2)
label value seroconFixed hivstatus_broad

label define seroconFixed 1 "Negative" 2 "Positive" 3 "Unknown"
label value seroconFixed seroconFixed
ta negAbove5 seroconFixed, mi


save "alpha_uMkhanyakude-v4.dta", replace //4th version 20180725
*/

use "alpha_uMkhanyakude-v4.dta" //4th version


// issue with hivstatus
//keep idno_original last_neg_date frst_pos_date exit entry hivstatus_detail hivstatus_broad agegrp sex //original
keep idno_original last_neg_date frst_pos_date exit entry hivstatus_detail hivstatus_broad agegrp sex hivstale5y neg_yr seroDur negAbove5


//gen neg_yr= (exit-last_neg_date)/365.25
//drop if neg_yr<0
//drop if last_neg_date==.
//gen negAbove5=1 if neg_yr>5
//replace negAbove5=0 if negAbove5==.
//summ negAbove5
//ta negAbove5 hivstatus_broad, mi //original
ta negAbove5 hivstale5y, mi
ta hivstatus_broad hivstale5y

browse if hivstale5y!=hivstatus_broad

//gen negAbove1=1 if neg_yr>1 & neg_yr!=.
//replace negAbove1=0 if neg_yr<1
ta negAbove1 hivstatus_broad, mi

browse if neg_yr<1 & hivstatus_broad==3 //seroconverter

browse if hivstatus_detail==7 & hivstatus_broad==1 //seroconverter negative

browse if neg_yr==1 //19131

browse if neg_yr==2 //1452

browse if neg_yr==3 //790

browse if neg_yr==4 //40

browse if neg_yr==5 // 14

browse if neg_yr==6 //11, doesn't really happen

browse if idno_original==47 //236


//format last_neg_date frst_pos_date %td

//drop va_* sr_* study_* clinic_* retro_* c_*
//Exploratory

//codebook //separately saved as codebook.log
//describe

/*
drop age
stsplit age, at(1 2 (1) 122)
browse if idno_original==12 //age check
browse if idno_original==17 //date format checks 
*/

//Data cleaning

//Analysis

//Graphing & Tabulating


//save "alpha_uMkhanyakude-170601.dta", replace

stset
ta hivstatus_broad, missing

recode hivstatus_broad (.=3)

strate hivstatus_broad

/*
//backup
ta negAbove5 hivstale5y, mi
replace hivstale5y=1 if neg_yr <= 5 & neg_yr >0 & hivstale5y==. & frst_pos_date==. //23352
//replace hivstale5y=3 if neg_yr > 5 & hivstale5y==. & frst_pos_date==. //8981 becareful of greater than sign
replace hivstale5y=3 if neg_yr!=. & neg_yr > 5 & hivstale5y==. & frst_pos_date==. //8981 becareful of greater than sign
ta negAbove5 hivstale5y, mi

//seroconverters
//tab hivstale5y if hivstatus_detail== 7 //seroconverters
gen seroDur = (frst_pos_date - last_neg_date)/365.25
replace hivstale5y = 1 if hivstatus_detail==7 & (neg_yr <= seroDur/2)
replace hivstale5y = 2 if hivstatus_detail==7 & (neg_yr > seroDur/2)
//tab hivstatus_broad hivstale5y
//browse if hivstatus_detail== 7
//ta negAbove5 hivstale5y, mi
*/
