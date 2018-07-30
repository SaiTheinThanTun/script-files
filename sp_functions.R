#functions for summer project
####function to extract a record####
reveal <- function(x){
  # View(liteumk[liteumk$idno_original==x,])
  View(dat[dat$idno_original==x,])
}
reveal2 <- function(x,y=NA){
  # View(liteumk[liteumk$idno_original==x,])
  if(is.na(y)){
    View(dat[dat$idno_original==x,])
  }
  else
    View(dat[dat$idno_original==x,which(names(dat) %in% y)])
}
####function to identify frequentest value####
tableR <- function(x){
  #function to identify the most frequent 'residence' value
  #taking the earliest one if there are 2 values with equal frequency
  unique_r <- unique(x)
  tr <- rbind(label=unique_r, count=sapply(unique_r,function(y)sum(x==y, na.rm = T)))
  tr[1,which.max(tr[2,])]
}

#mortality rate (lifetable) extraction####
lt <- function(x, per = 1, ci95 = FALSE, nax = .5, ageint){
  #x is a pyears object with unabridged (single age) ASMR
  #ageint=1 default has been removed
  tableNames <- c('Person-years','Event', 'Rate', 'Low 95%CI', 'High 95%CI', 'nqx', 'npx', 'lx', 'ndx', 'nLx', 'Tx', 'ex')
  py <- x$pyears
  event <- x$event
  rate <- event/py
  lo <- exp(log(rate)-(1.96/sqrt(event)))
  hi <- exp(log(rate)+(1.96/sqrt(event)))
  nqx <- (rate*ageint)/(1+(1-nax)*rate*ageint)
  nqx[length(event)] <- 1 #last nqx for open ended interval must be 1!
  npx <- 1-nqx
  
  lx <- NA
  for(i in 1:length(event)){
    if(i==1){lx[i] <- 1} #radix of 1 so that ex can be directly got from Tx
    else {lx[i] <- lx[i-1]*npx[i-1]}
  }
  
  ndx <- lx*nqx
  
  nLx <- NA
  for(i in 1:length(event)){
    if(i<length(event)){nLx[i] <- (lx[i+1]*ageint)+(ageint*nax*ndx[i])}
    else {
      if(rate[i]!=0) nLx[i] <- lx[i]/rate[i]
      else nLx[i] <- 0
      }
  }
  
  Tx <- NA
  for(i in 1:length(event)){
    Tx[i] <- sum(nLx[i:length(event)])
  }
  
  ex <- Tx/lx
  
  if(ci95==TRUE){
    y <- as.data.frame(cbind(py, event, rate*per, lo*per, hi*per, nqx, npx, lx, ndx, nLx, Tx, ex))
    colnames(y) <- tableNames
    y
  }
  else{
    y <- as.data.frame(cbind(py, event, rate*per, nqx, npx, lx, ndx, nLx, Tx, ex))
    colnames(y) <- tableNames[-c(4,5)]
    y
  }
}

#lt(asmr.all, ageint = 1)

#function for associated single decrement life table (ASDT) and cause deleted life table 
asdt <- function(allcause, i_cause, ageint, deletion=TRUE){
  #Associated single decrement life table
  #allcause and i_cause are dataframe resulted from `lt` function
  #ageint defines age interval
  #deletion specifies if othercause (allcause `minus` i_cause) is used
  tableNames <- c('Person-years','Event', 'Rate', 'nqx', 'npx', 'lx', 'ndx', 'nLx', 'Tx', 'ex')
  
  #check if the dateset input are correct
  if(sum(colnames(allcause) %in% tableNames)!=length(tableNames)) stop("allcause dataset incorrect")
  if(sum(colnames(i_cause) %in% tableNames)!=length(tableNames)) stop("i_cause dataset incorrect")
  
  #subsetting only necessary sections
  allcause <- allcause[,which(colnames(allcause) %in% tableNames)]
  i_cause <- i_cause[,which(colnames(i_cause) %in% tableNames)]
  
  if(deletion){
    new_nmx <- allcause$Rate - i_cause$Rate
  }
  else {new_nmx <- i_cause$Rate}
  
  s_npx <- exp(-ageint*new_nmx) #s denotes star or asterisk*
  s_npx[length(allcause$Event)] <- 0
  
  s_lx <- NA
  for(i in 1:length(allcause$Event)){
    if(i==1){s_lx[i] <- 1} #radix of 1 so that ex can be directly got from Tx
    else {s_lx[i] <- s_lx[i-1]*s_npx[i-1]}
  }
  
  s_ndx <- s_lx*(1-s_npx)
  
  s_nLx <- s_ndx/new_nmx
  
  s_Tx <- NA
  for(i in 1:length(allcause$Event)){
    s_Tx[i] <- sum(s_nLx[i:length(allcause$Event)])
  }
  
  s_ex <- s_Tx/s_lx
  y <- as.data.frame(cbind(allcause, new_nmx, s_npx, s_lx, s_ndx, s_nLx, s_Tx, s_ex))
  y
}

#asdt(allcause = lt.all, i_cause = lt.ext, ageint = 1)

#function for associated single decrement life table (ASDT) and cause deleted life table 
decom <- function(allcause.A, i_cause.A, allcause.B, i_cause.B, ageint){
  #decompose by age and cause
  #data needed: LT of 2 groups eg. HIV+ & -, their all cause mortality and # of i cause deaths
  #allcause and i_cause are dataframe resulted from `lt` function
  #ageint defines age interval
  tableNames <- c('Person-years','Event', 'Rate', 'nqx', 'npx', 'lx', 'ndx', 'nLx', 'Tx', 'ex')
  
  #check if the dateset input are correct
  if(sum(colnames(allcause.A) %in% tableNames)!=length(tableNames)) stop("allcause.A dataset incorrect")
  if(sum(colnames(i_cause.A) %in% tableNames)!=length(tableNames)) stop("i_cause.A dataset incorrect")
  if(sum(colnames(allcause.B) %in% tableNames)!=length(tableNames)) stop("allcause.B dataset incorrect")
  if(sum(colnames(i_cause.B) %in% tableNames)!=length(tableNames)) stop("i_cause.B dataset incorrect")
  
  #subsetting only necessary sections
  allcause.A <- allcause.A[,which(colnames(allcause.A) %in% tableNames)]
  i_cause.A <- i_cause.A[,which(colnames(i_cause.A) %in% tableNames)]
  allcause.B <- allcause.B[,which(colnames(allcause.B) %in% tableNames)]
  i_cause.B <- i_cause.B[,which(colnames(i_cause.B) %in% tableNames)]
  
  ndeltax <- 
}