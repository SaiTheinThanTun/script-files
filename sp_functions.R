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

#function for associated single decrement life table (ASDT) and cause deleted life table#### 
decom <- function(allcause.A, allcause.B, i_cause.A, i_cause.B, ageint){
  #B-A is done!
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
  
  lx.A <- allcause.A$lx
  nLx.A <- allcause.A$nLx
  Tx.A <- allcause.A$Tx
  
  lx.B <- allcause.B$lx
  nLx.B <- allcause.B$nLx
  Tx.B <- allcause.B$Tx
  
  nmx.A <- allcause.A$Rate
  nRxi.A <- i_cause.A$Event/allcause.A$Event
  nmx.B <- allcause.B$Rate
  nRxi.B <- i_cause.B$Event/allcause.B$Event
  
  ndeltax <- NA
  for(i in 1:length(allcause.A$Event)){
    if(i<length(allcause.A$Event)){
      ndeltax[i] <- ((lx.A[i]/lx.A[1])*((nLx.B[i]/lx.B[i])-(nLx.A[i]/lx.A[i]))) + ((Tx.B[i+1]/lx.A[1])*((lx.A[i]/lx.B[i])-(lx.A[i+1]/lx.B[i+1])))
    }
    else ndeltax[i] <- (lx.A[i]/lx.A[1])*((nLx.B[i]/lx.B[i])-(nLx.A[i]/lx.A[i]))
  }
  ndeltax.i <-ndeltax*(nRxi.B*nmx.B-nRxi.A*nmx.A)/(nmx.B-nmx.A)
  y <- as.data.frame(cbind(lx.A, nLx.A, Tx.A, lx.B, nLx.B, Tx.B, ndeltax,ndeltax.i))
  row.names(y) <- row.names(allcause.A)
  y
}

#function for associated single decrement life table (ASDT) and cause deleted life table, for a list of death events#### 
decomList <- function(allcause.A, allcause.B, i_cause.A, i_cause.B, ageint){
  #B-A is done!
  #decompose by age and cause
  #data needed: LT of 2 groups eg. HIV+ & -, their all cause mortality and # of i cause deaths
  #allcause are dataframe resulted from `lt` function
  #i_cause are 'list' of no. of death events in a population (there will be multiple dataset from each injury cause)
  #ageint defines age interval
  tableNames <- c('Person-years','Event', 'Rate', 'nqx', 'npx', 'lx', 'ndx', 'nLx', 'Tx', 'ex')
  
  #check if the dateset input are correct
  if(sum(colnames(allcause.A) %in% tableNames)!=length(tableNames)) stop("allcause.A dataset incorrect")
  if(sum(colnames(allcause.B) %in% tableNames)!=length(tableNames)) stop("allcause.B dataset incorrect")
  if(length(i_cause.A)!=length(i_cause.B)) stop("i_cause.A and i_cause.B have different length")
  
  #subsetting only necessary sections
  allcause.A <- allcause.A[,which(colnames(allcause.A) %in% tableNames)]
  allcause.B <- allcause.B[,which(colnames(allcause.B) %in% tableNames)]
  
  lx.A <- allcause.A$lx
  nLx.A <- allcause.A$nLx
  Tx.A <- allcause.A$Tx
  
  lx.B <- allcause.B$lx
  nLx.B <- allcause.B$nLx
  Tx.B <- allcause.B$Tx
  
  nmx.A <- allcause.A$Rate
  nmx.B <- allcause.B$Rate
  
  ndeltax <- NA
  for(i in 1:length(allcause.A$Event)){
    if(i<length(allcause.A$Event)){
      ndeltax[i] <- ((lx.A[i]/lx.A[1])*((nLx.B[i]/lx.B[i])-(nLx.A[i]/lx.A[i]))) + ((Tx.B[i+1]/lx.A[1])*((lx.A[i]/lx.B[i])-(lx.A[i+1]/lx.B[i+1])))
    }
    else ndeltax[i] <- (lx.A[i]/lx.A[1])*((nLx.B[i]/lx.B[i])-(nLx.A[i]/lx.A[i]))
  }
  y <- as.data.frame(cbind(lx.A, nLx.A, Tx.A, lx.B, nLx.B, Tx.B, ndeltax))
  row.names(y) <- row.names(allcause.A)
  
  for(j in 1:length(i_cause.A)){
    i_cause.A.event <- i_cause.A[[j]]
    i_cause.B.event <- i_cause.B[[j]]
    nRxi.A <- i_cause.A.event/allcause.A$Event
    nRxi.B <- i_cause.B.event/allcause.B$Event
    ndeltax.i <-ndeltax*(nRxi.B*nmx.B-nRxi.A*nmx.A)/(nmx.B-nmx.A)
    y <- cbind(y,ndeltax.i)
  }
  y
}


#pyears, event rates and CI function####
pyears2 <- function(x, per= 1){
  #x is a pyears object
  x.py <- x$pyears
  x.event <- x$event
  x.rate <- x.event/x.py
  x.lo <- exp(log(x.rate)-(1.96/sqrt(x.event)))
  x.hi <- exp(log(x.rate)+(1.96/sqrt(x.event)))
  y <- as.data.frame(cbind(x.event, x.py, x.rate*per, x.lo*per, x.hi*per))
  colnames(y) <- c("Event", "Person-years", "Rate", "Low 95%CI", "High 95%CI")
  y
}

# allFixed.rate <- pyears(Surv(time=time0, time2 = timex, event = fail2) ~ allFixed, data=dat, scale = 1)
# pyears2(allFixed.rate, per = 10000)
