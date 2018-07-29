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
lt <- function(x, per = 1, ci95 = FALSE, nax = .5, ageint = 1){
  #x is a pyears object with unabridged (single age) ASMR
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
    y <- cbind(py, event, rate*per, lo*per, hi*per, nqx, npx, lx, ndx, nLx, Tx, ex)
    colnames(y) <- tableNames
    y
  }
  else{
    y <- cbind(py, event, rate*per, nqx, npx, lx, ndx, nLx, Tx, ex)
    colnames(y) <- tableNames[-c(4,5)]
    y
  }
}
#things to fix
#last nqx should be 1!
lt(asmr.all)
