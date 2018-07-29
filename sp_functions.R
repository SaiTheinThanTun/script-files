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
lt <- function(x, per=1, ci95 = FALSE){
  #x is a pyears object
  tableNames <- c('Person-years','Event', 'Rate', 'Low 95%CI', 'High 95%CI')
  py <- x$pyears
  event <- x$event
  rate <- event/py
  lo <- exp(log(rate)-(1.96/sqrt(event)))
  hi <- exp(log(rate)+(1.96/sqrt(event)))
  
  if(ci95==TRUE){
    y <- cbind(py, event, rate*per, lo*per, hi*per)
    colnames(y) <- tableNames
    y
  }
  else{
    y <- cbind(py, event, rate*per)
    colnames(y) <- tableNames[-c(4,5)]
    y
  }
}