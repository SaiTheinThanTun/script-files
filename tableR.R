tableR <- function(x){
  #function to identify the most frequent 'residence' value
  #taking the earliest one if there are 2 values with equal frequency
  unique_r <- unique(x)
  tr <- rbind(label=unique_r, count=sapply(unique_r,function(y)sum(x==y, na.rm = T)))
  tr[1,which.max(tr[2,])]
}


numbers <- c(4,23,4,23,5,43,54,56,657,67,67,435,
             453,435,324,34,456,56,567,65,34,435, 67, NA)
tableR(numbers)

numbers <- c(4,23,4,23,5,43,54,56,657,435,67,67,
             453,435,324,34,456,56,567,65,34,435, 67)
tableR(numbers)