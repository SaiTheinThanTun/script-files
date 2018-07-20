#survival exercise
#https://bioconnector.org/r-survival.html
library(dplyr)
library(survival)
library(survminer)

lung <- as_tibble(lung)

s <- Surv(lung$time, lung$status)
class(s)
s
head(lung)

sfit <- survfit(Surv(time, status)~1, data=lung)
sfit
summary(sfit)

sfit <- survfit(Surv(time, status)~sex, data=lung)
sfit
summary(sfit)

summary(sfit, times=seq(0, 1000, 100))

ggsurvplot(sfit)


#another package called Epi and popEpi####
xcoh <- structure( list( id = c("A", "B", "C"),
                         birth = c("14/07/1952", "01/04/1954", "10/06/1987"),
                         entry = c("04/08/1965", "08/09/1972", "23/12/1991"),
                         exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                         fail = c(1, 0, 1) ),
                   .Names = c("id", "birth", "entry", "exit", "fail"),
                   row.names = c("1", "2", "3"),
                   class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh <- cal.yr( xcoh, format="%d/%m/%Y", wh=2:4 )
# See how it looks
xcoh
str( xcoh )
# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=entry ),
               exit = list( per=exit, age=exit-birth ),
               exit.status = fail,
               data = xcoh )
Lcoh
# Using character states may have undesired effects:

#tinkering
#another package called Epi and popEpi####
xcoh <- structure( list( id = c("A", "B", "C", "B"),
                         birth = c("14/07/1952", "01/04/1954", "10/06/1987", "01/04/1954"),
                         entry = c("04/08/1965", "08/09/1972", "23/12/1991", "23/05/1995"),
                         exit = c("27/06/1997", "23/05/1995", "24/07/1998", "23/10/1995"),
                         fail = c(1, 0, 1, 1) ),
                   .Names = c("id", "birth", "entry", "exit", "fail"),
                   row.names = c("1", "2", "3", "4"),
                   class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh <- cal.yr( xcoh, format="%d/%m/%Y", wh=2:4 )
# See how it looks
xcoh
str( xcoh )
# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( per=entry ),
               exit = list( per=exit, age=exit-birth ),
               exit.status = fail,
               data = xcoh )
Lcoh

if(FALSE){
  #viewing the code for Lexis function
function (entry, exit, duration, entry.status = 0, exit.status = 0, 
          id, data, merge = TRUE, states, tol = .Machine$double.eps^0.5, 
          keep.dropped = FALSE) 
{
  nmissing <- missing(entry) + missing(exit) + missing(duration)
  if (nmissing > 2) 
    stop("At least one of the arguments exit and duration must be supplied")
  only.exit <- missing(entry.status) && !missing(exit.status)
  if (!missing(data)) {
    if (!missing(entry)) {
      entry <- eval(substitute(entry), data, parent.frame())
    }
    if (!missing(exit)) {
      exit <- eval(substitute(exit), data, parent.frame())
    }
    if (!missing(duration)) {
      duration <- eval(substitute(duration), data, parent.frame())
    }
    entry.status <- eval(substitute(entry.status), data, 
                         parent.frame())
    exit.status <- eval(substitute(exit.status), data, parent.frame())
    if (!missing(id)) {
      id <- eval(substitute(id), data, parent.frame())
    }
    if (merge) {
      data <- as.data.frame(data)
    }
  }
  wh.miss <- any(is.na(entry.status)) + 2 * any(is.na(exit.status))
  if (wh.miss > 0) 
    stop("Missing values in ", switch(wh.miss, "entry status", 
                                      "exit status", "entry AND exit status"))
  if (only.exit) {
    if (is.logical(exit.status)) {
      entry.status <- FALSE
      cat("NOTE: entry.status has been set to FALSE for all.\n")
    }
    if (is.character(exit.status)) {
      exit.status <- factor(exit.status)
    }
    if (is.factor(exit.status)) {
      entry.status <- factor(rep(levels(exit.status)[1], 
                                 length(exit.status)), levels = levels(exit.status), 
                             labels = levels(exit.status))
      cat("NOTE: entry.status has been set to", paste("\"", 
                                                      levels(exit.status)[1], "\"", sep = ""), "for all.\n")
    }
    if (is.numeric(exit.status)) {
      entry.status <- rep(0, length(exit.status))
      cat("NOTE: entry.status has been set to 0 for all.\n")
    }
  }
  if (is.character(entry.status)) 
    entry.status <- factor(entry.status)
  if (is.character(exit.status)) 
    exit.status <- factor(exit.status)
  if (is.factor(entry.status) || is.factor(exit.status)) {
    if (is.factor(entry.status) && is.factor(exit.status)) {
      if (!identical(levels(entry.status), levels(exit.status))) {
        all.levels = union(levels(entry.status), levels(exit.status))
        entry.status <- factor(entry.status, levels = all.levels)
        exit.status <- factor(exit.status, levels = all.levels)
        cat("Incompatible factor levels in entry.status and exit.status:\n", 
            "both lex.Cst and lex.Xst now have levels:\n", 
            all.levels, "\n")
      }
    }
    else {
      stop("Incompatible classes for entry and exit status")
    }
  }
  else {
    if (mode(entry.status) != mode(exit.status)) {
      stop("Incompatible mode for entry and exit status")
    }
  }
  if (nmissing == 2) {
    if (!missing(exit)) {
      if (length(exit) > 1) 
        stop("If 'entry' is omitted, only one timescale can be specified.")
      else {
        entry <- exit
        entry[[1]] <- 0 * entry[[1]]
        cat("NOTE: entry is assumed to be 0 on the", 
            names(exit), "timescale.\n")
      }
    }
    else if (!missing(duration)) {
      if (length(duration) > 1) 
        stop("If 'entry' is omitted, only one timescale can be specified")
      else {
        entry <- duration
        entry[[1]] <- 0 * entry[[1]]
        cat("NOTE: entry is assumed to be 0 on the", 
            names(duration), "timescale.\n")
      }
    }
    else stop("Either exit or duration must be supplied.")
  }
  if (!missing(entry)) {
    entry <- as.data.frame(entry)
    if (is.null(names(entry))) 
      stop("entry times have no names")
    if (any(substr(names(entry), 1, 4) == "lex.")) 
      stop("names starting with \"lex.\" cannot be used for time scales")
  }
  if (!missing(exit)) {
    exit <- as.data.frame(exit)
    if (is.null(names(exit))) 
      stop("exit times have no names")
    if (any(substr(names(exit), 1, 4) == "lex.")) 
      stop("names starting with \"lex.\" cannot be used for time scales")
  }
  if (!missing(duration)) {
    duration <- as.data.frame(duration)
    if (is.null(names(duration))) 
      stop("duration have no names")
    if (any(substr(names(duration), 1, 4) == "lex.")) 
      stop("names starting with \"lex.\" cannot be used for time scales")
  }
  if (missing(entry)) {
    entry <- exit - duration
  }
  if (missing(duration)) {
    full.time.scales <- intersect(names(entry), names(exit))
    if (length(full.time.scales) == 0) {
      stop("Cannot calculate duration from entry and exit times")
    }
    duration <- exit[, full.time.scales[1]] - entry[, full.time.scales[1]]
  }
  if (missing(exit)) {
    all.time.scales <- names(entry)
  }
  else {
    all.time.scales <- unique(c(names(entry), names(exit)))
    entry.missing <- setdiff(all.time.scales, names(entry))
    if (length(entry.missing) > 0) {
      entry <- cbind(entry, exit[, entry.missing, drop = FALSE] - 
                       duration)
    }
    dura <- exit - entry[, names(exit), drop = FALSE]
    if (missing(duration)) {
      duration <- dura[, 1]
    }
    ok <- sapply(lapply(dura, all.equal, duration), isTRUE)
    if (!all(ok)) {
      stop("Duration is not the same on all time scales")
    }
  }
  if (missing(id)) {
    id <- 1:nrow(entry)
  }
  else if (any(duplicated(id))) {
  }
  if (is.data.frame(duration)) 
    duration <- duration[, 1]
  lex <- data.frame(entry, lex.dur = duration, lex.Cst = entry.status, 
                    lex.Xst = exit.status, lex.id = id)
  if (!missing(states)) {
    st.lev <- sort(unique(as.character(c(lex$lex.Cst, lex$lex.Xst))))
    lex$lex.Cst <- factor(as.character(lex$lex.Cst), levels = st.lev, 
                          labels = states)
    lex$lex.Xst <- factor(as.character(lex$lex.Xst), levels = st.lev, 
                          labels = states)
  }
  if (!missing(data) && merge) {
    duplicate.names <- intersect(names(lex), names(data))
    if (length(duplicate.names) > 0) {
      stop("Cannot merge data with duplicate names:", paste(duplicate.names, 
                                                            collapse = " "))
    }
    lex <- cbind(lex, data)
  }
  short.dur <- lex$lex.dur <= tol
  if (any(short.dur)) {
    warning("Dropping ", sum(short.dur), " rows with duration of follow up < tol\n", 
            if (keep.dropped) 
              "  The dropped rows are in the attribute 'dropped'\n", 
            if (keep.dropped) 
              "  To see them type attr(Obj,'dropped'),\n", 
            if (keep.dropped) 
              "  to get rid of them type: attr(Obj,'dropped') <- NULL\n", 
            if (keep.dropped) 
              "  - where 'Obj' is the name of your Lexis object")
    lex <- subset(lex, !short.dur)
    if (keep.dropped) 
      attr(lex, "dropped") <- subset(data, short.dur)
  }
  attr(lex, "time.scales") <- all.time.scales
  attr(lex, "time.since") <- rep("", length(all.time.scales))
  breaks <- vector("list", length(all.time.scales))
  names(breaks) <- all.time.scales
  attr(lex, "breaks") <- breaks
  class(lex) <- c("Lexis", class(lex))
  return(lex)
}
<bytecode: 0x1139925f0>
  <environment: namespace:Epi>
  }