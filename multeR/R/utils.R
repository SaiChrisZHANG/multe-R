## Auxiliary functions called by multe.R

## Get rid of AsIs class attribute for making sort work
unAsIs <- function(X) {
  if(inherits(X, "AsIs")) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

## Data Constructor
multeData <- function(d,base_val) {
  # Null the global variables, for r CMD check notes
  control <- TbyC_ct <- NULL

  # replace all empty strings to NA
  d[d==""] <- NA

  # no NAs in treatment variable
  if ( FALSE %in% complete.cases(d[[2]]) )
      #message(paste0("Drop ", round(length(d[!complete.cases(d[[2]]),
      #                                       ][[1]])),
      #               " observations with treatment variable missing"))
      d <- d[complete.cases(d[[2]]),]

  # no NAs in control strata indicators
  if ( FALSE %in% complete.cases(d[,-1:-2]) )
      #message(paste0("Drop ", round(length(d[!complete.cases(d[,-1:-2]),
      #                                       ][[1]])),
      #               " observations with control variables missing"))
      d <- d[complete.cases(d[,-1:-2]),]

  # generate a strata indicator variable: index by control variables
  d <- data.table::setDT(d)[,control:=.GRP,by=names(d[,-1:-2])]
  ## generate control by treatment group indicator
  d <- data.table::setDT(d)[,TbyC_ct:=data.table::uniqueN(get(names(d[,2]))
  ), by = control]

  # generate an index column
  d$df_index <- as.numeric(rownames(d))

  # make sure there are no errors in the data
  ## within a control strata, treatment levels smaller than overall
  stopifnot("Found more levels within a control strata than levels overall" =
              max(d$TbyC_ct) <= length(unique(d[[2]])))

  # drop strata with weak overlapping (extreme propensity score)
  if (min(d$TbyC_ct) < length(unique(d[[2]])))
    message(paste0("Drop ",
                   round(length(unique(d[d$TbyC_ct<length(unique(d[[2]])),
                   ]$control)),0),
                   " control strata without sufficient overlapping (",
                   round(length(d[d$TbyC_ct<length(unique(d[[2]])),
                   ]$TbyC_ct), 0), " observations dropped)"))

  d <- d[d$TbyC_ct == length(unique(d[[2]])),]

  # dependent variable
  Y <- d[[1]]
  # treatment
  Tr <- d[[2]]
  # control
  Ctrl <- d[['control']]
  # index
  index <- d$df_index

  ## no observations
  stopifnot("Insufficient number of observations for estimation" = length(Y)>0)

  ## dummify treatment and controls
  if (missing(base_val)) {
    base_val <- unique(Tr)[1]
  }
  Tm <- model.matrix(~relevel(as.factor(Tr),base_val)+0) # treatment, baseline value as first column
  Tlevels <- as.character(c(base_val,unique(Tr)[!unique(Tr)==base_val]))
  Cm <- model.matrix(~as.factor(Ctrl)+0) # control groups

  df <- list(Y=Y, Tr=Tr, Ctrl=Ctrl, var.names=c(names(d)[1:2],"control"),
             Tm = Tm, Tlevels = Tlevels, Cm = Cm, index=index)
  message(paste0(round(length(df$Y)), " observations without missing values ",
                 "in the treatment variable or control variables kept"))
  df
}
