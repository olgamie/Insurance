PVcashflow = function(qxt, ageStart, omegaAge, pensionAge, valyear, ir, type = "annuity"){
  if(!is.matrix(qxt))
    stop("qxt need to be a matrix")
  if(!is.numeric(ageStart) || length(ageStart)!=1 || ageStart<0)
    stop("ageStart has to be nonnagative, numeric of length 1")
  if(!is.numeric(omegaAge) || length(omegaAge)!=1 || omegaAge<0)
    stop("omegaAge has to be nonnagative, numeric of length 1")
  if(!is.numeric(pensionAge) || length(pensionAge)!=1 || pensionAge<0)
    stop("pensionAge has to be nonnegative, numeric of length 1")
  if(!is.numeric(valyear) || length(valyear)!=1 || valyear<0)
    stop("valyear has to be numeric of length 1")
  if(!is.numeric(ir) || length(ir)!=1)
    stop("ir has to be numeric of length 1")
  if((type!="premium" && type!="annuity") || length(type)!=1)
    stop("Type needs be either premium or annuity")
  if(ageStart > omegaAge)
    stop("ageStart cannot be larger than omegaAge")
  if(pensionAge > omegaAge)
    stop("pensionAge cannot be larger than omegaAge")
  
  if (type == "annuity"){
    CF <- c(rep(0,max(0,(pensionAge-ageStart+1))), rep(1,(omegaAge-max(pensionAge, ageStart-1))))
  }
  if(type=="premium"){
    CF <- c(rep(1,max(0,(pensionAge-ageStart+1))), rep(0,(omegaAge-max(pensionAge, ageStart-1))))
  }
  
  if(is.null(colnames(qxt)) || is.null(rownames(qxt)))
    stop("Matrix qxt needs to have colnames and rownames, at least one of the elements is missing")
  if(!toString(ageStart) %in% rownames(qxt))
    stop("ageStart value needs to have corresponding row in the qxt matrix")
  if(!toString(omegaAge) %in% rownames(qxt))
    stop("omegaAge value needs to have corresponding row in the qxt matrix")
  if(!toString(valyear) %in% colnames(qxt))
    stop("valyear value needs to have corresponding column in the qxt matrix")
  if(!toString(valyear+(omegaAge-ageStart)) %in% colnames(qxt))
    stop("Qxt matrix has insufficient number of columns to find probabilities up to omageAge")
  
  t_start <- which(colnames(qxt) == toString(valyear))
  x_start <- which(rownames(qxt) == toString(ageStart))
  
  qx <- diag(qxt[x_start:(x_start+(omegaAge-ageStart)), t_start:(t_start+(omegaAge-ageStart))])
  pxt <- cumprod(1- qx)
  k <- 0:(omegaAge-ageStart)
  interestRates <- rep(ir, (omegaAge-ageStart)+1)
  v <- (1 + interestRates)^-k
  pvCF <- sum(pxt*v*CF)
  pvCF
}

kannistoExtrapolation = function(qxt, ages, years, maxAge=120, nObs=15){
  stopifnot(is.matrix(qxt), is.numeric(ages), is.numeric(years))
  stopifnot(is.numeric(maxAge), length(maxAge)==1, maxAge>0, length(nObs)==1, is.numeric(nObs), nObs>0)
  stopifnot(length(ages)>=nObs)
  stopifnot(colnames(qxt)==as.character(years))
  
  ## calculate ages to be extrapolated
  agesExtrap <- (tail(ages,1)+1):maxAge
  
  ## extract the ages to be used for the regression
  obsAges <- tail(ages, nObs)
  
  #extrapolate for each year
  extrapolateQxt <- sapply(years, function(x){
    ## extract the mortality rates and ages which are going to be used for the regression
    obsQx <- pmin(tail(qxt[,toString(x)], nObs),1)
    obsQxt  <- obsQx
    ## transform the mortality rates
    obsQx <- -log(1 - obsQx)
    obsQx <- log(obsQx/(1 - obsQx))
    
    if (length(na.omit(obsQx))>2){
      ## apply the regression and obtain phi1 and phi2
      model <- lm(formula = obsQx~obsAges)
      phi1 <- exp(model$coefficient[1])
      phi2 <- model$coefficient[2]
      
      #extrapolate mortality rates and add to the overall matrix
      1-exp(-phi1*exp(phi2*agesExtrap)/(1+phi1*exp(phi2*agesExtrap)))
    }else{
      seq(obsQxt[nObs], to = 1, length = 30)
    }
  })
  
  qxt <- rbind(qxt,  extrapolateQxt)
  ages <- c(ages, agesExtrap)
  rownames(qxt) <- as.character(ages)
  list(qxt = qxt, ages.fit = ages)
}
