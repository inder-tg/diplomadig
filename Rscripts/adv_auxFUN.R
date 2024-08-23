

# ===============================
# --- core_FPCA_tryCatch STEP ---
# ===============================

dbind_test <- function (M1, M2) {
  if (is.null(M1)) {
    p <- M2
  }
  else {
    if (is.null(M2)) {
      p <- M1
    }
    else {
      r1 <- dim(M1)[1]
      r2 <- dim(M2)[1]
      c1 <- dim(M1)[2]
      c2 <- dim(M2)[2]
      r <- r1 + r2
      c <- c1 + c2
      p <- matrix(rep(0, r * c), ncol = c)
      p[(1:r1), (1:c1)] <- M1
      p[((r1 + 1):r), ((c1 + 1):c)] <- M2
    }
  }
  
  return(p)
}

spectralnorm_test <- function (M) {
  result <- sqrt(max(svd(t(M) %*% M)$d))
  return(result)
}

core_FPCA_tryCatch_test <- function(DATA, basis, k, corr, num_times){
  
  q <- c(2, 2)
  q.trend <- q[1]
  q.ef <- q[2]
  iterations <- ncol(DATA) * 10 # num.cases * 10
  convergence <- NULL
  MoD.all <- NULL
  ti <- NULL
  di <- NULL
  x <- seq(0, 1, length = nrow(DATA))
  nd <- floor(nrow(DATA)/8)  # floor(num_times/8)
  f.d1 <- matrix(nrow=nrow(DATA), ncol=ncol(DATA))
  f.d2 <- matrix(nrow=nrow(DATA), ncol=ncol(DATA))
  
  n <- num_times # nrow(DATA)
  message("Running STEP 1 * Computing Initial Values")
  q.trend1 <- 2 # why using 2?
  
  f.X <- basis$eigenvectorsQR[, 1:q.trend1]
  f.Z <- basis$eigenvectorsQR[, (q.trend1 + 1):nrow(DATA)] %*% 
    diag(1/sqrt(basis$eigenvalues[(q.trend1 + 1):nrow(DATA)]))
  f.D <- diag(basis$eigenvalues)
  f.N <- basis$eigenvectorsQR
  
  One.Vec.All <- rep(1, nrow(DATA) * ncol(DATA))
  X.pop <- kronecker(matrix(1, nrow = ncol(DATA)), f.X)
  Z.pop <- kronecker(matrix(1, nrow = ncol(DATA)), f.Z)
  
  y <- c(as.matrix(DATA))
  
  est <- nlme::lme(y ~ -1 + X.pop,
                   random = list(One.Vec.All = nlme::pdIdent(~Z.pop - 1)),
                   correlation = NULL)
  
  sigma2.e <- est$sigma^2
  sigma2.b.pop <- (est$sigma^2) * exp(2 * unlist(est$modelStruct))
  lambda <- sigma2.e/sigma2.b.pop
  
  f <- matrix(rep(est$fitted[1:nrow(DATA), 2], ncol(DATA)), 
              nrow=nrow(DATA), ncol=ncol(DATA))
  di <- DATA - f
  
  message("Running STEP 2 * Selecting k (Omited)")
  svd <- svd(t(di), nu = ncol(DATA), nv = nrow(DATA))
  u <- svd$u
  v <- svd$v
  d <- svd$d
  y <- t(di) %*% v
  ef <- v[, 1:k]
  if (k == 1 & (min(ef) < 0)) {
    ef <- -ef
  }
  
  message("Running STEP 3 * Updating Components")
  One.Vec.All <- rep(1, nrow(DATA) * ncol(DATA))
  time <- rep(x, ncol(DATA)) # num.cases)
  casesId <- rep(1:ncol(DATA), each = nrow(DATA))
  y <- c(as.matrix(DATA))
  
  # basis <- basis
  # f.X <- basis$eigenvectorsQR[, 1:q.trend]
  # f.Z <- basis$eigenvectorsQR[, (q.trend + 1):nrow(DATA)] %*% 
  #   diag(1/sqrt(basis$eigenvalues[(q.trend + 1):nrow(DATA)]))
  # f.D <- diag(basis$eigenvalues)
  # f.N <- basis$eigenvectorsQR
  
  f.X <- basis$eigenvectorsQR[, 1:q.trend]
  f.Z <- basis$eigenvectorsQR[, (q.trend + 1):nrow(DATA)] %*% 
    diag(1/sqrt(basis$eigenvalues[(q.trend + 1):nrow(DATA)]))
  f.D <- diag(basis$eigenvalues)
  f.N <- basis$eigenvectorsQR
  
  X.pop <- kronecker(matrix(1, nrow = ncol(DATA)), f.X)
  Z.pop <- kronecker(matrix(1, nrow = ncol(DATA)), f.Z)
  
  for (ITER in 1:iterations) {
    message(paste("Iteration ", ITER))
    f.Z.ef <- ef
    Z.cases.ef <- kronecker(matrix(1, nrow = ncol(DATA)), f.Z.ef)
    
    est <- nlme::lme(y ~ -1 + X.pop, 
                     random = list(One.Vec.All = nlme::pdIdent(~Z.pop - 1), 
                                   casesId = nlme::pdDiag(~Z.cases.ef - 1)), 
                     correlation = NULL)
    
    sigma2.e <- est$sigma^2
    sigma2.b.pop <- (est$sigma^2) * 
      exp(2 * unlist(est$modelStruct)[k + 1])
    lambda <- sigma2.e/sigma2.b.pop
    alpha <- t(est$coefficients$random$casesId)
    fi <- matrix(est$fitted[, 3], nrow=nrow(DATA), ncol=ncol(DATA))
    di <- ef %*% alpha
    f <- fi - di
    si <- NULL
    ei <- DATA - f
    ef.c <- NULL
    ef.d1.c <- NULL
    ef.d2.c <- NULL
    thetaf.c <- NULL
    if (k == 1) {
      alpha <- matrix(alpha, ncol = 1)
    } else {
      alpha <- matrix(alpha, ncol = ncol(alpha))
      alpha <- t(alpha)
    }
    
    for (i in 1:k) {
      if (k == 1) {
        di.m <- 0
      } else {
        di.m <- ef[, -i] %*% t(alpha[, -i])
      }
      r <- DATA - f - di.m
      r <- c(as.matrix(r))
      X.ef <- NULL
      Z.ef <- NULL
      N.ef <- NULL
      D.ef <- NULL
      
      for (j in 1:ncol(DATA)) {
        # basisd <- basis
        basisd <- basis
        alphamat <- diag(rep(alpha[j, i], nrow(DATA)))
        X.ef.c <- alphamat %*% basisd$eigenvectorsQR[, 1:q.ef]
        X.ef <- dbind_test(X.ef, X.ef.c)
        Z.ef.c <- alphamat %*% basisd$eigenvectorsQR[, (q.ef + 1):nd] %*% 
          diag(1/sqrt(basisd$eigenvalues[(q.ef + 1):nd]))
        
        Z.ef <- dbind_test(Z.ef, Z.ef.c)
        D.ef.c <- diag(basisd$eigenvalues[1:nd])
        D.ef <- dbind_test(D.ef, D.ef.c)
        N.ef.c <- alphamat %*% basisd$eigenvectorsQR[, 1:nd]
        N.ef <- dbind_test(N.ef, N.ef.c)
      }
      
      One.Vec.All <- rep(1, nrow(DATA) * ncol(DATA))
      if (is.null(corr)) {
        est.r <- nlme::lme(r ~ -1 + X.ef, 
                           random = list(One.Vec.All = nlme::pdIdent(~Z.ef - 1)), 
                           correlation = NULL)
      } else {
        est.r <- nlme::lme(r ~ -1 + X.ef, 
                           random = list(One.Vec.All = nlme::pdIdent(~Z.ef - 1)), 
                           correlation = NULL) # correlation=corr?
      }
      sigma2.e <- est.r$sigma^2
      sigma2.b <- (est.r$sigma^2) * exp(2 * unlist(est.r$modelStruct))
      lambda <- sigma2.e/sigma2.b
      add1 <- 0
      add2 <- 0
      X.ef <- basisd$eigenvectorsQR[, 1:q.ef]
      Z.ef <- basisd$eigenvectorsQR[, (q.ef + 1):nd] %*% 
        diag(1/sqrt(basisd$eigenvalues[(q.ef + 1):nd]))
      N.ef <- basisd$eigenvectorsQR[, 1:nd]
      D.ef <- diag(basisd$eigenvalues[1:nd])
      
      for (l in 1:ncol(DATA)) {
        r.c <- matrix(r, ncol = ncol(DATA))[, l]
        add1 <- add1 + (alpha[l, i]^2) * 
          diag(ncol(N.ef)) + (lambda/ncol(DATA)) * D.ef
        add2 <- add2 + (alpha[l, i]) * (t(N.ef) %*% r.c)
      }
      
      ef.c <- cbind(ef.c, N.ef %*% (solve(add1) %*% add2))
      thetaf.c <- cbind(thetaf.c, solve(add1) %*% add2)
    }
    
    convergence.c <- spectralnorm_test(ef - ef.c)
    convergence <- c(convergence, convergence.c)
    threshold <- 5/100
    
    if( convergence[length(convergence)] < (threshold) ){
      message("Convergence achieved")
      if (ITER == 1) {
        ef <- ef.c
        ef.d1 <- ef.d1.c
        ef.d2 <- ef.d2.c
        thetaf <- thetaf.c
        if (k > 1) {
          ef <- ef.c
          ef.d1 <- ef.d1.c
          ef.d2 <- ef.d2.c
          thetaf <- thetaf.c
        }
      }
      break
    } else {
      if (ITER == iterations) {
        stop("No convergence achieved")
      }
    }
    
    if (k == 1) {
      ef.c <- ef.c/sqrt(sum(ef.c^2))
    } else {
      thetaf.c <- qr(thetaf.c)
      thetaf.c <- qr.Q(thetaf.c)
      ef.c <- N.ef %*% thetaf.c
    }
    ef <- ef.c
    thetaf <- thetaf.c
  }
  
  message("Running STEP 4 * Getting Estimators for Components")
  if (is.null(corr)) {
    corcoefs <- "Model fitted with white noise remainder"
    est <- nlme::lme(y ~ -1 + X.pop, 
                     random = list(One.Vec.All = nlme::pdIdent(~Z.pop - 1), 
                                   casesId = nlme::pdDiag(~Z.cases.ef - 1)), 
                     correlation = NULL)
  } else {
    est <- nlme::lme(y ~ -1 + X.pop, 
                     random = list(One.Vec.All = nlme::pdIdent(~Z.pop - 1), 
                                   casesId = nlme::pdDiag(~Z.cases.ef - 1)), 
                     correlation = corr)
  }
  
  sigma2.e <- est$sigma^2
  sigma2.b.pop <- (est$sigma^2) * exp(2 * unlist(est$modelStruct)[k + 1])
  lambda <- sigma2.e/sigma2.b.pop
  cr.Zef <- est$coefficients$random$casesId
  alpha <- cr.Zef
  
  fi <- matrix(est$fitted[, 3], nrow=nrow(DATA), ncol=ncol(DATA))
  di <- ef %*% t(alpha)
  f <- fi - di
  error <- DATA - fi
  nf <- ncol(X.pop)
  nr <- ncol(Z.pop)
  beta <- est$coefficients$fixed[1:nf]
  mu <- est$coefficients$random$One.Vec.All[1:nr]
  samples <- t(MASS::mvrnorm(n = 10000, f.X %*% beta, sigma2.b * 
                               f.Z %*% t(f.Z)))
  upp <- apply(samples, 1, max)
  low <- apply(samples, 1, min)
  
  list(est = est, f = f[, 1])
}

# ==========================
# --- FPCA_tryCatch STEP ---
# ==========================

show_condition_test <- function(code) {
  tryCatch(code, error = function(c) "ERROR")
}

FPCA_tryCatch_test <- function(DATA, k, corr = NULL, basis) {
  
  if( missing(k) ){
    stop("k must be an integer between 1 and the sample size")
  }
  
  df_DATA <- data.frame(DATA)
  
  if( ncol(df_DATA) < 2 ){
    stop("Data matrix must contain at least 2 columns")
  }
  
  if( missing(basis) ){
    basis <- drbasis(nn = nrow(df_DATA), qq = 2)
    message("Running STEP 0 * Computing DR Basis")
  } else {
    message("Running STEP 0 * Loading DR Basis")
  }
  
  auxiliary <- show_condition_test(core_FPCA_tryCatch_test(DATA=DATA, basis=basis,
                                                           corr=corr, k=k,
                                                           num_times=nrow(DATA)))
  
  out <- list()
  if(inherits(auxiliary, "character")){
    out$status <- "Failure"
    out$est <- NA
    out$f <- NA
  } else {
    out$status <- "Success"
    out$est <- auxiliary$est
    out$f <- auxiliary$f
  }
  
  list(est=out$est, f=out$f, status=out$status)
}

# ====================
# --- getFPCA STEP ---
# ====================

getFPCA_test <- function(series, seriesByCluster, clusterSize, k, mAug, basis, corr){
  
  flag <- FALSE # 
  
  if(missing(series)) {
    seriesInCluster <- lapply(1:length(unique(seriesByCluster)), 
                              function(s) which(seriesByCluster == s))
    
    if( length(seriesInCluster[[1]]) >= clusterSize ) {
      SERIES <- seriesInCluster[[1]]
    }
    
    if( length(seriesInCluster[[2]]) >= clusterSize ) {
      SERIES <- seriesInCluster[[2]]
    }
    
    AUX <- ifelse(any(table(seriesByCluster) >= clusterSize), TRUE, FALSE)
    if( !AUX ) {
      SERIES <- 1:length(seriesByCluster)
      flag <- TRUE
    }
  } else {
    SERIES <- series
    flag <- ifelse( length(SERIES) > clusterSize, TRUE, FALSE )
  }
  
  # maug <- mAug[,SERIES]
  out_put <- FPCA_tryCatch_test(DATA=mAug[,SERIES], k=k, basis=basis, corr=corr)
  
  if( out_put$status == "Failure" ){
    out_put <- list()
    out_put$f <- NA
    flag <- NA
  }
  
  # --- usedTotal is not updated, must report how many series were used!
  
  list(fpca=out_put$f, usedTotal=flag, seriesUsed=SERIES, status=out_put$status)
}


# =============================
# --- getPhenopar_tryCatch ---
# =============================

check_phenoCond_test <- function(GU_Mat, SoS_EoS, Sen, Dor, f2der){
  
  x_GU <- GU_Mat$max
  
  x_SoS <- SoS_EoS$max
  
  x_Mat <- GU_Mat$min
  
  x_Sen <- Sen$x_opt # Sen$mins[ which.max(f2der(Sen$mins)) ]
  
  x_EoS <- SoS_EoS$min
  
  x_Dorm <- Dor$x_opt # Dorm$maxs[ which.min(f2der(Dorm$maxs)) ]
  
  params <- list(GU= ifelse(x_GU%%1 < 0.5, floor(x_GU), ceiling(x_GU)), 
                 SoS=ifelse(x_SoS%%1 < 0.5, floor(x_SoS), ceiling(x_SoS)), 
                 Mat=ifelse(x_Mat%%1 < 0.5, floor(x_Mat), ceiling(x_Mat)), 
                 Sen=ifelse(x_Sen%%1 < 0.5, floor(x_Sen), ceiling(x_Sen)), 
                 EoS=ifelse(x_EoS%%1 < 0.5, floor(x_EoS), ceiling(x_EoS)), 
                 Dor=ifelse(x_Dorm%%1 < 0.5, floor(x_Dorm), ceiling(x_Dorm)))
  
  num_of_na <- sum(is.na(unlist(params)))
  
  flag <- ifelse( num_of_na==0, "Full", 
                  ifelse(num_of_na==1, "Partial-5", 
                         ifelse(num_of_na==2, "Partial-4",
                                ifelse(num_of_na==3, "Partial-3",
                                       ifelse(num_of_na==2, "Partial-4",
                                              ifelse(num_of_na==1, "Partial-5", "Failure")
                                       )
                                ) 
                         ) 
                  ) 
  )
  
  list(params=params, flag=flag)
}


getPhenoDatesUpdate_test <- function(fit, numFreq=3, delta=0, L=365, 
                                     interval=seq(1,365,length=365)){
  
  fzero <- ndvi_derivatives_test(amp=fit$amplitude, pha=fit$phase, degree=0, L=L)
  
  fprime <- ndvi_derivatives_test(amp=fit$amplitude, pha=fit$phase, degree=1, L=L)
  
  fbiprime <- ndvi_derivatives_test(amp=fit$amplitude, pha=fit$phase, degree=2, L=L)
  
  fthrprime <- ndvi_derivatives_test(amp=fit$amplitude, pha=fit$phase, degree=3, L=L)
  
  ffthprime <- ndvi_derivatives_test(amp=fit$amplitude, pha=fit$phase, degree=4, L=L)
  
  GU_Mat <- global_min_max_test(f=fbiprime, f1der=fthrprime, f2der=ffthprime, D=interval)
  
  # Sen <- get_xsen(f=fbiprime, f1der=fthrprime, f2der=ffthprime, x0=GU_Mat$min)
  Sen <- local_min_max_test(f=fbiprime, f1der=fthrprime, f2der=ffthprime, what="min", x0=GU_Mat$min, 
                            D=interval)
  
  SoS_EoS <- global_min_max_test(f=fprime, f1der=fbiprime, f2der=fthrprime, D=interval)
  
  # Dorm <- get_xdorm(f=fbiprime, f1der=fthrprime, f2der=ffthprime, x0=GU_Mat$max)
  Dorm <- local_min_max_test(f=fbiprime, f1der=fthrprime, f2der=ffthprime, what="max", 
                             x0=GU_Mat$max, D=interval)
  
  parPheno <- check_phenoCond_test(SoS_EoS=SoS_EoS, GU_Mat=GU_Mat, Sen=Sen, Dor=Dorm,
                                   f2der=fbiprime)
  
  list(dates=unlist(parPheno$params),
       status_phenodates=parPheno$flag, 
       zder=fzero, fder=fprime, sder=fbiprime, 
       thrder=fthrprime, fthder=ffthprime)
}

getPhenopar_tryCatch_test <- function(x, frequency=23, method=c("OLS", "WLS"),
                                       sigma, numFreq, delta=0, 
                                       distance, clusterSize,
                                       samples, basis, corr, k, trace){
  
  # x=test
  # frequency=12; method="OLS"
  # sigma=NULL; numFreq=3; delta=0
  # distance="dtw2"; clusterSize=6
  # samples=100; basis=BASIS_fourthRound
  # corr=NULL; k=3; trace=FALSE
  
  method <- match.arg(method)
  
  status <- "Failure"
  
  xMat <- vecToMatrix_test(x=x, lenPeriod=frequency)
  
  if(method == "OLS"){
    Sigma <- sigma
  } else {
    Sigma <- hetervar(m=xMat, lenPeriod = frequency, method = "standard")
  }
  
  mat_smooth_sigma <- getSmoothing_test(m=xMat, sigma=Sigma, numFreq=numFreq, 
                                        delta=delta, method=method)
  
  clust_distance <- tsclusters_centroid_test(m=mat_smooth_sigma, seed=11,
                                             distance=distance, trace=trace)
  
  mat_augSmooth <- getHarmonicFit_test(samples=samples, m=mat_smooth_sigma,
                                       sigma=Sigma, numFreq=numFreq,
                                       delta=delta, method=method)
  
  x_fpca <- getFPCA_test(seriesByCluster=clust_distance@cluster, 
                         clusterSize=clusterSize, k=k,
                         mAug=mat_augSmooth, basis=basis, corr=corr)
  
  if(x_fpca$status=="Failure"){
    fit <- list()
    fit$a.coef <- NA
    fit$b.coef <- NA
    fit$amplitude <- NA
    fit$phase <- NA
    
    phenoDates <- list()
    phenoDates$zder <- NA
    phenoDates$fder <- NA
    phenoDates$sder <- NA
    phenoDates$thrder <- NA
    phenoDates$fthder <- NA
    
    phenoDates$status_phenodates <- "Failure" #status #"Failure"
    
    message("FPCA model non inverted, try increasing 'samples'") # MAX iterations was achieved
  } else {
    fit <- haRmonics(y=x_fpca$fpca, numFreq=numFreq, delta=delta)
    
    # --- this seems to be another source of issues
    phenoDates <- getPhenoDatesUpdate_test(fit=fit, numFreq=numFreq)
    
  }
  
  if(phenoDates$status_phenodates == "Full"){
    status <- "Success"
  }
  
  # if(phenoDates$status_phenodates == "Partial"){
  #   status <- "Partial"
  # }
  
  if(unlist(strsplit(phenoDates$status_phenodates, "-"))[1] == "Partial"){
    status <- "Partial"
  }
  
  list(x=x,
       freq=frequency,
       sigma=Sigma, 
       x_smooth=mat_smooth_sigma,
       m_aug_smooth=t(mat_augSmooth),
       clustering=clust_distance,
       clusterSize=clusterSize,
       fpca=x_fpca$fpca, 
       fpca_harmfit_params=list(a.coef=fit$a.coef, b.coef=fit$b.coef, 
                                amplitude=fit$amplitude, phase=fit$phase),
       fpca_fun_0der=phenoDates$zder,
       fpca_fun_1der=phenoDates$fder,
       fpca_fun_2der=phenoDates$sder,
       fpca_fun_3der=phenoDates$thrder,
       fpca_fun_4der=phenoDates$fthder,
       phenoparams=phenoDates$dates,
       status=status,
       series=x_fpca$seriesUsed)
}
# ---

# =================================
# --- phenoparPolygon_auxiliar ---
# =================================

save_output <- function(y,frequency,method,sigma,trace,numFreq,delta,distance,clusterSize,
                        samples,basis,k,corr,objectsToSave,fileToSaveErrors,dirToSaveErrors,i){
  
  get_phenoparams <- show_condition_test( getPhenopar_tryCatch_test(x=y, 
                                                                    frequency=frequency, 
                                                                    method=method,
                                                                    sigma=sigma, 
                                                                    trace=trace,
                                                                    numFreq=numFreq, 
                                                                    delta=delta, 
                                                                    distance=distance, 
                                                                    clusterSize=clusterSize,
                                                                    samples=samples, 
                                                                    basis=basis, 
                                                                    k=k, 
                                                                    corr=corr) )
  
  if(inherits(get_phenoparams, "character")){
    
    texto <- paste0("Found ERROR on ROW: ", i)
    write(texto, file = fileToSaveErrors, append = TRUE)
    
    s <- "Error" 
    
    # errorsFromFile <- LoadToEnvironment(paste0( dirToSaveErrors, "/phenopar_polygon_errors.RData" ))$errorObject
    # 
    # errorObject <- c(errorsFromFile, i)
    # 
    # # write(i, file=errorsFile, append = TRUE)
    # 
    # save(errorObject, file = errorsFromFile)
    
  } else {
    s <- what_objects_test(x=get_phenoparams,
                           whatToExtract=objectsToSave)
  }
  
  s
}

phenoparPolygon_auxiliar_test <- function(data, product, frequency, method, sigma,
                                           numFreq, delta, distance, clusterSize, samples,
                                           basis, corr, k, numCores, 
                                           objectsToSave, dirToSave,
                                           reportFileBaseName, trace) {
  
  reportFileName <- paste0( dirToSave, "/", reportFileBaseName, ".txt" )
  file.create(path = reportFileName, showWarnings = FALSE)
  write( "--- SEPHORA began at ---", file = reportFileName,
         append = TRUE )
  write( as.character(Sys.time()[1]), file = reportFileName,
         append = TRUE )
  
  # --- Added on May 31, 2024
  reportErrorsFileName <- paste0( dirToSave, "/", reportFileBaseName, "_errors.txt" )
  file.create(path = reportErrorsFileName, showWarnings = FALSE)
  write( "--- SEPHORA began at ---", file = reportErrorsFileName,
         append = TRUE )
  write( as.character(Sys.time()[1]), file = reportErrorsFileName,
         append = TRUE )
  
  # errorObject <- c()
  # # errorObject <- c(errorObject,0)
  # errorsFromFile <- paste0( dirToSave, "/phenopar_polygon_errors.RData" )
  # # paste0( dirToSaveErrors, "/phenopar_polygon_errors.RData" )
  # save(errorObject, file = errorsFromFile)
  # # file.create(path = errorsFromFile, showWarnings = FALSE)

  # ---
  
  dataPOLYGON <- data 
  
  kluster <- parallel::makeCluster(spec=numCores, outfile="")
  registerDoParallel(kluster)
  
  output <- foreach(i=1:nrow(dataPOLYGON), #.combine="rbind",
                    .export=c("vecFromData_test", "vecToMatrix_test",
                              "fill_initialgap_MOD13Q1_test",
                              "getPhenopar_tryCatch_test",
                              "hetervar", "getSmoothing_test",
                              "tsclusters_centroid_test",
                              "getHarmonicFit_test",
                              "getFPCA_test", "FPCA_tryCatch_test",
                              "core_FPCA_tryCatch_test",
                              "show_condition_test", "getPhenoDatesUpdate_test",
                              "harmonicFit_test", "fitHarmonic_test",
                              "spectralnorm_test", "dbind_test",
                              "global_min_max_test", "check_phenoCond_test",
                              "local_min_max_test", "ndvi_derivatives_test",
                              "what_objects_test", "save_output"),
                    .packages=c("dtwclust", "TSclust", "geoTS", "bigmemory",
                                "rootSolve")) %dopar% {
                                  
                                  x <- vecFromData_test(data = dataPOLYGON, numRow = i, 
                                                        product = product, lenPeriod=frequency)
                                  
                                  # objectsToSave,fileToSaveErrors,dirToSaveErrors,i,errorObject)
                                  
                                  s <- save_output(y=x$vec, 
                                                   frequency=frequency, 
                                                   method=method,
                                                   sigma=sigma, 
                                                   trace=trace,
                                                   numFreq=numFreq, 
                                                   delta=delta, 
                                                   distance=distance, 
                                                   clusterSize=clusterSize,
                                                   samples=samples, 
                                                   basis=basis, 
                                                   k=k, 
                                                   corr=corr,
                                                   objectsToSave=objectsToSave,
                                                   fileToSaveErrors=reportErrorsFileName,
                                                   dirToSaveErrors=dirToSave,
                                                   i=i) #,
                                                   # errorObject=errorsObj)
                                  
                                  chunkSize <- floor(nrow(dataPOLYGON) * 0.2)
                                  
                                  if( i %% chunkSize == 0 ){
                                    texto <- paste0("Working on ROW: ", i)
                                    write(texto, file = reportFileName, append = TRUE)
                                  }
                                  
                                  return(s)
                                }
  stopCluster(kluster)
  
  write( as.character(Sys.time()[1]), file = reportFileName, append = TRUE)
  write(paste("OUTPUT will be saved at ", dirToSave,  " shortly"),
        file = reportFileName, append = TRUE)
  write( "--- SEPHORA ended at ---", file = reportFileName, append = TRUE)
  write( as.character(Sys.time()[1]), file = reportFileName, append = TRUE)
  
  write( "--- SEPHORA ended at ---", file = reportErrorsFileName, append = TRUE)
  write( as.character(Sys.time()[1]), file = reportErrorsFileName, append = TRUE)
  
  output
}

# =============================
# --- phenopar_polygon STEP ---
# =============================

phenopar_polygon_test <- function(path=NULL, product=c("MOD13Q1", "independent"), data,
                                  frequency=23, method=c("OLS", "WLS"), 
                                  save_all=FALSE, save_some, 
                                  sigma=NULL, numFreq, delta=0, distance, clusterSize, 
                                  basis, samples, corr=NULL,
                                  k, trace=FALSE, numCores=20, dirToSave,
                                  reportFileBaseName="phenopar_progress",
                                  outputFileBaseName="polygon"){
  
  # product = "independent"
  # data = DATA_rTp_available_fourthRound$values
  # frequency = 12
  # method="OLS"
  # save_all=FALSE
  # sigma=NULL
  # numFreq = 1
  # distance = "dtw2"
  # clusterSize = 6
  # basis = BASIS_fourthRound
  # samples=125
  # k = 3
  # numCores = 3
  # # save_all = TRUE,
  # save_some = c("clustering", "fpca", "phenoparams", "status")
  # reportFileBaseName="phenopar_progress_fourthRound"
  # outputFileBaseName="polygon_phenoParams_fourthRound"
  # dirToSave = paste0( getwd(), "/RData/chunk_1_2/output_numFreq_2_k_3" )

  
  product <- match.arg(product)
  if(product == "MOD13Q1"){
    if(is.null(path)){
      stop('path must be provided')
    }
    
    dataEnv <- LoadToEnvironment(RData=path)
    data <- get(x="poly", envir=dataEnv)[[1]]
    
  } else {
    if(missing(data)){
      stop('data must be provided')
    }
  }
  
  method <- match.arg(method)
  
  if(save_all){
    objectsToSave <- validObjects
  } else {
    if(missing(save_some)){
      stop("save_some must provided")
    } else {
      toExtract <- base::intersect(x=validObjects, y=save_some)
      objectsToSave <- toExtract
    }
  }
  
  if(missing(numFreq)){
    stop("numFreq must provided")
  }
  
  if(missing(distance)){
    stop("distance must provided")
  }
  
  if(missing(clusterSize)){
    stop("clusterSize must provided")
  }
  
  if(missing(basis)){
    stop("basis must provided")
  }
  
  if(missing(samples)){
    if(missing(basis)){
      stop("samples must provided")
    } else {
      samples <- length(basis[["x"]])
    }
  } else {
    if(samples != length(basis[["x"]])){
      stop("Length of samples and length(basis$x) must be equal")
    }
  }
  
  if(missing(k)){
    stop("k must be provided")
  }
  
  if(missing(dirToSave)){
    stop('dirToSave must be provided')
  }
  
  outPut <- phenoparPolygon_auxiliar_test(data=data, product=product, frequency=frequency, 
                                           method=method, sigma=sigma, numFreq=numFreq, delta=delta, 
                                           distance=distance, clusterSize=clusterSize,
                                           samples=samples, basis=basis, corr=corr, k=k, 
                                           trace=trace, numCores=numCores, objectsToSave=objectsToSave,
                                           dirToSave=dirToSave, 
                                           reportFileBaseName=reportFileBaseName)
  
  save(outPut, file = paste0(dirToSave, "/",
                             ifelse( !is.null(path), 
                                     tools::file_path_sans_ext(basename(path)), 
                                     outputFileBaseName ),
                             "_phenoParams.RData"))
  
}
