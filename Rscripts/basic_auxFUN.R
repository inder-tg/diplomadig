

# --- March, 18, 2022

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}

# --- Added on Feb 24, 2024

spRast_valueCoords <- function(spRaster, na_rm=FALSE){
  
  spPoints <- as.points(spRaster, na.rm=na_rm)
  
  spValues <- terra::extract(spRaster, spPoints)
  
  DIM <- dim(spValues)
  
  spRasterToPoints <- as.matrix(spValues[1:DIM[1],2:DIM[2]])
  
  spCoords <- crds(spRaster, na.rm=na_rm)
  
  list(values=spRasterToPoints, coords=spCoords)  
}


getSmoothing_test <- function(m, method=c("OLS", "WLS"), numFreq=3, delta=0, 
                               sigma=NULL){
  
  method <- match.arg(method)
  
  if(method=="OLS"){
    method <- "harmR"
  } else {
    method <- "wls_harmR"
    if(is.null(sigma)){
      stop("sigma must be provided when WLS method is used")
    }
  }
  
  smooth_m <- matrix(nrow = nrow(m), ncol = ncol(m))
  
  for(i in 1:nrow(m)){
    smooth_m[i,] <- haRmonics(y=m[i,], method=method, sigma=sigma, 
                              numFreq=numFreq, delta=delta)$fitted
  }
  
  smooth_m
}


fitHarmonic_test <- function(amp,pha,L,t){
  sum(sapply( 0:(length(amp)-1),
              function(k) amp[k+1] * cos( 2 * pi * k * t / L - pha[k+1]/(180/pi) )  
  )
  )
}

getHarmonicFit_test <- function(samples=100, m, method=c("OLS", "WLS"), sigma=NULL, 
                                 numFreq=4, delta=0){

  method <- match.arg(method)
  
  if(method == "OLS"){
    # method <- "ols_harmR"
    method <- "harmR"
  } else {
    method <- "wls_harmR"
    if(is.null(sigma)){
      stop("sigma must be provided when WLS method is used")
    }
  }
  
  if ( length(numFreq) == 1 ){
    numFreq <- rep(numFreq, nrow(m))
  }
  
  smooth_m_aug <- matrix(nrow=samples, ncol=nrow(m))
  for(i in 1:ncol(smooth_m_aug)){
    TEMP <- haRmonics(y=m[i,], method=method, sigma=sigma, numFreq=numFreq[i], 
                      delta=delta)
    
    smooth_m_aug[,i] <- sapply(seq(0, (ncol(m)-1), length=samples),
                               function(s) fitHarmonic_test(amp=TEMP$amplitude,
                                                            pha=TEMP$phase,
                                                            t=s, L=ncol(m)))
  }
  
  smooth_m_aug
}


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

# --- Added May 29, 2023

validObjects <- c("x", "freq", "sigma", "x_smooth",
                  "m_aug_smooth", "clustering",
                  "clusterSize", "fpca",
                  "fpca_harmfit_params",
                  "fpca_fun_0der", "fpca_fun_1der",
                  "fpca_fun_2der", "fpca_fun_3der",
                  "fpca_fun_4der", "phenoparams",
                  "status", "series")

vecFromData_test <- function(product=c("MOD13Q1", "independent"), data, 
                             numRow, lenPeriod=23){
  
  data_vec <- as.numeric(data[numRow,])
  data_mat <- vecToMatrix_test(x=data_vec, lenPeriod=lenPeriod)
  
  product <- match.arg(product)
  
  if( product == "MOD13Q1" ){
    data_mat[1,1:3] <- fill_initialgap_MOD13Q1_test(m=data_mat)
  }
  
  list(mat=data_mat, vec=c(t(data_mat)))

}


vecToMatrix_test <- function(x, lenPeriod=23){
  
  if(length(x) %% lenPeriod !=0){
    stop("Length of 'x' must be a multiple of 'lenPeriod'")
  }
  
  output <- matrix(nrow=length(x)/lenPeriod, ncol=lenPeriod)
  
  for(i in seq_len(nrow(output))){
    output[i,] <- x[((i-1) * lenPeriod + 1):(i * lenPeriod)]
  }
  
  output
}

fill_initialgap_MOD13Q1_test <- function(m, fun=stats::median){
  apply( m[-1,1:3], MARGIN=2, FUN=fun)
}

tsclusters_centroid_test <- function(m, type='h', seed, distance="L2",
                                     trace){
  
  tsclust(series=m, type = type, seed = seed, distance = distance, 
          centroid = dtwclust::shape_extraction, trace = trace, 
          control = dtwclust::hierarchical_control(method = "average"))
}

fitHarmonic_test <- function(amp,pha,L,t){
  sum(sapply( 0:(length(amp)-1),
              function(k) amp[k+1] * cos( 2 * pi * k * t / L - pha[k+1]/(180/pi) )  
  )
  )
}

harmonicFit_test <- function(amp,pha,L,t){
  sapply(1:length(t), function(s) fitHarmonic_test(amp=amp, pha=pha, L=L, t=s) )
}

global_min_max_test <- function(f, f1der, f2der, D){
  
  crtPts <- uniroot.all(f=f1der, interval= c(D[1],D[length(D)]),
                        tol = .Machine$double.eps)
  
  maxPts <- which(f2der(crtPts)<0)
  
  minPts <- which(f2der(crtPts)>0)
  
  MAX <- crtPts[maxPts[which.max(f(crtPts[maxPts]))]]
  
  MIN <- crtPts[minPts[which.min(f(crtPts[minPts]))]]
  
  list(min=MIN, max=MAX, mins=crtPts[minPts], maxs=crtPts[maxPts])
}

local_min_max_test <- function(f, f1der, f2der, what=c("min", "max"), x0, D){
  if(x0 < D[1] | x0 > D[length(D)]){
    stop("x0 must satisfy that D[1] < x0 < D[length(D)]")
  }
  
  what <- match.arg(what)
  
  x0 <- ifelse(x0%%1 < 0.5, floor(x0), ceiling(x0))
  
  x_int1 <- uniroot.all(f=f1der, interval = c(D[1], x0-1))
  x_int2 <- uniroot.all(f=f1der, interval = c(x0+1, D[length(D)]))
  
  if(what == "min"){
    x_int1_valid <- which(f2der(x_int1)>0)
    x_int2_valid <- which(f2der(x_int2)>0)
    
    valid_x_int1 <- x_int1[x_int1_valid[which.min(f(x_int1[x_int1_valid]))]]
    valid_x_int2 <- x_int2[x_int2_valid[which.min(f(x_int2[x_int2_valid]))]]
    
    if( (length(valid_x_int1) != 0) & (length(valid_x_int2) != 0) ){
      
      xopt <- c(valid_x_int1, valid_x_int2)[which.min(c(f(valid_x_int1), f(valid_x_int2)))]
      
    } else {
      if( length(valid_x_int1) == 0 ){
        xopt <- ifelse( f(valid_x_int2) > f(x0), valid_x_int2, NA )
      } else {
        xopt <- ifelse( f(valid_x_int1) > f(x0), valid_x_int1, NA )
      }
    }
    
  } else {
    x_int1_valid <- which(f2der(x_int1)<0)
    x_int2_valid <- which(f2der(x_int2)<0)
    
    valid_x_int1 <- x_int1[x_int1_valid[which.max(f(x_int1[x_int1_valid]))]]
    valid_x_int2 <- x_int2[x_int2_valid[which.max(f(x_int2[x_int2_valid]))]]
    
    if( (length(valid_x_int1) != 0) & (length(valid_x_int2) != 0) ){
      
      xopt <- c(valid_x_int1, valid_x_int2)[which.max(c(f(valid_x_int1), f(valid_x_int2)))]
      
    } else {
      if( length(valid_x_int1) == 0 ){
        xopt <- ifelse( f(valid_x_int2) < f(x0), valid_x_int2, NA )
      } else {
        xopt <- ifelse( f(valid_x_int1) < f(x0), valid_x_int1, NA )
      }
    }
    
  }
  
  list(x_opt=xopt, locals=c(valid_x_int1, valid_x_int2),
       crtPts=list(x_d1=x_int1[x_int1_valid], x_d2=x_int2[x_int2_valid]),
       type=what)
  
}

ndvi_derivatives_test <- function(amp, pha, degree, L){
  
  vp <- 0:(length(amp)-1)
  
  if( degree == 0 ){
    core_derivative <- function(t) cos( vp * (2*pi/L) * t - pha/(180/pi) )
  } else {
    if( degree %% 2 == 0 ){
      core_derivative <- function(t) (-1)^(degree/2) * (2*pi/L)^(degree/2) * vp^degree *
        cos( vp * (2*pi/L) * t - pha/(180/pi) )
    } else {
      core_derivative <- function(t) (-1)^( (degree-1)/2 + 1) * (2*pi/L)^degree *
        vp^degree * sin( vp * (2*pi/L) * t - pha/(180/pi) )
    }
  }
  
  function(t) sapply( t, function(s)  sum( amp * core_derivative( s ) ) )
}

what_objects_test <- function(x, whatToExtract){
  
  output <- vector(mode="list", length=length(whatToExtract))
  output <- lapply(whatToExtract, function(s) x[[s]] )
  names(output) <- whatToExtract
  
  output
}

# Added on June 7, 2024

getErrorIndeces_phenoparPolygon <- function(fileName){
  my_data <- read.delim(file = fileName, header = FALSE, skip=2)
  
  out <- unlist(lapply(1:(nrow(my_data)-2), 
                       function(s){ temp <- strsplit(my_data[s,], " ")[[1]]; as.numeric(temp[length(temp)])} ))
  
  out
}

# Added on Aug 23, 2024

extractPhenoParams <- function(data, L){
  
  out <- vector("list", L)
  
  badIndices <- c()
  
  for(i in 1:length(data) ){
    if( length(data[[i]]) == 3 ){
      out[[i]] <- data[[i]]$phenoparams
    } else {
      badIndices <- c(badIndices, i)
    }
  }
  
  list(phenoParams=out, badIndices=badIndices)
}

getParam <- function(x, phenoparam){
  
  sapply(1:length(x), function(s) x[[s]][ names( x[[s]] )  == phenoparam ] )
  
}

