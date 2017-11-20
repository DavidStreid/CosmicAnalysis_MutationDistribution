# ADD RELATIVE PATH
setwd('PATH/TO/CosmicAnalysis_MutationDistribution/analysis')

# Plots observed data, @sample: name of sample (E.g. ARID1A)
plot_data <- function(sample){
  # Reading in Data
  df <- read.table(sprintf("../data/cleanedData/%s.csv", sample), 
                   header = FALSE,
                   sep = ",")
  smooth_x <- unlist(df['V1'])
  smooth_y <- unlist(df['V2'])
  yLim = c(0,max(smooth_y))
  plot(smooth_x,smooth_y,ylim=yLim,
       ylab="Frequency", xlab="Mutation Positions",
       main=sprintf("%s Mutation Position versus Frequency", sample))
}

# Performs constant-span running mean 
# @y: observed y values, @k: span size
csrm <- function(y,k){
  est_y <- vector(mode="numeric", length=0);
  l <- ((k-1)/2);
  n <- length(y);
  for(i in 1:l){
    est_y <- c(est_y,mean(y[1:(i+l)]))
  }
  for(i in (l+1):(n-l)){
    est_y <- c(est_y,mean(y[(i-l):(i+l)]))
  }
  for(i in (n-l+1):n){
    est_y <- c(est_y,mean(y[(i - l):n]))
  }
  return(est_y);
}

# Calculates CVRSS value for observed & estimated y values
# @obs_y: observed y values, @est_y: estimated y values, @k: size of span
cvrss <- function(obs_y, est_y, k){
  n <- length(obs_y);
  l <- ((k-1)/2);
  cvrss <- 0;
  if(n != length(est_y)){
    print("Mismatching dimensions");
    return(-1);
  }
  
  for(i in 1:l){
    S_ii <- 1/(i+l);
    cvrss <- cvrss + ((obs_y[i] - est_y[i])/(1-S_ii))**2;
  }
  S_ii <- 1/k;
  for(i in (l+1):(n-l)){
    cvrss <- cvrss + ((obs_y[i] - est_y[i])/(1-S_ii))**2;
  }
  for(i in (n-l+1):n){
    S_ii <- 1/(n - i + 1 + l);
    cvrss <- cvrss + ((obs_y[i] - est_y[i])/(1-S_ii))**2;
  }
  
  return(cvrss/n);
}

# Plots the smoother with optimal span & returns optimal k
# @obs_y: Observed response variables
# @mean: boolean (T -> Mean, F -> Median)
plot_cvrss <- function(obs_y,mean){
  cvrss_vals <- vector(mode="numeric", length=0);
  k_vals <- vector(mode="numeric", length=0);
  min_cvrss <- 9999;
  min_k <- -1;
  for(i in 1:11){
    k <- (2*i + 1);
    if(mean){
      est_y <- csrm(obs_y,k);  
    }
    else{
      est_y <- csrmed(obs_y,k);
    }
    cvrss_i <- cvrss(obs_y, est_y, k);
    
    k_vals <- c(k_vals,k);
    cvrss_vals <- c(cvrss_vals,cvrss_i);
    
    if(cvrss_i < min_cvrss){
      min_cvrss <- cvrss_i;
      min_k <- k;
    }
  }
  
  # Output best k
  print(min_cvrss);
  print(min_k);
  plot(k_vals,cvrss_vals,ylab='spans',main='CVRSS v. k spans');
  return(min_k);
}

plot_csrm <- function(sample, overlay){
  # Reading in Data
  df <- read.table(sprintf("../data/cleanedData/%s.csv", sample), 
                   header = FALSE,
                   sep = ",")
  smooth_x <- unlist(df['V1'])
  smooth_y <- unlist(df['V2'])
  
  # Optimizing Constant-Span Running Mean
  opt_k <- plot_cvrss(smooth_y,TRUE) 
  est_y <- csrm(smooth_y,opt_k)
  yLim = c(0,max(smooth_y))
  
  plot(smooth_x,est_y,type="l",ylim=yLim,
       ylab="predicted y", xlab="x",
       main=sprintf("%s CSRM Smoother, k=%d", sample, opt_k))
  par(new=TRUE)
  if(overlay){
    plot(smooth_x,smooth_y,ylim=yLim,
         ylab="", xlab="",
         main="")
  }
}

plot_data('ARID1A')
plot_csrm('ARID1A',FALSE)
plot_csrm('ARID1A',TRUE)

plot_data('ARID1B')
plot_csrm('ARID1B',FALSE)
plot_csrm('ARID1B',TRUE)

plot_data('ARID2')
plot_csrm('ARID2',FALSE)
plot_csrm('ARID2',TRUE)
