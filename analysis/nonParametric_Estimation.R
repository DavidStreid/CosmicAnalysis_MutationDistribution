# ADD RELATIVE PATH
setwd('PATH/TO/CosmicAnalysis_MutationDistribution/analysis')

# UTILITY FUNCTIONS
# read_data: Reads in data
# @sample: ARID1A, ARID1B, ARID2; @filled: boolean - filled data or not
read_data <- function(sample, filled){
  if(filled){
    return(read.table(sprintf("../data/cleanedData/%s_filled.csv", sample), 
                      header = FALSE,
                      sep = ","))
  }
  else{
    return(read.table(sprintf("../data/cleanedData/%s.csv", sample), 
                      header = FALSE,
                      sep = ","))
  }
}
# make_hist: Makes histogram with set number of bins
# @data: Input data, @nBins: Number of bins
make_hist <- function(data, sample, nBins){
  minVal <- min(data);
  maxVal <- max(data);
  print(maxVal);
  interval <- ((maxVal-minVal)/nBins);
  hist(data, breaks=seq(minVal,maxVal,interval),
       main=sprintf("Histogram - %s", sample))
}

# CSRM Functions
# Plots observed data, @sample: name of sample (E.g. ARID1A)
# @sample: sample data to be plotted, @filled: Missing data points filled in (with 0)
plot_data <- function(sample, filled){
  # Reading in Data
  df <- read_data(sample, filled)
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

# Plots CVRSS to determine optimal k value
# @obs_y: Observed response variables
# @mean: boolean (T -> Mean, F -> Median) when calculating CSRM
# @plot_graph: boolean (T -> Plot, F -> Don't Plot)
plot_cvrss <- function(obs_y,mean,plot_graph){
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
  if(plot_graph){
    plot(k_vals,cvrss_vals,ylab='spans',main='CVRSS v. k spans');  
  }
  return(min_k);
}

# Plots the smoother with optimal span & returns optimal k
# @obs_y: Observed response variables
# @mean: boolean (T -> Mean, F -> Median)
plot_csrm <- function(sample, overlay, filled){
  # Reading in Data
  df <- read_data(sample,filled)
  
  smooth_x <- unlist(df['V1'])
  smooth_y <- unlist(df['V2'])
  
  # Optimizing Constant-Span Running Mean (Don't Plot)
  opt_k <- plot_cvrss(smooth_y,TRUE,FALSE) 
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

# (sample, filled)
plot_data('ARID1A', FALSE)        
plot_data('ARID1A', TRUE)

plot_data('ARID1B', FALSE)        
plot_data('ARID1B', TRUE)

plot_data('ARID2', FALSE)        
plot_data('ARID2', TRUE)

# (sample, overlay, filled)
plot_csrm('ARID1A',FALSE,TRUE)    # FILLED      NOT OVERLAYED
plot_csrm('ARID1A',TRUE,TRUE)     #             OVERLAYED
plot_csrm('ARID1A',FALSE,FALSE)   # NOT FILLED  NOT OVERLAYED
plot_csrm('ARID1A',TRUE,FALSE)    #             OVERLAYED

plot_csrm('ARID1B',FALSE,TRUE)    # FILLED      NOT OVERLAYED
plot_csrm('ARID1B',TRUE,TRUE)     #             OVERLAYED
plot_csrm('ARID1B',FALSE,FALSE)   # NOT FILLED  NOT OVERLAYED
plot_csrm('ARID1B',TRUE,FALSE)    #             OVERLAYED

plot_csrm('ARID2',FALSE,TRUE)    # FILLED      NOT OVERLAYED
plot_csrm('ARID2',TRUE,TRUE)     #             OVERLAYED
plot_csrm('ARID2',FALSE,FALSE)   # NOT FILLED  NOT OVERLAYED
plot_csrm('ARID2',TRUE,FALSE)    #             OVERLAYED


# NON-PARAMETRIC BOOTSTRAPPING
# Takes empirical data and desired number of t statistics
# @data: sample, @n: Number of t statistics
bootstrap_mean <- function(data,n){
  b_means <- vector(mode="numeric", length=0);
  for(i in 1:n){
    b_mu <- mean(as.numeric(sample(data,length(data), replace=T)));   
    b_means <- c(b_means, b_mu);
  }
  return(b_means);
}

# calculates variance
# @bootstrap_estimates: bootstrap estimates of mean
calc_varianceEstimate <- function(bootstrap_estimates){
  trueEstimate <- mean(bootstrap_estimates);
  var <- 0;
  for(i in 1:length(bootstrap_estimates)){
    var <- var + (trueEstimate - bootstrap_estimates[i])**2;
  }
  return(var/(length(bootstrap_estimates)-1));
}

# Plots observed data, @sample: name of sample (E.g. ARID1A)
generate_bootstrapStatistics <- function(sample){
  # Reading in data
  df <- read.table(sprintf("../data/cleanedData/%s_filled.csv", sample), 
                   header = FALSE,
                   sep = ",")
  freq_at_positions <- unlist(df['V2'])
  
  # Generating bootstrap Means
  bs_means <- bootstrap_mean(freq_at_positions,1000) 
  bs_means_sorted <- sort(bs_means) 
  
  # Estimate mean and variance
  bs_meanEst <- mean(bs_means)                          
  bs_varEst <- calc_varianceEstimate(bs_means_sorted);  
  
  # Calculating Bootstrap-t C.I. of 95%
  bl <- mean(as.numeric(freq_at_positions)) - (bs_meanEst - bs_means_sorted[25]) * bs_varEst;   
  bu <- mean(as.numeric(freq_at_positions)) + (bs_means_sorted[975] - bs_meanEst) * bs_varEst; 
  
  print(sprintf("Bootstrap Mean: %f", bs_meanEst));
  print(sprintf("Bootstrap Variance: %f", bs_varEst));
  print(sprintf("Bootstrap .95 Confidence Interval: [%f,%f]", bl, bu));
}

generate_bootstrapStatistics('ARID1A')
generate_bootstrapStatistics('ARID1B')
generate_bootstrapStatistics('ARID2')

# Reading in data
# ECDF Creation
plot_ecdf <- function(sample){
  df <- read.table(sprintf("../data/cleanedData/%s_filled.csv", sample),
                   header = FALSE,
                   sep = ",")  
  freq <- unlist(df['V2'])
  l <- length(freq)
  
  saltus <- seq(from = 1/l, to = 1, by = 1/l)
  plot(sort(freq), saltus, main=sprintf("ECDF Plot - %s",sample),
       ylab="ECDF(x)",xlab="x")
}

plot_ecdf('ARID1A')
plot_ecdf('ARID1B')
plot_ecdf('ARID2')