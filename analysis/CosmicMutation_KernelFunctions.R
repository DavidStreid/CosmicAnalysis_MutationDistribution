## Kernel Functions

# Requires the install of package "RCurl" by typing in cmd line: install.packages("RCurl")
require(RCurl)

# Import data
arid1a <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID1A.csv"))
arid1b <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID1B.csv"))
arid2 <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID2.csv"))

# Arranging data
arid1a.x.data <- arid1a$X26696416
arid1a.y.data <- arid1a$X1
sort.data <- sort(arid1a.x.data, index.return=TRUE)
arid1a.x <- arid1a.x.data[sort.data$ix]
arid1a.y <- arid1a.y.data[sort.data$ix]

arid1b.x.data <- arid1b$X156778152
arid1b.y.data <- arid1b$X1
sort.data <- sort(arid1b.x.data, index.return=TRUE)
arid1b.x <- arid1b.x.data[sort.data$ix]
arid1b.y <- arid1b.y.data[sort.data$ix]

arid2.x.data <- arid2$X45729909
arid2.y.data <- arid2$X1
sort.data <- sort(arid2.x.data, index.return=TRUE)
arid2.x <- arid2.x.data[sort.data$ix]
arid2.y <- arid2.y.data[sort.data$ix]

# Constant Span Running Mean Smoother

################################
# CVRSS for Span Selection
################################

k <- seq(3,31,by=2)

cvrss <- function(k,y){
  n <- length(y)
  S <- matrix(0,n,n)
  a <- (k-1)/2
  
  if(k > 1){
    for(i in 1:a){
      S[n-i+1,(n-a-i+1):n] <- 1/(k-a+i-1)
      S[i,1:(a+i)] <- 1/(k-a+i-1)
    }
    for(i in (a+1):(n-a)){
      S[i,(i-a):(i+a)] <- 1/k
    }
  }
  
  if(k==1){
    S <- diag(1,n)
  }
  
  s.hat <- S%*%y
  out <- sum(((y-s.hat)/(1-diag(S)))^2)
  return(out)
}

onea.cvrss.nums <- rep(0,length(k))
oneb.cvrss.nums <- rep(0,length(k))
two.cvrss.nums <- rep(0,length(k))

for(i in 1:14){
  onea.cvrss.nums[i] <- cvrss(k[i],arid1a.y)
  oneb.cvrss.nums[i] <- cvrss(k[i],arid1b.y)
  two.cvrss.nums[i] <- cvrss(k[i],arid2.y)
}

plot(k,onea.cvrss.nums,type='b',main="ARID1A CVRSS Plot")
plot(k,oneb.cvrss.nums,type='b',main="ARID1B CVRSS Plot")
plot(k,two.cvrss.nums,type='b',main="ARID2 CVRSS Plot")

# From analysis of the CVRSS plot, it would be best to use a k value of 13+.
# Arbitraily choosing, value of 15 was chosen for the kernel smoother.

######################################
# Using truncated neighborhood average
######################################

trunc <- function(k,y){
  n <- length(y)
  S <- matrix(0,n,n)
  a <- (k-1)/2
  
  if(k > 1){
    for(i in 1:a){
      S[i,1:(a+i)] <- 1/(k-a+i-1)
      S[n-i+1,(n-a-i+1):n] <- 1/(k-a+i-1)
    }
    
    for(i in (a+1):(n-a)){
      S[i,(i-a):(i+a)] <- 1/k
    }
  }
  
  if(k==1){
    S <- diag(1,n)
  }
  
  out = S%*%y
  return(out)
}

onea.s15 <- trunc(15,arid1a.y)
oneb.s15 <- trunc(15,arid1b.y)
two.s29 <- trunc(29,arid2.y)

s <- function(x){
  (x^3) * sin((x+3.4)/2)
}

################### ARID1A CSRM Plot ###############################
arid1a.x.plot <- seq(min(arid1a.x),max(arid1a.x),length.out=1000)
arid1a.y.plot <- s(arid1a.x.plot)
plot(arid1a.x,arid1a.y,xlab="Predictor",ylab="Response",main="Constant-Span Running Mean")
lines(arid1a.x,onea.s15,col="blue")

################### ARID1B CSRM Plot ###############################
arid1b.x.plot <- seq(min(arid1b.x),max(arid1b.x),length.out=1000)
arid1b.y.plot <- s(arid1b.x.plot)
plot(arid1b.x,arid1b.y,xlab="Predictor",ylab="Response",main="ARID1B Constant-Span Running Mean")
lines(arid1b.x,oneb.s15,col="black")

################### ARID2 CSRM Plot ###############################
arid2.x.plot <- seq(min(arid2.x),max(arid2.x),length.out=1000)
arid2.y.plot <- s(arid2.x.plot)
plot(arid2.x,arid2.y,xlab="Predictor",ylab="Response",main="ARID1B Constant-Span Running Mean")
lines(arid2.x,two.s29,col="black")


############################
# Kernel smoother
############################


h <- 100
h1 <- 250
h2 <- 300

kz0 = function(z,h){
  dnorm((z-arid1a.x)/h)/h
}

kz1 = function(z,h1){
  dnorm((z-arid1b.x)/h1)/h1
}

kz2 = function(z,h2){
  dnorm((z-arid2.x)/h2)/h2
}

kernel = function(h,y,x){
  n = length(y)
  s.hat = rep(0,n)
  
  for(i in 1:n){
    a = kz0(x[i],h)
    s.hat[i] = sum(y * a/sum(a))
  }
  
  return(s.hat)
}

kernel1 = function(h,y,x){
  n = length(y)
  s.hat = rep(0,n)
  
  for(i in 1:n){
    a1 = kz1(x[i],h1)
    s.hat[i] = sum(y * a1/sum(a1))
  }
  
  return(s.hat)
}

kernel2 = function(h,y,x){
  n = length(y)
  s.hat = rep(0,n)
  
  for(i in 1:n){
    a2 = kz2(x[i],h2)
    s.hat[i] = sum(y * a2/sum(a2))
  }
  
  return(s.hat)
}

onea.kernel.s = kernel(h,arid1a.y,arid1a.x)
oneb.kernel.s = kernel1(h,arid1b.y,arid1b.x)
two.kernel.s = kernel2(h,arid2.y,arid2.x)
  
## OUTPUT PLOTS
s = function(x){
  (x^3) * sin((x+3.4)/2)
}
  
################### ARID1A Kernel Plot ########################
onea.x.plot = seq(min(arid1a.x),max(arid1a.x),length.out=1000)
onea.y.plot = s(onea.x.plot)
plot(arid1a.x,arid1a.y,xlab="Predictor",ylab="Response",main="ARID1A Kernel")
lines(arid1a.x,onea.kernel.s,type="l",col="blue")

################### ARID1B CSRM Plot ###############################
oneb.x.plot = seq(min(arid1b.x),max(arid1b.x),length.out=1000)
oneb.y.plot = s(oneb.x.plot)
plot(arid1b.x,arid1b.y,xlab="Predictor",ylab="Response",main="ARID1B Kernel")
lines(arid1b.x,oneb.kernel.s,type="l",col="blue")

################### ARID2 CSRM Plot ###############################
two.x.plot = seq(min(arid2.x),max(arid2.x),length.out=1000)
two.y.plot = s(two.x.plot)
plot(arid2.x,arid2.y,xlab="Predictor",ylab="Response",main="ARID2 Kernel")
lines(arid2.x,two.kernel.s,type="l",col="blue")