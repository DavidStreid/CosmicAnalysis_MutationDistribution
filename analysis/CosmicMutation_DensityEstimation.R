# Density Estimation

# Requires the install of package "RCurl" by typing in cmd line: install.packages("RCurl")
require(RCurl)

# Import data
arid1a <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID1A.csv"))
arid1b <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID1B.csv"))
arid2 <- read.csv(text=getURL("https://raw.githubusercontent.com/DavidStreid/CosmicAnalysis_MutationDistribution/master/data/cleanedData/ARID2.csv"))

# Function for f.hat
f.hat = function(x,h,data){
  mean(dnorm((x-data)/h))/h
}

# Silverman's Rule bandwidth selection as pilot for Sheather-Jones
h <- function(data,n){
  data.sd <- sd(data)
  out <- ((4/(3*n))^(1/5))*data.sd
  return(out)
}

# Sheather-Jones Bandwidth Estimate using the Normal as the roughness function

#############################
##### ARID1A
#############################

arid1a.x.data <- arid1a$X1
n0 <- length(arid1a.x.data)
arid1a.x = seq(min(arid1a.x.data)-sd(arid1a.x.data), max(arid1a.x.data)+sd(arid1a.x.data), length.out=n0)

sd0 <- sd(arid1a.x.data)

L <- f.hat(arid1a.x,h(arid1a.x.data,n0),arid1a.x.data)
h0 <- ((1/(2*sqrt(pi)))/(n0*(sd0^4)*L))^(1/5)
print(h0)

d0 <- matrix(0,length(arid1a.x),length(h0))
d0 <- mapply(f.hat,arid1a.x,h0,arid1a.x.data)

hist(arid1a.x.data,breaks=20,freq=FALSE, xlim=c(min(arid1a.x),max(arid1a.x)),main="ARID1A Density Estimate")
lines(arid1a.x,d0,col="blue")

#############################
##### ARID1B
#############################

arid1b.x.data <- arid1b$X1
n1 <- length(arid1b.x.data)
arid1b.x = seq(min(arid1b.x.data)-sd(arid1b.x.data), max(arid1b.x.data)+sd(arid1b.x.data), length.out=n1)

sd1 <- sd(arid1b.x.data)

L1 <- f.hat(arid1b.x,h(arid1b.x.data,n1),arid1b.x.data)
h1 <- ((1/(2*sqrt(pi)))/(n1*(sd1^4)*L))^(1/5)
print(h1)

d1 <- matrix(0,length(arid1b.x),length(h1))
d1 <- mapply(f.hat,arid1b.x,h1,arid1b.x.data)

hist(arid1b.x.data,breaks=20,freq=FALSE, xlim=c(min(arid1b.x),max(arid1b.x)),main="ARID1B Density Estimate")
lines(arid1b.x,d1,col="blue")

#############################
##### ARID2
#############################

arid2.x.data <- arid2$X1
n2 <- length(arid2.x.data)
arid2.x = seq(min(arid2.x.data)-sd(arid2.x.data), max(arid2.x.data)+sd(arid2.x.data), length.out=n2)

sd2 <- sd(arid2.x.data)

L2 <- f.hat(arid2.x,h(arid2.x.data,n2),arid2.x.data)
h2 <- ((1/(2*sqrt(pi)))/(n2*(sd2^4)*L))^(1/5)
print(h2)

d2 <- matrix(0,length(arid2.x),length(h2))
d2 <- mapply(f.hat,arid2.x,h2,arid2.x.data)

hist(arid2.x.data,breaks=20,freq=FALSE, xlim=c(min(arid2.x),max(arid2.x)),main="ARID2 Density Estimate")
lines(arid2.x,d2,col="blue")