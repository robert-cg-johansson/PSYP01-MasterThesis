#clear workspace
graphics.off()
rm(list=ls(all=TRUE))


#Set working directory
setwd("C:/R.working.directory")

#Packages
library(rethinking)

# detect cores and leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

#Get the data
PF_Data <- read.table("PF_nopauses.txt", quote="\"", comment.char="")
GP_Data <- read.table("GP_nopauses.txt", quote="\"", comment.char="")

PF_Data <- PF_Data[-c(1:16),]
names(PF_Data) = c("trial", "Foreperiod", "RT", "brightness", "loudness", "Something", "RF")

GP_Data <- GP_Data[-c(1:16),]
names(GP_Data) = c("trial", "Foreperiod", "RT", "brightness", "loudness", "Something", "RF")


#Transform RT data from seconds to milliseconds for interpretability
PF_Data[,"RT"] <- PF_Data[,"RT"]*1000
GP_Data[,"RT"] <- GP_Data[,"RT"]*1000

#Remove outliers
PF_Data <- PF_Data[ which(PF_Data$RT<1000), ]  #Misses
GP_Data <- GP_Data[ which(GP_Data$RT<1000), ]  #Misses

PF_Data <- PF_Data[ which(PF_Data$RT>100 ), ]  #Anticipations
GP_Data <- GP_Data[ which(GP_Data$RT>100 ), ]  #Anticipations


#Standardize the data
PF_Data$ZRT<- with(PF_Data, (RT-mean(RT))/sd(RT))
PF_Data$ZRF<- with(PF_Data, (RF-mean(RF))/sd(RF))
GP_Data$ZRT<- with(GP_Data, (RT-mean(RT))/sd(RT))
GP_Data$ZRF<- with(GP_Data, (RF-mean(RF))/sd(RF))

#############################################
# Relations (correlation) between RT and RF #
#############################################

PF_m1 <- map2stan( # from McElreath 2nd edition
  alist(
    ZRF ~ normal( mu, sigma ),
    mu <- b0 + b1*ZRT,
    c(b0,b1) ~ dnorm(0, 1),
    sigma ~ dcauchy(0,1)
  ) ,
  data=PF_Data, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(PF_m1,depth=3,prob=0.95,corr=TRUE)

GP_m1 <- ulam(
  alist(
    ZRF ~ normal( mu, sigma ),
    mu <- b0 + b1*ZRT,
    c(b0,b1) ~ dnorm(0, 1),
    sigma ~ dcauchy(0,1)
  ) ,
  data=GP_Data, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(GP_m1,depth=3,prob=0.95,corr=TRUE)


# 2 panel plot. 1 Scatterplot for PF and 1 Scatterplot for GP
windows(height=10,width=15) #use X11 for Mac OS
par( mar=0.5+c(5,4,1,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

with(PF_Data, plot(ZRF~ZRT, main="PF", xlab=expression(paste("Z(RT [ms])")),ylab=expression(paste("Z(RF [Newtons])")), xlim=c(-2, 8), ylim=c(-2, 5),axes=T,bty="n"))
abline(PF_m1)
reg<-lm(ZRF ~ ZRT, data = PF_Data)
abline(reg)

with(GP_Data, plot(ZRF~ZRT, main="GP", xlab=expression(paste("Z(RT [ms])")),ylab=expression(paste("Z(RF [Newtons])")), xlim=c(-2, 8), ylim=c(-2, 5),axes=T,bty="n"))
reg2<-lm(ZRF ~ ZRT, data = GP_Data)
abline(reg2)

abline(GP_m1)

################################################
# Relations (correlation) between RT and trial #
################################################


PF_m2 <- map2stan( # from McElreath 2nd edition
  alist(
    RF ~ normal( mu, sigma ),
    mu <- b0 + b1*trial,
    c(b0,b1) ~ dnorm(0, 1),
    sigma ~ dcauchy(0,1)
  ) ,
  data=PF_Data, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(PF_m2,depth=3,prob=0.95,corr=TRUE)


GP_m2 <- ulam(
  alist(
    RF ~ normal( mu, sigma ),
    mu <- b0 + b1*trial,
    c(b0,b1) ~ dnorm(0, 1),
    sigma ~ dcauchy(0,1)
  ) ,
  data=GP_Data, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(GP_m1,depth=3,prob=0.95,corr=TRUE)


# 4 panel plot. 2 Scatterplots for PF and 2 Scatterplots for GP
windows(height=10,width=15) #use X11 for Mac OS
par( mar=0.5+c(5,4,1,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

with(PF_Data, plot(ZRF~ZRT, main="PF", xlab=expression(paste("Z(RT [ms])")),ylab=expression(paste("Z(RF [Newtons])")), xlim=c(-2, 8), ylim=c(-2, 5),axes=T,bty="n"))
reg1 <- lm(ZRF~ZRF, data=PF_Data)
abline(reg1)

with(GP_Data, plot(ZRF~ZRT, main="GP", xlab=expression(paste("Z(RT [ms])")),ylab=expression(paste("Z(RF [Newtons])")), xlim=c(-2, 8), ylim=c(-2, 5),axes=T,bty="n"))
reg2 <- lm(ZRF~ZRF, data=GP_Data)
abline(GP_m1)


with(PF_Data, plot(RF~trial, main="PF", xlab=expression(paste("Trial Number")),ylab=expression(paste("RF (Newtons)")), xlim=c(0, 410), ylim=c(0, 15),axes=T,bty="n"))
reg2 <- lm(RF~trial, data=PF_Data)
abline(reg2)

with(GP_Data, plot(RF~trial, main="GP", xlab=expression(paste("Trial Number")),ylab=expression(paste("RF (Newtons)")), xlim=c(0, 410), ylim=c(0, 15),axes=T,bty="n"))
reg3 <- lm(RF~trial, data=GP_Data)
abline(reg3)
