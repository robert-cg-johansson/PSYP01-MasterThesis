#clear workspace
graphics.off()
rm(list=ls(all=TRUE))


#Set WD and chain parameters
setwd("C:/R.working.directory")
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

#Packages
library(BayesFactor)
library(rethinking)

#######################
### Import PFs Data ###
#######################

#Store the participants data in a preliminary dataframe
Data <- read.table("PF_nopauses.txt", quote="\"", comment.char="")

Data <- Data[17:416,]
names(Data) = c("trial", "Foreperiod", "RT", "brightness", "loudness", "Something", "RF")


# Invert RT in seconds to represent speed
Data$speed <-  1/Data$RT


# compute power
Data$Power <-  with(Data, RF*speed)


#Transform data from seconds to milliseconds for interpretability
Data[,3] <- Data[,3]*1000

#Remove outliers
Data <- Data[ which(Data$RT<1000  #Misses
), ]

Data <- Data[ which(Data$RT>100   #Anticipations
), ]



#Store the different trial types in separate dataframes for later analyses
D <- Data[ which(Data$brightness==90 
                 & Data$loudness==0.00), ]

B <- Data[ which(Data$brightness==240
                 & Data$loudness==0.00), ]

S <- Data[ which(Data$brightness==0
                 & Data$loudness==0.01), ]

L <- Data[ which(Data$brightness==0
                 & Data$loudness==1.0), ]

SD <- Data[ which(Data$brightness==90
                  & Data$loudness==0.01), ]

SB <- Data[ which(Data$brightness==240
                  & Data$loudness==0.01), ]

LD <- Data[ which(Data$brightness==90
                  & Data$loudness==1.0), ]

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]


#set trials per case for all trial types
tpcS  <- nrow(S)
tpcL  <- nrow(L)
tpcD  <- nrow(D)
tpcB  <- nrow(B)
tpcSD <- nrow(SD)
tpcLD <- nrow(LD)
tpcSB <- nrow(SB)
tpcLB <- nrow(LB)



#########################
### pRF Regression PF ###
#########################

###Create a dataset for RF regression

RF = c(S[,7], L[,7], D[,7], B[,7],SD[,7], LD[,7], SB[,7], LB[,7])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

RF.regression.PF <- cbind(RF, redundancy, loudness, brightness)
RF.regression.PF <- as.data.frame(RF.regression.PF)
names(RF.regression.PF) = c("RF","redundancy", "loudness", "brightness")



PF.redundancy <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=RF.regression.PF,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant PF, pRF ~ Redundancy, Stan estimates")
precis(PF.redundancy,prob=0.95,corr=TRUE)

lmBF(RF~redundancy, data=RF.regression.PF)

########################
### RT Regression PF ###
########################

###Create a dataset for RT regression

RT = c(S[,3], L[,3], D[,3], B[,3],SD[,3], LD[,3], SB[,3], LB[,3])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

RT.regression.PF <- cbind(RT, redundancy, loudness, brightness)
RT.regression.PF <- as.data.frame(RT.regression.PF)
names(RT.regression.PF) = c("RT","redundancy", "loudness", "brightness")



PF.redundancy <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=RT.regression.PF,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



print("Participant PF, RT ~ Redundancy, Stan estimates")
precis(PF.redundancy,prob=0.95,corr=TRUE)

lmBF(RT~redundancy, data=RT.regression.PF)


###########################
### Power Regression PF ###
###########################

###Create a dataset for power regression

Power = c(S[,9], L[,9], D[,9], B[,9],SD[,9], LD[,9], SB[,9], LB[,9])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

Power.regression.PF <- cbind(Power, redundancy, loudness, brightness)
Power.regression.PF <- as.data.frame(Power.regression.PF)
names(Power.regression.PF) = c("Power","redundancy", "loudness", "brightness")



PF.power <- ulam(
  alist(
    Power ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=Power.regression.PF,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



print("Participant PF, Power ~ Redundancy, Stan estimates")
precis(PF.power,prob=0.95,corr=TRUE)

lmBF(Power~redundancy, data=Power.regression.PF)


#######################
### Import GPs Data ###
#######################

#Store the participants data in a preliminary dataframe
Data <- read.table("GP_nopauses.txt", quote="\"", comment.char="")

Data <- Data[17:416,]
names(Data) = c("trial", "Foreperiod", "RT", "brightness", "loudness", "Something", "RF")

# Invert RT in seconds to represent speed
Data$speed <-  1/Data$RT

# compute power
Data$Power <-  with(Data, RF*speed)

#Transform data from seconds to milliseconds for interpretability
Data[,3] <- Data[,3]*1000

#Remove outliers
Data <- Data[ which(Data$RT<1000  #Misses
), ]

Data <- Data[ which(Data$RT>100   #Anticipations
), ]


#Store the different trial types in separate dataframes for later analyses
D <- Data[ which(Data$brightness==90 
                 & Data$loudness==0.00), ]

B <- Data[ which(Data$brightness==240
                 & Data$loudness==0.00), ]

S <- Data[ which(Data$brightness==0
                 & Data$loudness==0.01), ]

L <- Data[ which(Data$brightness==0
                 & Data$loudness==1.0), ]

SD <- Data[ which(Data$brightness==90
                  & Data$loudness==0.01), ]

SB <- Data[ which(Data$brightness==240
                  & Data$loudness==0.01), ]

LD <- Data[ which(Data$brightness==90
                  & Data$loudness==1.0), ]

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]


#set trials per case for all trial types
tpcS  <- nrow(S)
tpcL  <- nrow(L)
tpcD  <- nrow(D)
tpcB  <- nrow(B)
tpcSD <- nrow(SD)
tpcLD <- nrow(LD)
tpcSB <- nrow(SB)
tpcLB <- nrow(LB)

#######################################
### Multiple Regression for RF Data ###
#######################################

###Create a dataset for RF regression

RF = c(S[,7], L[,7], D[,7], B[,7],SD[,7], LD[,7], SB[,7], LB[,7])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

RF.regression.GP <- cbind(RF, redundancy, loudness, brightness)
RF.regression.GP <- as.data.frame(RF.regression.GP)
names(RF.regression.GP) = c("RF","redundancy", "loudness", "brightness")



GP.redundancy <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=RF.regression.GP,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant GP, redundancy effects, Stan estimates")
precis(GP.redundancy,prob=0.95,corr=TRUE)

lmBF(RF~redundancy, data=RF.regression.GP)


#######################################
### Multiple Regression for RT data ###
#######################################


###Create a dataset for RT regression

RT = c(S[,3], L[,3], D[,3], B[,3],SD[,3], LD[,3], SB[,3], LB[,3])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

RT.regression.GP <- cbind(RT, redundancy, loudness, brightness)
RT.regression.GP <- as.data.frame(RT.regression.GP)
names(RT.regression.GP) = c("RT","redundancy", "loudness", "brightness")



GP.redundancy <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=RT.regression.GP,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



print("Participant GP, RT ~ Redundancy, Stan estimates")
precis(GP.redundancy,prob=0.95,corr=TRUE)

lmBF(RT~redundancy, data=RT.regression.GP)

###########################
### Power Regression GP ###
###########################

###Create a dataset for power regression

Power = c(S[,9], L[,9], D[,9], B[,9],SD[,9], LD[,9], SB[,9], LB[,9])

redundancy = c(rep(0, tpcS+tpcL+tpcD+tpcB), rep(1, tpcSD+tpcLD+tpcSB+tpcLB))

loudness = c(rep(1, tpcS), rep(2, tpcL), rep(0, tpcD), rep(0, tpcB),
             rep(1, tpcSD), rep(2, tpcLD), rep(1, tpcSB), rep(2, tpcLB))

brightness= c(rep(0, tpcS), rep(0, tpcL), rep(1, tpcD), rep(2, tpcB),
              rep(1, tpcSD), rep(1, tpcLD), rep(2, tpcSB), rep(2, tpcLB))

Power.regression.GP <- cbind(Power, redundancy, loudness, brightness)
Power.regression.GP <- as.data.frame(Power.regression.GP)
names(Power.regression.GP) = c("Power","redundancy", "loudness", "brightness")



GP.power <- ulam(
  alist(
    Power ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * redundancy ,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=Power.regression.GP,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



print("Participant GP, Power ~ Redundancy, Stan estimates")
precis(GP.power,prob=0.95,corr=TRUE)

lmBF(Power~redundancy, data=Power.regression.GP)
