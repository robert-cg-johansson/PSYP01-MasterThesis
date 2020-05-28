####################
### SFT Analysis ###
####################


graphics.off()
rm(list=ls(all=TRUE))

#Set working directory
setwd("C:/R.working.directory/")

#Packages
library(sft) 

#Get the data
PF_Data <- read.table("PF_nopauses.txt", quote="\"", comment.char="")
GP_Data <- read.table("GP_nopauses.txt", quote="\"", comment.char="")

#Prepare Data
PF_Data <- PF_Data[17:416,]
names(PF_Data) = c("trial", "foreperiod", "RT", "brightness", "loudness", "button2", "RF")

GP_Data <- GP_Data[17:416,]
names(GP_Data) = c("trial", "foreperiod", "RT", "brightness", "loudness", "button2", "RF")

#Transform data from seconds to milliseconds for interpretability
PF_Data[,3] <- PF_Data[,3]*1000

GP_Data[,3] <- GP_Data[,3]*1000

#Remove outliers
PF_Data <- PF_Data[ which(PF_Data$RT<1000), ]  #Misses
PF_Data <- PF_Data[ which(PF_Data$RT>100), ]   #Anticipations

GP_Data <- GP_Data[ which(GP_Data$RT<1000), ]  #Misses
GP_Data <- GP_Data[ which(GP_Data$RT>100), ]   #Anticipations

#####################################
### Survivor Interaction Contrast ###
#####################################

windows(height=15,width=20)
par(mfrow=c(1,2))

#Store PFs bimodal RT data in frames for SIC

Data <- PF_Data

SD <- Data[ which(Data$brightness==90
                  & Data$loudness==0.01), ]

SB <- Data[ which(Data$brightness==240
                  & Data$loudness==0.01), ]

LD <- Data[ which(Data$brightness==90
                  & Data$loudness==1.0), ]

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]

#Enter PFs bimodal data for SIC analysis
hh <- LB[,3]
hl <- LD[,3]
lh <- SB[,3]
ll <- SD[,3]

#Run the SIC analysis
SerialAND <- sic(hh,hl,lh,ll)

#Plot the SIC-curve
plot(SerialAND$SIC, do.p=FALSE, ylab="SIC(t)", xlab="time (milliseconds)", main="Survivor Interaction Contrast (Participant PF)", xlim=c(200,450))



#Store GPs bimodal RT data in frames for SIC

Data <- GP_Data

SD <- Data[ which(Data$brightness==90
                  & Data$loudness==0.01), ]

SB <- Data[ which(Data$brightness==240
                  & Data$loudness==0.01), ]

LD <- Data[ which(Data$brightness==90
                  & Data$loudness==1.0), ]

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]

#Enter GPs bimodal data for SIC analysis
hh <- LB[,3]
hl <- LD[,3]
lh <- SB[,3]
ll <- SD[,3]

#Run the SIC analysis
SerialAND <- sic(hh,hl,lh,ll)

#Plot the SIC-curve
plot(SerialAND$SIC, do.p=FALSE, ylab="SIC(t)", xlab="time (milliseconds)", main="Survivor Interaction Contrast (Participant GP)", xlim=c(200,450))


#########################
### Capacity Analysis ###
  #########################
  
  windows(height=15,width=20)
  par(mfrow=c(2,2))
  
  ###PFs Data ###
  
  #Import Data
  Data <- PF_Data
  
  #Store PFs RT data in frames for capacity analysis
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
  
  LB <- Data[ which(Data$brightness==240
                    & Data$loudness==1.0), ]
  
  #Input single source response times (soft and dim)
  pa <- S[,3]
  ap <- D[,3]
  
  #Input multiple source response times (soft + dim)
  pp <- SD[,3]
  
  #Run the capacity analysis (soft+dim/soft and dim)
  cap_strong <- capacity.or(list(pp, pa, ap))
  cap_strong$Ctest
  plot(cap_strong$Ct, xlim=c(220,450), ylim=c(0,7), ylab = "C(t)", xlab="Time (msec)", main="Soft and Dim Signals")
  abline(h=1, lty=2)
  #Input single source response times
  pa <- L[,3]
  ap <- B[,3]
  
  #Input multiple source response times
  pp <- LB[,3]
  
  #Run the capacity analysis
  cap_strong <- capacity.or(list(pp, pa, ap))
  cap_strong$Ctest
  plot(cap_strong$Ct, xlim=c(220,450), ylim=c(0,7), ylab = "C(t)", xlab="Time (msec)", main="Loud and Bright Signals")
  abline(h=1, lty=2)
  mtext("Participant PF", outer = TRUE, line=-2, cex = 1.5)
  mtext("Participant GP", outer = TRUE, line=-26, cex = 1.5)
  
  
  ###GPs Data###
  
  #Import Data
  Data <- GP_Data
  
  #Store PFs RT data in frames for capacity analysis
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
  
  LB <- Data[ which(Data$brightness==240
                    & Data$loudness==1.0), ]
  
  #Input single source response times (soft and dim)
  pa <- S[,3]
  ap <- D[,3]
  
  #Input multiple source response times (soft + dim)
  pp <- SD[,3]
  
  #Run the capacity analysis (soft+dim/soft and dim)
  cap_strong <- capacity.or(list(pp, pa, ap))
  cap_strong$Ctest
  plot(cap_strong$Ct, xlim=c(220,450), ylim=c(0,9), ylab = "C(t)", xlab="Time (msec)", main="Soft and Dim Signals")
  abline(h=1, lty=2)
  
  #Input single source response times
  pa <- L[,3]
  ap <- B[,3]
  
  #Input multiple source response times
  pp <- LB[,3]
  
  #Run the capacity analysis
  cap_strong <- capacity.or(list(pp, pa, ap))
  cap_strong$Ctest
  plot(cap_strong$Ct, xlim=c(220,450), ylim=c(0,9), ylab = "C(t)", xlab="Time (msec)", main="Loud and Bright Signals")
  abline(h=1, lty=2)

#############################
### Race Model Inequality ###
#############################

#RMI algorithms for RT data as devised by Ulrich et al. (2007) 
#https://link.springer.com/article/10.3758/BF03193160?LI=true

#RMI algorithm kindly translated from MatLab to R by Y. Lin & L. Piwek via StackExchange. 
#https://stats.stackexchange.com/questions/48581/testing-the-race-model-inequality-in-r

#Probspace function
probSpace <- function(len){ 
  P <- numeric(len);
  for(i in 1:len){
    P[i] <- (i - .5) / len;
  }
  return(P)
}


###Specify Functions###

#CDF function
cdf.ulrich <- function(data=NULL, maximum=3000){
  # Create a container, whose length is the longest data vector
  # data is the output from `ties` function
  U <- data[,1]
  R <- data[,2]
  C <- data[,3]
  G <- numeric(maximum);
  
  
  # Length of processed data vector, trimming off ties, if there is any.
  k <- length(U);  # U contains data in millisecond, e.g., 320 ms etc. 
  
  # The last element of the cumulative frequency is thelength of the data vector.
  n <- C[k]
  
  for(i in 1:k) { U[i] <- round(U[i]) }
  
  # from 1 ms to the first value of the data set, their probability should be 0.
  for(t in 1:U[1]) { G[t] <- 0 }   
  
  for(t in U[1]:U[2]){
    G[t] <- ( R[1]/2 + (R[1]+R[2]) / 2*(t-U[1]) / (U[2] - U[1]) ) / n;
  }
  
  for(i in 2:(k-1)){
    for(t in U[i]:U[i+1]){
      G[t] <- (C[i-1] + R[i] / 2+(R[i] +R[i+1]) / 2*(t-U[i]) / (U[i+1] - U[i])) / n;
    }
  }
  
  for(t in U[k]:maximum){
    G[t] <- 1;
  }
  return(G)
}

#Ties function
ties <- function(W){
  # Count number k of unique values and store these values in U.
  U <- NULL; W <- sort(W); n = length(W); k = 1; U[1] <- W[1]
  for (i in 2:n) {
    if (W[i] != W[i-1]) {
      k <- k+1;
      U <- cbind(U, W[i])
    }
  }
  U <- U[1,]
  
  # Determine number of replications R
  # k is the length of the vector, after trimming off the ties
  R <- numeric(k) 
  for (i in 1:k){
    for (j in 1:n){
      if (U[i] == W[j]) R[i] <- R[i] + 1;
    }
  }
  
  # Determine the cumlative frequency
  C <- numeric(k) 
  C[1] <- R[1]
  for(i in 2:k){
    C[i] <- C[i-1] + R[i];
  }
  res <- list(U, R, C)
  names(res) <- c("U", "R", "C")
  return(as.data.frame(res))
}

##Get percentile function
GetPercentile <- function( P, G, tmax ){
  
  # Determine minimum of |G(Tp[i]) - P[i]|
  np <- length(P);
  Tp <- numeric(np)
  for( i in 1:np) {
    cc <- 100;
    for(t in 1:tmax) {
      if ( abs(G[t] - P[i]) < cc ) {
        c <- t;
        cc <- abs(G[t] - P[i]);
      }
    }
    
    if( P[i] > G[c] ){
      Tp[i] <-  c + (P[i] - G[c]) / (G[c+1] - G[c]);
    } else {
      Tp[i] <- c + (P[i] - G[c]) / (G[c] - G[c-1]);
    }
  }
  return( Tp )
}

###Compute and plot###

windows(height=15,width=20)

#Import PFs Data
Data <- PF_Data

#Store PFs RT data in frames for race model analysis
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

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]


#cx=channel 1, cy=channel 2, cz=redundant trials

#Soft and Dim signals
cx <- S[,3]
cy <- D[,3]
cz <- SD[,3]


psq <- probSpace(10); psq

dfx <- ties(cx)
dfy <- ties(cy)
dfz <- ties(cz)
tmax <- max(cx,cy,cz)

gx <- cdf.ulrich(data=dfx, maximum=tmax)
gy <- cdf.ulrich(data=dfy, maximum=tmax)
gz <- cdf.ulrich(data=dfz, maximum=tmax)

b <- gx + gy
xp <- GetPercentile(psq, gx, tmax)
yp <- GetPercentile(psq, gy, tmax)
zp <- GetPercentile(psq, gz, tmax);
bp <- GetPercentile(psq, b, tmax);

PF_RMI_1_data<- data.frame(RT =c(zp,bp), CDF =rep(psq, 2),
                  Condition =rep(c("Soft+Dim","Bounding sum"), each=length(xp)))

plot(y=PF_RMI_1_data$CDF, x=PF_RMI_1_data$RT)

PF_RMI_plot1 <- ggplot(PF_RMI_1_data,title="RMI", aes(x = RT, y = CDF, group=Condition,
                                       colour=Condition, shape=Condition)) + 
  geom_point() + geom_line() + theme_bw() +
  ggtitle("Soft and Dim Signals (PF)") +
  theme(plot.title = element_text(hjust = 0.5))

plot1 <- PF_RMI_plot1 + coord_cartesian(xlim = c(200, 410), ylim=c(-.01,1.01)) +
  theme(legend.position= c(.85, .20),  
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


#Loud and Bright signals
cx <- L[,3]
cy <- B[,3]
cz <- LB[,3]


psq <- probSpace(10); psq

dfx <- ties(cx)
dfy <- ties(cy)
dfz <- ties(cz)
tmax <- max(cx,cy,cz)

gx <- cdf.ulrich(data=dfx, maximum=tmax)
gy <- cdf.ulrich(data=dfy, maximum=tmax)
gz <- cdf.ulrich(data=dfz, maximum=tmax)

b <- gx + gy
xp <- GetPercentile(psq, gx, tmax)
yp <- GetPercentile(psq, gy, tmax)
zp <- GetPercentile(psq, gz, tmax);
bp <- GetPercentile(psq, b, tmax);

PF_RMI_2_data<- data.frame(RT =c(zp,bp), CDF =rep(psq, 2),
                           Condition =rep(c("Loud+Bright","Bounding sum"), each=length(xp)))

PF_RMI_plot2 <- ggplot(PF_RMI_2_data,title="RMI", aes(x = RT, y = CDF, group=Condition,
                                                      colour=Condition, shape=Condition)) + 
  geom_point() + geom_line() + theme_bw() +
  ggtitle("Loud and Bright Signals (PF)") +
  theme(plot.title = element_text(hjust = 0.5))

plot2 <- PF_RMI_plot2 + coord_cartesian(xlim = c(200, 410), ylim=c(-.01,1.01)) +
  theme(legend.position= c(.85, .20),  
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


#Import GPs Data
Data <- GP_Data

#Store GPs RT data in frames for race model analysis
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

LB <- Data[ which(Data$brightness==240
                  & Data$loudness==1.0), ]


#cx=channel 1, cy=channel 2, cz=redundant trials

#Soft and Dim signals
cx <- S[,3]
cy <- D[,3]
cz <- SD[,3]


psq <- probSpace(10); psq

dfx <- ties(cx)
dfy <- ties(cy)
dfz <- ties(cz)
tmax <- max(cx,cy,cz)

gx <- cdf.ulrich(data=dfx, maximum=tmax)
gy <- cdf.ulrich(data=dfy, maximum=tmax)
gz <- cdf.ulrich(data=dfz, maximum=tmax)

b <- gx + gy
xp <- GetPercentile(psq, gx, tmax)
yp <- GetPercentile(psq, gy, tmax)
zp <- GetPercentile(psq, gz, tmax);
bp <- GetPercentile(psq, b, tmax);

GP_RMI_1_data<- data.frame(RT =c(zp,bp), CDF =rep(psq, 2),
                           Condition =rep(c("Soft+Dim","Bounding sum"), each=length(xp)))

GP_RMI_plot1 <- ggplot(GP_RMI_1_data,title="RMI", aes(x = RT, y = CDF, group=Condition,
                                                      colour=Condition, shape=Condition)) + 
  geom_point() + geom_line() + theme_bw() +
  ggtitle("Soft and Dim Signals (GP)") +
  theme(plot.title = element_text(hjust = 0.5))

plot3 <- GP_RMI_plot1 + coord_cartesian(xlim = c(200, 410), ylim=c(-.01,1.01)) +
  theme(legend.position= c(.85, .20),  
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


#Loud and Bright signals
cx <- L[,3]
cy <- B[,3]
cz <- LB[,3]


psq <- probSpace(10); psq

dfx <- ties(cx)
dfy <- ties(cy)
dfz <- ties(cz)
tmax <- max(cx,cy,cz)

gx <- cdf.ulrich(data=dfx, maximum=tmax)
gy <- cdf.ulrich(data=dfy, maximum=tmax)
gz <- cdf.ulrich(data=dfz, maximum=tmax)

b <- gx + gy
xp <- GetPercentile(psq, gx, tmax)
yp <- GetPercentile(psq, gy, tmax)
zp <- GetPercentile(psq, gz, tmax);
bp <- GetPercentile(psq, b, tmax);

GP_RMI_2_data<- data.frame(RT =c(zp,bp), CDF =rep(psq, 2),
                           Condition =rep(c("Loud+Bright","Bounding sum"), each=length(xp)))


plot(GP_RMI_2_data)
GP_RMI_plot2 <- ggplot(GP_RMI_2_data,title="RMI", aes(x = RT, y = CDF, group=Condition,
                                                      colour=Condition, shape=Condition)) + 
  geom_point() + geom_line() + theme_bw() +
  ggtitle("Loud and Bright Signals (GP)") +
  theme(plot.title = element_text(hjust = 0.5))

plot4 <- GP_RMI_plot2 + coord_cartesian(xlim = c(200, 410), ylim=c(-.01,1.01)) +
  theme(legend.position= c(.85, .20),  
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

###open all 4 plots in a window###
grid.arrange(plot1, plot2,
             plot3, plot4, ncol=2, nrow=2)







#Loud and Bright signals
cx <- L[,7]
cy <- B[,7]
cz <- LB[,7]


psq <- probSpace(10); psq

dfx <- ties(cx)
dfy <- ties(cy)
dfz <- ties(cz)
tmax <- max(cx,cy,cz)

gx <- cdf.ulrich(data=dfx, maximum=tmax)
gy <- cdf.ulrich(data=dfy, maximum=tmax)
gz <- cdf.ulrich(data=dfz, maximum=tmax)

b <- gx + gy
xp <- GetPercentile(psq, gx, tmax)
yp <- GetPercentile(psq, gy, tmax)
zp <- GetPercentile(psq, gz, tmax);
bp <- GetPercentile(psq, b, tmax);

GP_RMI_2_data<- data.frame(RT =c(zp,bp), CDF =rep(psq, 2),
                           Condition =rep(c("Loud+Bright","Bounding sum"), each=length(xp)))


plot(GP_RMI_2_data)
GP_RMI_plot2 <- ggplot(GP_RMI_2_data,title="RMI", aes(x = RT, y = CDF, group=Condition,
                                                      colour=Condition, shape=Condition)) + 
  geom_point() + geom_line() + theme_bw() +
  ggtitle("Loud and Bright Signals (GP)") +
  theme(plot.title = element_text(hjust = 0.5))

GP_RMI_plot2 + coord_cartesian(xlim = c(0, 20), ylim=c(-.01,1.01)) +
  theme(legend.position= c(.85, .20),  
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))


