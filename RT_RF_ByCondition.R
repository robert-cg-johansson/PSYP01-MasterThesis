#clear workspace
graphics.off()
rm(list=ls(all=TRUE))


#####################################################################
#
#  TAKES ABOUT 20 MINS TO RUN
#  You must have Utilities.R in your working directory
#
######################################################################


#Set working directory
setwd("C:/R.working.directory/")

# detect cores and leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

#Get the data
PF_Data <- read.table("PF_nopauses.txt", quote="\"", comment.char="")
GP_Data <- read.table("GP_nopauses.txt", quote="\"", comment.char="")

source("Utilities.R")
PF_Data=prep(PF_Data)
GP_Data=prep(GP_Data)


######################################################
# Participant PF analysis
##################################################

PF_d_uni=prep_unimodal(PF_Data)
PF_d_bi=prep_bimodal(PF_Data)

##############################################
# PF Reaction time

require(rethinking)

# unclear how ulam works with factors so use numeric coding instead

PF_RT_Unim1 <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(400, 1e4),
    b1 ~ dnorm(-1, 1e4),
    b2 ~ dnorm(-40, 1e4),
    b3 ~ dnorm(-20, 1e4),
    sigma ~ dcauchy(100,1e4)
  ) ,
  data=PF_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


PF_RT_Bim2 <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(300, 1e4),
    b1 ~ dnorm(-20, 1e4),
    b2 ~ dnorm(-20, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(100,1e4)
  ) ,
  data=PF_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



print("Participant PF, unimdal RT, Stan estimates")
precis(PF_RT_Unim1,prob=0.95,corr=TRUE)



print("Participant PF, Bimodal RT, Stan estimates")
precis(PF_RT_Bim2,prob=0.95,corr=TRUE)

print("Participant PF, unimdal / Bimodal RT, maximum likelihood and BFs")
BF_anal(PF_d_uni, PF_d_bi, DV="RT")


##############################################
# PF Response force

PF_RF_Unim1a <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=PF_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



PF_RF_Bim2a <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=PF_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)
print("Participant PF, unimdal RF, Stan estimates")
precis(PF_RF_Unim1a,prob=0.95,corr=TRUE)

print("Participant PF, Bimodal RF, Stan estimates")
precis(PF_RF_Bim2a,prob=0.95,corr=TRUE)

print("Participant PF, unimdal / Bimodal RF, maximum likelihood and BFs")
BF_anal(PF_d_uni, PF_d_bi, DV="RF")


##############################################
# PF Power

PF_Pow_Unim1b <- ulam(
  alist(
    Pow ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(10, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=PF_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



PF_Pow__Bim2b <- ulam(
  alist(
    Pow ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(10, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=PF_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant PF, unimdal Power, Stan estimates")
precis(PF_Pow_Unim1b, prob=0.95, corr=TRUE)

print("Participant PF, Bimodal Power, Stan estimates")
precis(PF_Pow__Bim2b, prob=0.95, corr=TRUE)

print("Participant PF, unimdal / Bimodal Power, maximum likelihood and BFs")
BF_anal(PF_d_uni, PF_d_bi, DV="Pow")


#################################################################################
# Repeat for GP
#################################################################################

GP_d_uni=prep_unimodal(GP_Data)
GP_d_bi=prep_bimodal(GP_Data)

##############################################
# PF  Reaction time

require(rethinking)

# unclear how ulam works with factors so use numeric coding instead

GP_RT_Unim1 <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(400, 1e4),
    b1 ~ dnorm(-1, 1e4),
    b2 ~ dnorm(-40, 1e4),
    b3 ~ dnorm(-20, 1e4),
    sigma ~ dcauchy(100,1e4)
  ) ,
  data=GP_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



GP_RT_Bim2 <- ulam(
  alist(
    RT ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(300, 1e4),
    b1 ~ dnorm(-20, 1e4),
    b2 ~ dnorm(-20, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(100,1e4)
  ) ,
  data=GP_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant GP, unimdal RT, Stan estimates")
precis(GP_RT_Unim1, prob=0.95, corr=TRUE)

print("Participant GP, Bimdal RT, Stan estimates")
precis(GP_RT_Bim2, prob=0.95,corr=TRUE)

print("Participant GP, unimdal / Bimodal RT, maximum likelihood and BFs")
BF_anal(GP_d_uni, GP_d_bi, DV="RT")
library(BayesFactor)

lm(formula = RT ~ modality + intensity + modality:intensity, data = PF_d_uni)
##############################################
# GP Response force

GP_RF_Unim1a <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=GP_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



GP_RF_Bim2a <- ulam(
  alist(
    RF ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(4, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=GP_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant GP, unimdal RF, Stan estimates")
precis(GP_RF_Unim1a,prob=0.95,corr=TRUE)

print("Participant GP, bimdal RF, Stan estimates")
precis(GP_RF_Bim2a,prob=0.95,corr=TRUE)


print("Participant GP, unimdal / Bimodal RF, maximum likelihood and BFs")
BF_anal(GP_d_uni, GP_d_bi, DV="RF")


##############################################
# GP Power

GP_Pow_Unim1b <- ulam(
  alist(
    Pow ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * mod_N + b2 * int_N + b3 * mod_N * int_N,
    b0 ~ dnorm(10, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=GP_d_uni,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)



GP_Pow_Bim2b <- ulam(
  alist(
    Pow ~ dnorm( mu, sigma ),
    mu <- b0 + b1 * aud_N + b2 * vis_N + b3 * aud_N * vis_N,
    b0 ~ dnorm(10, 1e4),
    b1 ~ dnorm(0, 1e4),
    b2 ~ dnorm(0, 1e4),
    b3 ~ dnorm(0, 1e4),
    sigma ~ dcauchy(0,1e4)
  ) ,
  data=GP_d_bi,
  constraints = list(sigma = "lower=0"),
  chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)


print("Participant GP, unimdal RF, Stan estimates")
precis(GP_Pow_Unim1b,prob=0.95,corr=TRUE)

print("Participant GP, bimdal RF, Stan estimates")
precis(GP_Pow_Bim2b,prob=0.95,corr=TRUE)


print("Participant GP, unimdal / Bimodal Pow, maximum likelihood and BFs")
BF_anal(GP_d_uni, GP_d_bi, DV="Pow")





###############################################################
# Plot results
#######################################################

windows(height=15,width=20)
par( mar=0.5+c(5,4,6,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(2,4))



# 1st panel unimodal RT, participant PF
modality=0 # visual
muVis <- link(PF_RT_Unim1, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(PF_RT_Unim1, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(100,1000), main="Reaction Time (RT)", xlab="Signal intensity", ylab="RT (msec)", type = "n", xaxt="n", bty="l")
with(PF_d_uni, points(RT ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
#shade(x,y) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)

legend(0.1, 950, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)
mtext(text="Unimodal",side=3,line=3.5)


##############################################################################
# 2nd panel biimodal RT, participant PF

auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(PF_RT_Bim2, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(PF_RT_Bim2, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(100,1000), main="Reaction Time (RT)", xlab="Auditory component intensity", ylab="RT (msec)", type = "n", xaxt="n", bty="l")
with(PF_d_bi, points(RT ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)


legend(0, 950, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)
mtext(text="Bimodal",side=3, line=3.5)


##############################################################################
# 3rd panel unimodal RF, participant PF

modality=0 # visual
muVis <- link(PF_RF_Unim1a, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(PF_RF_Unim1a, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,20), main="Response Force (RF)", xlab="Signal intensity", ylab="RF (Newtons)", type = "n", xaxt="n", bty="l")
with(PF_d_uni, points(RF ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
#shade(muVis.HPDI,0:1) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)

#legend(0.1, 19, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)
mtext(text="Unimodal",side=3,line=3.5)


#########################################################################
# 4th panel biimodal RF, participant PF

auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(PF_RF_Bim2a, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(PF_RF_Bim2a, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,20), main="Response Force (RF)", xlab="Auditory component intensity", ylab="RF (Newtons)", type = "n", xaxt="n", bty="l")
with(PF_d_bi, points(RF ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)

#legend(0, 19, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)
mtext(text="Bimodal",side=3, line=3.5)

mtext("Participant PF", outer = TRUE, cex = 1.5)



#############################################################
# GP
# 5th panel unimodal RT, participant GP
modality=0 # visual
muVis <- link(GP_RT_Unim1, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(GP_RT_Unim1, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(100,1000), main="Reaction Time (RT)", xlab="Signal intensity", ylab="RT (msec)", type = "n", xaxt="n", bty="l")
with(GP_d_uni, points(RT ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
#shade(x,y) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)


#legend(0.1, 950, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)
mtext(text="Unimodal",side=3,line=3.5)


#######################################################
# 6th panel biimodal RT, participant GP
auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(GP_RT_Bim2, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(GP_RT_Bim2, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(100,1000), main="Reaction Time (RT)", xlab="Auditory component intensity", ylab="RT (msec)", type = "n", xaxt="n", bty="l")
with(GP_d_bi, points(RT ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)

#legend(0, 950, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)
mtext(text="Bimodal",side=3, line=3.5)


####################################################
# 7th panel unimodal RF, participant GP
modality=0 # visual
muVis <- link(GP_RF_Unim1a, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(GP_RF_Unim1a, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,20), main="Response Force (RF)", xlab="Signal intensity", ylab="RF (Newtons)", type = "n", xaxt="n", bty="l")
with(GP_d_uni, points(RF ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
#shade(x,y) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)

#legend(0.1, 19, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)
mtext(text="Unimodal",side=3,line=3.5)


##########################################################
# 8th panel biimodal RF, participant GP

auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(GP_RF_Bim2a, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(GP_RF_Bim2a, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,20), main="Response Force (RF)", xlab="Auditory component intensity", ylab="RF (Newtons)", type = "n", xaxt="n", bty="l")
with(GP_d_bi, points(RF ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)

#legend(0, 19, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)
mtext(text="Bimodal",side=3, line=3.5)

mtext("Participant GP", outer = TRUE, line=-28, cex = 1.5)


########################################################################
#
# Seperate 4 panel plot for Power
#
########################################################################

windows(height=15,width=15)
par( mar=0.5+c(5,4,6,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(2,2))


# 1st panel unimodal Power, participant PF
modality=0 # visual
muVis <- link(PF_Pow_Unim1b, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(PF_Pow_Unim1b, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,70), main="Unimodal", xlab="Signal intensity", ylab="Power", type = "n", xaxt="n", bty="l")
with(PF_d_uni, points(Pow ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
# shade(x,y) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)

legend(0.1, 65, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)


# 2nd panel biimodal Power, participant PF

auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(PF_Pow__Bim2b, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(PF_Pow__Bim2b, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,70), main="Bimodal", xlab="Auditory component intensity", ylab="Power", type = "n", xaxt="n", bty="l")
with(PF_d_bi, points(Pow ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)

legend(0.1, 65, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)

mtext("Participant PF", outer = TRUE, cex = 1.5)

#############################
# Participant GP
# 3rd panel unimodal Power, participant GP
modality=0 # visual
muVis <- link(GP_Pow_Unim1b, data=data.frame(mod_N=0, int_N=c(0,1)))
muVis.mean <- apply(muVis,2,mean)
muVis.HPDI <- apply(muVis,2,HPDI, prob=0.95)

modality=1 # auditory
muAud <- link(GP_Pow_Unim1b, data=data.frame(mod_N=1, int_N=c(0,1)))
muAud.mean <- apply(muAud,2,mean)
muAud.HPDI <- apply(muAud,2,HPDI,prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,70), main="Unimodal", xlab="Signal intensity", ylab="Power", type = "n", xaxt="n", bty="l")
with(GP_d_uni, points(Pow ~ int_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVis.HPDI[2,1], muVis.HPDI[2,2],muVis.HPDI[1,2],muVis.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)
# shade(x,y) # throws a warning

x<-c(0.02,0.99,0.99,0.02)
y<-c(muAud.HPDI[2,1], muAud.HPDI[2,2],muAud.HPDI[1,2],muAud.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVis.mean, lty=1, lwd=2)
lines(0:1, muAud.mean, lty=2, lwd=2)

#legend(0.1, 65, legend=c("Visual", "Auditory"), lty=1:2, cex=0.8)

# 4th panel biimodal Power, participant GP

auditory=c(0,1) # auditory high/low
visual=0 # visual low
muVL <- link(GP_Pow_Bim2b, data=data.frame(aud_N=auditory, vis_N=visual))
muVL.mean <- apply(muVL, 2, mean)
muVL.HPDI <- apply(muVL, 2, HPDI, prob=0.95)

auditory=c(0,1) # auditory high/low
visual=1 # visual high
muVH <- link(GP_Pow_Bim2b, data=data.frame(aud_N=auditory, vis_N=visual))
muVH.mean <- apply(muVH, 2, mean)
muVH.HPDI <- apply(muVH, 2, HPDI, prob=0.95)

plot(0, 0, xlim = c(-0.5,1.5), ylim = c(0,70), main="Bimodal", xlab="Auditory component intensity", ylab="Power", type = "n", xaxt="n", bty="l")
with(GP_d_bi, points(Pow ~ aud_N))
axis(1, at=c(0,1), labels = c("Low", "High"))

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVL.HPDI[2,1], muVL.HPDI[2,2],muVL.HPDI[1,2],muVL.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

x<-c(0.02,0.99,0.99,0.02)
y<-c(muVH.HPDI[2,1], muVH.HPDI[2,2],muVH.HPDI[1,2],muVH.HPDI[1,1])
polygon(x,y,col=gray(0.9), border = NA)

lines(0:1, muVL.mean, lty=1, lwd=2)
lines(0:1, muVH.mean, lty=2, lwd=2)

#legend(0.1, 65, legend=c("Visual - Low", "Visual - High"), lty=1:2, cex=0.8)

mtext("Participant GP", outer = TRUE, line=-23, cex = 1.5)
