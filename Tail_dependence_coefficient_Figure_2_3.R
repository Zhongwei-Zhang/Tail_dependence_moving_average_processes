library(ghyp)
library(Bessel)

#set your working directory
setwd("your working directory")

# Model: X_1 = Y_1 + 0.3 Y_2, X_2 = Y_1 + a_{22} Y_2, with a_{22} > 0.3
a_12 <- 0.3
a_22 <- seq(a_12+0.001, 0.99999, length.out=500)

temp_lambda <- c(-0.5, 1, 5, 30) #different lambda
chi.a22.lambda <- matrix(NA, nrow=length(temp_lambda), ncol=length(a_22))
temp_chi_psi <- c(0.5, 1, 5, 30) #different chi, psi
chi.a22.chi <- matrix(NA, nrow=length(temp_chi_psi), ncol=length(a_22))
chi.a22.psi <- matrix(NA, nrow=length(temp_chi_psi), ncol=length(a_22))

#run the following code for different values of lambda, chi, and psi respectively
for (i in 1:length(temp_lambda)){
  lambda <- temp_lambda[i]
  psi <- 1
  mu <- gammaGH <- 0
  chi <- 1
  
  #moment generating function
  MGF_GH <- function(t){
    exp(mu*t)*(psi/(psi-2*gammaGH*t-t^2))^(lambda/2)*
      besselK(sqrt(chi*(psi-2*gammaGH*t-t^2)), nu=lambda)/besselK(sqrt(chi*psi),nu=lambda)
  }

  # useful constants
  c_1 <- log(MGF_GH(a_12*sqrt(psi)))/sqrt(psi)
  c_2 <- log(MGF_GH(a_22*sqrt(psi)))/sqrt(psi)
  c_3 <- (c_2 - c_1)/(a_22 - a_12)
  # tail dependence coefficient as a function of a_22
  ghyp.margin <- ghyp(lambda=lambda,chi=chi,psi=psi,mu=mu,gamma=gammaGH)
  TDC <- function(a){
    c3_a <- (log(MGF_GH(a*sqrt(psi)))/sqrt(psi) - c_1)/(a - a_12)
    const1 <- integrate(function(x) exp(a*sqrt(psi)*x)*dghyp(x,object = ghyp.margin), lower=-Inf, upper=c3_a)$value
    const2 <- integrate(function(x) exp(a_12*sqrt(psi)*x+dghyp(x,object = ghyp.margin,logvalue=TRUE)), lower=c3_a, upper=Inf)$value
    const1/MGF_GH(a*sqrt(psi)) + const2/MGF_GH(a_12*sqrt(psi))
  }

  chi.a22.lambda[i,] <- sapply(a_22, TDC)
}

colorpalette <- grey.colors(3, start=0.3, end=0.7)
a = 2.5

png("Figure/limitchi_a22.png", units="in", width=12, height=4, res=300)
par(mfrow=c(1,3), pty="s", mar = c(4, 5, 0.3, 0.3))
plot(a_22, chi.a22.lambda[1,], type="l", lty=1, ylab="",
     xlab="", ylim=c(0,1), cex.lab=a, cex.axis=a/1.1)
title(ylab=expression(chi), mgp=c(2.5,1,0), cex.lab=a)
title(xlab=expression(a[22]), mgp=c(3,1,0), cex.lab=a)
lines(a_22, chi.a22.lambda[2,], lty=2, col=colorpalette[1])
lines(a_22, chi.a22.lambda[3,], lty=3, col=colorpalette[2])
lines(a_22, chi.a22.lambda[4,], lty=4, col=colorpalette[3])

plot(a_22, chi.a22.chi[1,], type="l", lty=1, ylab="", yaxt='n',
     xlab="",ylim=c(0,1),cex.lab=a,cex.axis=a/1.1)
title(xlab=expression(a[22]), mgp=c(3,1,0), cex.lab=a)
lines(a_22, chi.a22.chi[2,], lty=2, col=colorpalette[1])
lines(a_22, chi.a22.chi[3,], lty=3, col=colorpalette[2])
lines(a_22, chi.a22.chi[4,], lty=4, col=colorpalette[3])

plot(a_22, chi.a22.psi[1,], type="l", ylab="", yaxt='n',
     xlab="",ylim=c(0,1),cex.lab=a,cex.axis=a/1.1)
title(xlab=expression(a[22]), mgp=c(3,1,0), cex.lab=a)
lines(a_22, chi.a22.psi[2,], lty=2, col=colorpalette[1])
lines(a_22, chi.a22.psi[3,], lty=3, col=colorpalette[2])
lines(a_22, chi.a22.psi[4,], lty=4, col=colorpalette[3])
dev.off()


########################################################################
# plot convergence of eta for non-Gaussian OU process
Time <- 4
a <- 0.2
delta <- c(0.4, 0.2, 0.05) # in decreasing order
N <- 1000
t <- seq(0.01, Time, length=N)
eta_t <- matrix(NA, nrow= length(delta), ncol=N)
for(j in 1:length(delta)){
  partition.point <- seq(0, Time, by=delta[j])
  for (i in 1:N){
    t_n2 <- max(partition.point[partition.point-t[i] <=0])
    eta_t[j, i] <- 1/(2-exp(-a*t_n2))
  }
}
eta_true <- 1/(2-exp(-a*(t)))

colorpalette <- grey.colors(3, start=0.3, end=0.7)
c <- 1.8 #scale parameter for axis lab

png("Figure/limiteta.png", units="in", width=6, height=6, res=300)
par(pty="s", mar=c(2.8, 4.5, 0.2, 0.2))
plot(t, eta_true, ylim=c(0.64,1), xlab = "", yaxt="n", lwd=1.5,
     ylab="", type="l", col=1, cex.lab=c, cex.axis=c/1.1)
axis(2, at=c(0.7,0.8,0.9,1),labels=c(0.7,0.8,0.9,1), las=2, cex.axis=c/1.1)
title(ylab=expression(eta), mgp=c(3.2,1,0), cex.lab=c)
title(xlab=expression(t), mgp=c(2.6,1,0), cex.lab=c)
lines(t, eta_t[1,], col=colorpalette[3], lty=4, lwd=1.5)
lines(t, eta_t[2,], col=colorpalette[1], lty=3, lwd=1.5)
lines(t, eta_t[3,], col=colorpalette[2], lty=2, lwd=1.5)
dev.off()


