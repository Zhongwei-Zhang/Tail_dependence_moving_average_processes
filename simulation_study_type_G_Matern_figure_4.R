library(INLA)
library(fields)
library(ngme2)
library(statmod)
library(parallel)

#set your working directory
setwd("your working directory")
#load Rcpp functions
Rcpp::sourceCpp("Code/functions.cpp")

##################################### Convergence of chi_u #############################################

################################################
################ NIG noise #############
# specify SPDE parameters
kappa <- 2
nu <- 1
alpha <- 2
# specify NIG noise: mu + gamma_nig*Gamma + sigma*sqrt(Gamma)*Z, Gamma~IG(nu_nig,nu_nig), Z~N(0,1)
sigma <- 1  #variance parameter
nu.nig <- 1
gamma.nig <- 0 #skewness parameter
mu.nig <- -gamma.nig  #location parameter

# generate observational locations
nx <- 10
x.loc <- y.loc <- seq(0.1, 1, length.out=nx)
coord <- as.matrix(expand.grid(x.loc,y.loc))
set.seed(12345)
coord <- coord - matrix(runif(200,0,0.1), nrow=nx^2)
n.obs <- dim(coord)[1]
plot(coord)
distmat <- rdist(coord)

# construct a lattice covering the observational loctions
nx.lattice <- 40
lattice <- inla.mesh.lattice(seq(-0.5, 1.5, length.out = nx.lattice), seq(-0.5, 1.5, length.out = nx.lattice))
# construct mesh based on the lattice with only outer extension
mesh <- inla.mesh.create(lattice = lattice, refine = list(max.edge = 5/(nx.lattice-1)*sqrt(2)+0.1))
mesh$n
plot(mesh)
points(mesh$loc, pch=19, cex=.5, col="red")
points(coord, pch=23, col="darkgreen")

# compute the finite emelent matrices based on the mesh
fem <- inla.mesh.fem(mesh)
h <- diag(fem$c0)
C <- fem$c1
G <- fem$g1
K <- kappa^2 * C + G

# simulation of NIG noise model
n.sample <- 1*10^5
# extract projection matrix
A <- inla.spde.make.A(mesh, coord)


######### simulation of data 
for (i in 1:10){
   set.seed(12347+i) #start from 12345
   V <- t(mapply(rinvgauss, shape=nu.nig*(h)^2,mean=h,n=n.sample))
   Z <- matrix(rnorm(n=mesh$n*n.sample), ncol=n.sample)
   Noise.nig <- mu.nig*h + gamma.nig*V + sigma*sqrt(V)*Z
    # obtain simulations at the observational sites
   X <- t(as.matrix(A%*%solve(K, Noise.nig)))
   save(X, file=paste("Data/Simu_NIG_",i,".RData",sep=""))
   rm(V,Z,Noise.nig,X)
}
#######################

Coeff.matrix <- A%*%solve(K)
Coeff.matrix.norm <- as.matrix(Coeff.matrix/apply(Coeff.matrix,1,max))
eta.fem <- eta_matrix_cpp(Coeff.matrix.norm)

#extract the simulated data at certain locations
locations <- c(1, 2, 21, 53, 97)
load("Data/Simu_NIG_1.RData")
X.sub <- X[,locations]
for (i in 2:10){
  load(paste("Data/Simu_NIG_",i,".RData",sep=""))
  X.sub <- rbind(X.sub, X[,locations])
}

#extremal dependence
u <- seq(0.8, 0.9999, by=0.0005)
emp_chi_loc <- array(NA, dim=c(length(u),rep(length(locations),2)))
for (i in 1:length(u)){
    emp_chi_loc[i,,] <- chieta_matrix(X.sub, u=u[i])$chi
}
save(emp_chi_loc, file="Data/emp_chi_locations.RData")
#load("Data/emp_chi_locations.RData")

# plot
u.start <- 301
a <- 1
colorpalette <- grey.colors(3, start=0.3, end=0.7)
png("Figure/emp_chi_wrt_u.png", units="in", height = 4, width = 8, res=300)
par(mfrow=c(1,2), pty="s")
par(mar = c(3, 0, 0.1, 0))
plot(mesh, main="")
points(x=coord[, 1][1], y=coord[, 2][1], pch=19, cex=1, col=1)
points(x=coord[, 1][locations[-1]], y=coord[, 2][locations[-1]], pch=17, cex=1, col=1)

par(mar = c(3, 3, 0.3, 0.2))
plot(u[u.start:400], emp_chi_loc[,1,2][u.start:400], xlab="", ylim=c(0,1),
     ylab="", type="l", lwd=1.5, col=1, cex.lab=a, cex.axis=a/1.1)
title(ylab=expression(chi(q)), mgp=c(2,1,0), cex.lab=a)
title(xlab="q", mgp=c(2,1,0), cex.lab=a)
for(i in 1:(length(locations)-2)){
  lines(u[u.start:400], emp_chi_loc[,1,i+2][u.start:400], col=colorpalette[i], 
        lwd=1.5, cex.lab=a, lty=i+1)
}
dev.off()




