library(INLA)
library(fields)
library(scales)
library(expm)

#set your working directory
setwd("~/Dropbox KAUST/Projects/NonGaussian_SPDE/R code")
#setwd("C:/Users/zhangzho/Dropbox (KAUST)/Projects/NonGaussian_SPDE/R code")
Rcpp::sourceCpp("Code/functions.cpp")

# generate observational locations
nx <- 15
x.loc <- y.loc <- seq(0.1, 1.5, length.out=nx)
coord <- as.matrix(expand.grid(x.loc,y.loc))
set.seed(123)
coord <- coord - matrix(runif(nx^2*2,0,0.1), nrow=nx^2)
n.obs <- dim(coord)[1]
plot(coord)
distmat <- rdist(coord)

# construct a lattice covering the observational loctions
nx.lattice <- 40
lattice <- inla.mesh.lattice(seq(-1, 2.5, length.out = nx.lattice), seq(-1, 2.5, length.out = nx.lattice))
# construct mesh based on the lattice with only outer extension
mesh <- inla.mesh.create(lattice = lattice, refine = list(max.edge = 5/(nx.lattice-1)*sqrt(2)+0.1))
mesh$n
mesh$loc
plot(mesh)
points(mesh$loc, pch=19, cex=.5, col="red")
points(coord, pch=23, cex=.7, col="darkgreen")

# SPDE parameters
#alpha defined later
d <- 2
kappa <- 2
# Compute the finite emelent matrices based on the mesh
fem <- inla.mesh.fem(mesh)
h <- diag(fem$c0)
C <- fem$c1
G <- fem$g1
K <- kappa^2 * C + G

# computing the principle square root of a matrix is slow
M.alpha1 <- expm::sqrtm(solve(K))%*%diag(1/sqrt(diag(C)), dim(C)[1], dim(C)[1])
M.alpha2 <- solve(K)
temp.mat <- M.alpha2%*%C
M.alpha3 <- temp.mat%*%M.alpha1
M.alpha4 <- temp.mat%*%M.alpha2
M.alpha5 <- temp.mat%*%M.alpha3
M.alpha6 <- temp.mat%*%M.alpha4

# Green's function
G <- function(x, alpha){
  const <- 2^{1-(alpha-d)/2}/(4*pi*gamma(alpha/2)*kappa^(alpha-d))
  return(const*(kappa*x)^((alpha-d)/2)*besselK(kappa*x, nu=(alpha-d)/2))
}
G0 <- function(alpha){
  gamma((alpha-d)/2)/(4*pi*gamma(alpha/2)*kappa^(alpha-d)) #finite if alpha>2
}

####################### Check approximation accuracy by correlation (Gaussian model) 
# compute correlation matrix for approximation models
A <- inla.spde.make.A(mesh, coord)
Coeff.matrix <- A%*%M.alpha5
Sigma <- Coeff.matrix%*%C%*%t(Coeff.matrix)
fem.cor <- as.matrix(cov2cor(Sigma)) #correlation of FEM approximations
# correlation for integral approximation
dist.obs.meshnode <- rdist(x1=coord, x2=mesh$loc[,1:2])
#A_coeff <- apply(dist_obs_meshnode, 1:2, G)
A.coeff <- structure(vapply(dist.obs.meshnode, G, alpha=5, numeric(1)), 
                     dim=dim(dist.obs.meshnode))
integral.cor <- cov2cor(A.coeff%*%t(A.coeff))

#correlations of true process
t <- seq(0.001, 2, by=0.001)
matern.cor <- inla.matern.cov(x=t, kappa=kappa, nu=5-d/2, d=2, corr=TRUE)

#plot
index <- order(distmat[upper.tri(distmat)])
plot(distmat[upper.tri(distmat)], integral.cor[upper.tri(integral.cor)], xlab="h", 
     ylim=c(0,1), ylab=expression(rho), pch=19, col="grey", 
     main=expression(alpha == 2))

points(distmat[upper.tri(distmat)], fem.cor[upper.tri(fem.cor)], pch=18, 
       col="darkgreen")
lines(t, matern.cor, col="red", lwd=4)


####################### Check approximation accuracy by eta (Type G model) 
M <- list(M.alpha1, M.alpha2, M.alpha3, M.alpha4, M.alpha5, M.alpha6)
eta.IA <- eta.fem <- list()
#SPDE projection matrix
A <- inla.spde.make.A(mesh, coord)

for(alpha in 1:6){
  # Integral approximation
  dist.obs.meshnode <- rdist(x1=coord, x2=mesh$loc[,1:2])
  A.coeff <- structure(vapply(dist.obs.meshnode, G, alpha=alpha, numeric(1)), 
                       dim=dim(dist.obs.meshnode))
  A.coeff.norm <- A.coeff/apply(A.coeff, 1, max)
  eta.IA[[alpha]] <- eta_matrix_cpp(A.coeff.norm)
  # FEM approximation
  Coeff.matrix <- A%*%M[[alpha]]
  Coeff.matrix.norm <- as.matrix(Coeff.matrix/apply(Coeff.matrix,1,max))
  eta.fem[[alpha]] <- eta_matrix_cpp(Coeff.matrix.norm)
}

save(eta.IA, eta.fem, file="Data/TwoApprox_eta_mesh40_allalpha.RData")

load("Data/TwoApprox_eta_mesh40_allalpha.RData")

#plot
t <- seq(0.001, 2.5, by=0.001)
colorpalette <- grey.colors(12, start=0.3, end=0.9)

png("Figure/Conv_eta_typeG_allalpha.png", units="in", height = 5, width = 20, res=300)
par(mfrow=c(1,4), pty="s")
alpha = 2
a = 2
par(mar = c(5, 5, 3, 0))
plot(distmat[upper.tri(distmat)], eta.IA[[alpha]][upper.tri(eta.IA[[alpha]])], xlab="", 
     ylim=c(0.5,1), ylab="", pch=19, col=colorpalette[4],
     main=bquote(alpha == .(alpha)), 
     cex.lab=a*1.1, cex.axis=a, cex.main=a, cex.sub=a
     )
title(ylab=expression(eta), mgp=c(3.5,1,0), cex.lab=a*1.1)
title(xlab="h", mgp=c(4,1,0), cex.lab=a*1.1)
points(distmat[upper.tri(distmat)], eta.fem[[alpha]][upper.tri(eta.fem[[alpha]])], 
       pch=17, col=alpha(colorpalette[10], 0.6))
lines(t, rep(1/2, length(t)), lwd=a, col=1)

for (alpha in 3:5){
  par(mar = c(5, 4, 3, 0))
  plot(distmat[upper.tri(distmat)], eta.IA[[alpha]][upper.tri(eta.IA[[alpha]])], xlab="", 
       ylim=c(0.5,1), ylab="", pch=19, col=colorpalette[4],
       main=bquote(alpha == .(alpha)), yaxt='n',
       cex.lab=a*1.1, cex.axis=a, cex.main=a, cex.sub=a
       )
       #main=expression(alpha == 4 ~Integral~Approx.))
  title(xlab="h", mgp=c(4,1,0), cex.lab=a*1.1)
  points(distmat[upper.tri(distmat)], eta.fem[[alpha]][upper.tri(eta.fem[[alpha]])], 
         pch=17, col=alpha(colorpalette[10], 0.6))
  if(alpha>2){
    lines(t, 1/2+sapply(t, G, alpha=alpha)/(2*G0(alpha=alpha)), lwd=a, col=1)
    lines(t, sapply(t/2, G, alpha=alpha)/G0(alpha), lty=2, lwd=a, col=1)
  }
}
dev.off()


