library(INLA)
library(inlabru)
library(fields)
library(ngme2)
library(ggplot2)

# set working directory
setwd("your working directory")
# load Rcpp functions
Rcpp::sourceCpp("functions.cpp")

data(argo_float)
head(argo_float)
# take longitude and latitude to build the mesh
max.edge    <- 1
bound.outer <- 5
loc_2d <- unique(cbind(argo_float$lon, argo_float$lat))
distmat <- rdist(loc_2d)
# nrow(loc) == nrow(dat) no replicate
argo_mesh <- inla.mesh.2d(loc = loc_2d,
                          # the inner edge and outer edge
                          max.edge = c(1,5),
                          cutoff = 0.3,
                          # offset extension distance inner and outer extenstion
                          offset = c(max.edge, bound.outer)
)
plot(argo_mesh)
argo_mesh$n

plot(argo_mesh, main="", col="lightgrey")
points(loc_2d[, 1], y = loc_2d[, 2], pch=19, cex=.8, col="darkgreen")


##################################### compute residual tail dependece
# SPDE parameters
d <- 2
kappa <- 2
# Compute the finite emelent matrices based on the mesh
fem <- inla.mesh.fem(argo_mesh)
h <- diag(fem$c0)
C <- fem$c1
G <- fem$g1
K <- kappa^2 * C + G

library(expm)
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
  if (x==0 & alpha>2){
    return(gamma((alpha-d)/2)/(4*pi*gamma(alpha/2)*kappa^(alpha-d)))
  } else {
  const <- 2^{1-(alpha-d)/2}/(4*pi*gamma(alpha/2)*kappa^(alpha-d))
  return(const*(kappa*x)^((alpha-d)/2)*besselK(kappa*x, nu=(alpha-d)/2))
  }
}
G0 <- function(alpha){
  gamma((alpha-d)/2)/(4*pi*gamma(alpha/2)*kappa^(alpha-d)) #finite if alpha>2
}

# Check approximation accuracy by eta (Type G model) 
M <- list(M.alpha1, M.alpha2, M.alpha3, M.alpha4, M.alpha5, M.alpha6)
eta.IA <- eta.fem <- list()
#SPDE projection matrix
A <- inla.spde.make.A(argo_mesh, loc_2d)

alpha = 3
# Integral approximation
dist.obs.meshnode <- rdist(x1=loc_2d, x2=argo_mesh$loc[,1:2])
A.coeff <- structure(vapply(dist.obs.meshnode, G, alpha=alpha, numeric(1)), 
                     dim=dim(dist.obs.meshnode))
A.coeff.norm <- A.coeff/apply(A.coeff, 1, max)
eta.IA[[alpha]] <- eta_matrix_cpp(A.coeff.norm)
# FEM approximation
Coeff.matrix <- A%*%M[[alpha]]
Coeff.matrix.norm <- as.matrix(Coeff.matrix/apply(Coeff.matrix,1,max))
eta.fem[[alpha]] <- eta_matrix_cpp(Coeff.matrix.norm)

#save(eta.IA, eta.fem, file="Data/Intro_eta.RData")

#plot
t <- seq(0.001, 5, by=0.01)
#only consider distances less than or equal to 5
index <- distmat[upper.tri(distmat)] <= 5
colorpalette <- grey.colors(12, start=0.3, end=0.9)

png("Figure/Intro_argo.png", units="in", height = 4, width = 8, res=300)
par(mfrow=c(1,2), pty="s")
par(mar = c(3, 0, 0.1, 0))
plot(argo_mesh, main="")
points(x=loc_2d[, 1], y=loc_2d[, 2], pch=17, cex=1, col=1)

par(mar = c(3, 3, 0.3, 0.2))
plot(distmat[upper.tri(distmat)][index], eta.IA[[alpha]][upper.tri(eta.IA[[alpha]])][index], 
     xlab="", ylim=c(0.5,1), ylab="", pch=19, col=colorpalette[4])
title(ylab=expression(eta), mgp=c(2,1,0))
title(xlab="h", mgp=c(2,1,0))
points(distmat[upper.tri(distmat)][index], eta.fem[[alpha]][upper.tri(eta.fem[[alpha]])][index], 
       pch=17, col=alpha(colorpalette[10], 0.6))
if(alpha>2){
    lines(t, 1/2+sapply(t, G, alpha=alpha)/(2*G0(alpha=alpha)), lwd=2.5, col=1)
}
dev.off()

