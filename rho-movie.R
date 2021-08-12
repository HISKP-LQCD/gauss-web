library("animation")
library("plotrix")
library("png")

etmclogo <- readPNG("~/daten/workdir/papers/rho211/plots/Logo_ETMC_RGB.png")

r0sqZmpi <- function(Zchi, Cchi, g, Mpi) {
  Zchi + Cchi*Mpi^2 -g*sqrt(Zchi)*Mpi^3
}
getMGamma <- function(Z) {
  return(list(M=Re(sqrt(Z)), Gamma=2*Im(sqrt(Z))))
}
tandelta <- function(M, Gamma, Mpi, Ecm) {
  g <- sqrt(6*pi*M^2*Gamma/sqrt(M^2/4-Mpi^2)^3)
  return(g^2/6/pi*sqrt(Ecm^2/4-Mpi^2)^3/Ecm/(M^2-Ecm^2))
}

load("feng_2011_eq19_fit_summary.Rdata")
Zchi <- feng_2011_eq19_fit_summary$r0_sq_pole_chiral_re_val[3] + 1i*feng_2011_eq19_fit_summary$r0_sq_pole_chiral_im_val[3]
Cchi <- feng_2011_eq19_fit_summary$coupling_chi_val[3] 
gfit <- feng_2011_eq19_fit_summary$g_val[3]
X <- read.csv("rho_phaseshift_exp.csv")

logoing_func<-function(logo, x, y, size){
  dims<-dim(logo)[1:2] #number of x-y pixels for the logo (aspect ratio)
  AR<-dims[1]/dims[2]
  par(usr=c(0, 1, 0, 1))
  rasterImage(logo, x-(size/2), y-(AR*size/2), x+(size/2), y+(AR*size/2), interpolate=TRUE)
}

r0 <- 0.475
saveVideo({
  ani.options(interval=0.1)

  for(Mpi in seq(0.450, 0.135, by=-0.005)) {
    plotwitherror(x=X[,1], y=X[,2]/180*pi, dy=X[,3]/180*pi, pch=21,
                  col="red", bg="red",
                  xlim=c(0.5, 1.2), ylim=c(0,pi),
                  xlab=expression(E[CM] ~ "[GeV]"), ylab=expression(delta[1] ~ "[rad]"))
    MGamma <- getMGamma(Z=r0sqZmpi(Zchi=Zchi, Cchi=Cchi, g=gfit, Mpi=r0*Mpi/0.198))
    Ecm <- seq(0.5, 1.2, 0.01)
    td1 <- tandelta(MGamma$M/r0*0.198, MGamma$Gamma/r0*0.198, Mpi, Ecm)
    d1 <- atan(td1) %% pi
    lines(x=Ecm, y=d1, col="blue", lwd=2)
    text(x=0.6, y=1.475, labels=expression(M[pi]), pos=2, cex=1.5)
    text(x=0.6, y=1.5, labels=paste0("= ", floor(1000*Mpi), " MeV"), pos=4, cex=1.5)
    legend("topleft",
           legend=c("Experiment (Protopopescu et al, Phys. Rev. D7 , 1279 (1973))"), pch=21, col="red", pt.bg="red",
           bty="n")
    legend("topright",
           legend=c("Lattice QCD"), lty=1, lwd=2, col="blue",
           bty="n")
    logoing_func(etmclogo, x=0.9, y=0.2, size=0.2)
  }
}, img.name="rho-movie", loop=TRUE,
ani.width = 600, ani.height = 600, navigator=FALSE)
