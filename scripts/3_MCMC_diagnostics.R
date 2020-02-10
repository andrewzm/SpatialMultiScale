###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Convergence analysis of Gibbs-sampler traces
##          This code comes with no warranty or guarantee of any kind.
###########################################################################
library("dplyr")
library("INLA")
library("parallel")
source("utils.R")

set.seed(1)

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load meshes
load(file.path(cache_folder, "meshes.rda"))

## Extract variable sizes
load(file.path(cache_folder, "etasamp0_2"))
r_eta0 <- length(current_samp)

load(file.path(cache_folder, "etasamp1_2"))
r_eta1 <- length(current_xsamp1)

load(file.path(cache_folder, "hypsamps_MCMC.rda"))
ntheta <- ncol(hyp_samp_log_rho1)

## Raw number of samples
N_raw <- nrow(hyp_samp_log_rho1)
if(!(N_raw  == 10000)) warning("Expecting 10000 iterations. Has the MCMC algorithm
                        finished?")

## Iterations to use: Discard first 5000 as burn-in, and thin by a factor of 50
samp_num <- seq(round(N_raw / 2), N_raw, by = 50)
N_thinned <- length(samp_num)

## Compile eta0 and eta1 into matrices
eta0 <- matrix(0, N_thinned, r_eta0)
eta1 <- matrix(0, N_thinned, r_eta1)
for(i in seq_along(samp_num)) {
  
  ## Load the samples
  load(file.path(cache_folder, paste0("etasamp0_", samp_num[i])))
  load(file.path(cache_folder, paste0("etasamp1_", samp_num[i])))
  
  ## Update the matrix with the thinned samples
  eta0[i, ] <- current_samp
  eta1[i, ] <- current_xsamp1
}

## Compile parameter samples into matrices
load(file.path(cache_folder, "hypsamps_MCMC.rda"))
logsigma0 <- matrix(hyp_samp0[1, ])
logrho0 <- matrix(hyp_samp0[2, ])
logsigma1 <- hyp_samp_log_sigma1
logrho1 <- hyp_samp_log_rho1
logsigma_eps <- hyp_samp_log_sigma_eps


## Plot all the trace plots. Most appear OK.
plot_par_traces <- function(obj, fname, panels = c(1, 1), ylab = NULL) {
  png(paste0("../img/", fname), width = panels[2] * 750, 
      height = panels[1] * 750, res = 150)
  par(mfrow = panels)
  ess <- rep(0, ntheta)
  for(i in 1 : ncol(obj))  
    plot(obj[, i], type = 'l', main = i, xlab = "Itreation", ylab = ylab)
  dev.off()
}

plot_par_traces(logsigma0[, , drop = FALSE], "traces_logsigma0.png",
                ylab = expression(log(sigma[0])))
plot_par_traces(logrho0[, , drop = FALSE], "traces_logrho0.png",
                ylab = expression(log(rho[0])))
plot_par_traces(logsigma_eps, "traces_logsigma_eps.png",
                panels = c(21, 10), 
                ylab = expression(bold(theta)[epsilon]))
plot_par_traces(logsigma1, "traces_logsigma1.png",
                panels = c(21, 10),
                ylab = expression(bold(theta)[sigma[1]]))
plot_par_traces(logrho1, "traces_logrho1.png",
                panels = c(21, 10),
                ylab = expression(bold(theta)[rho[1]]))

## Compute neff of all unknowns. Most are OK, some parameters have low
## effective sample sizes but shouldn't be of too much concern
eta0_ess <- compute_neff(eta0)
eta1_ess <- compute_neff(eta1)
logsigma0_ess <- compute_neff(logsigma0[samp_num, , drop = FALSE])
logrho0_ess <- compute_neff(logrho0[samp_num, , drop = FALSE])
logsigma1_ess <- compute_neff(logsigma1[samp_num, ])
logrho1_ess <- compute_neff(logrho1[samp_num, ])
logsigma_eps_ess <- compute_neff(logsigma_eps[samp_num, ])

cat("Median effective sample sizes: \n")
cat("eta0: "); cat(median(eta0_ess)); cat("\n")
cat("eta1: "); cat(median(eta1_ess)); cat("\n")
cat("logsigma0: "); cat(logsigma0_ess); cat("\n")
cat("logrho0: "); cat(logrho0_ess); cat("\n")
cat("logsigma1: "); cat(median(logsigma1_ess)); cat("\n")
cat("logrho1: "); cat(median(logrho1_ess)); cat("\n")
cat("logsigma_eps: "); cat(median(logsigma_eps_ess)); cat("\n")

## Now we compute the effective sample size of the process Y, which is our
## key quantity of interest when doing prediction.

## First form the base map
base_map0 <- plot_map(mesh, current_samp, plot = FALSE,
                        xlim = NULL, ylim =  c(-85, 90))
base_map1 <- plot_map(mesh_fine, current_xsamp1, plot = FALSE,
                        xlim = NULL, ylim =  c(-85, 90))
  
## Extract projections
proj0 <- base_map0$proj
proj1 <- base_map1$proj
  
## Find global field for each sample on the base maps
ngridpts <- nrow(base_map1$grid_data)
zsamp0 <- zsamp1 <- zsamp <- 
    matrix(0, N_thinned, ngridpts)
  
for (i in 1:N_thinned) {
  zsamp_mat0 <- inla.mesh.project(proj0, field = eta0[i, ])
  zsamp_mat1 <- inla.mesh.project(proj1, field = eta1[i, ])
  zsamp_mat <- zsamp_mat0 + zsamp_mat1
  
  zsamp0[i, ] <- reshape2::melt(zsamp_mat0)[,3]
  zsamp1[i, ] <- reshape2::melt(zsamp_mat1)[,3]
  zsamp[i, ] <- reshape2::melt(zsamp_mat)[,3]
}

## Remove all the land pixels
land_pixels <- which(apply(zsamp, 2, function(x) any(is.na(x))))
zsamp_sea <- zsamp[, -land_pixels]

## Compute the effective sample size of Y at these ~350000 ocean locations
zsamp_sea_ess <- compute_neff(zsamp_sea)

## Compute the empirical densities of effective sample sizes from purely
## random samples and AR1 traces
N_Sim <- 50000
Random_Matrix <- matrix(rnorm(N_Sim * N_thinned), N_thinned, N_Sim)
ess_empirical_Random <-  Random_Matrix %>% compute_neff()

a <- 0.1
AR1_Matrix <- matrix(0, N_thinned, N_Sim)
AR1_Matrix[1, ] <- rnorm(N_Sim, sd = sqrt(1 / (1 - a^2)))
for(i in 2:N_thinned) AR1_Matrix[i, ] <- a * AR1_Matrix[i - 1, ] + rnorm(N_Sim)
ess_empirical_AR1 <-  AR1_Matrix %>% compute_neff()

## Plot the empirical density of the effective sample size, as well as 
## those from random and AR1 traces.
png("../img/neff.png", width = 1200, height = 900, res = 200)
par(mar = c(4, 5, 1.5, 2) + 0.1)
plot(density(ess_empirical_Random, width = 10), xlim = c(0, 300), main = "",
     xlab = expression(hat(n)[eff]), lwd = 2, ylim = c(0, 0.022))
lines(density(zsamp_sea_ess, width = 10), col = "red", lty = 2, lwd = 2)
lines(density(ess_empirical_AR1, width = 10), col = "blue", lty = 3, lwd = 2)
dev.off()

## Plot a couple of traces just to show lack of autocorrelation
par(mfrow = c(2,2))
for (i in sample(1:ncol(zsamp_sea), 4L)) 
  plot(zsamp_sea[, i], type = 'l', ylab = paste0("Y", i))

