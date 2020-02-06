###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Compiles results of SPDE2 in a data frame for comparison
##          Also examines parameter summaries etc.
##          This code comes with no warranty or guarantee of any kind.
###########################################################################

## Load packages and functions
library("INLA")
source("utils.R")

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


## Posterior means of sigma0 and rho0
sigma0_mean <- mean(exp(hyp_samp0[1, samp_num]))
rho0_mean <- mean(exp(hyp_samp0[2, samp_num]))

## Posterior IQRs of sigma1 and rho1
quantile(colMeans(exp(hyp_samp_log_sigma1[samp_num, ])), c(0.25, 0.75))
quantile(colMeans(exp(hyp_samp_log_rho1[samp_num, ])), c(0.25, 0.75))

## See which basis were removed from these mesh (these  need to be re-added
## for mesh projection onto validation locations using the INLA package)
S <- inla.spde.make.A(mesh_lo,
                      mesh_fine$loc)
S_symb <- S
S_symb@x <- S_symb@x / S_symb@x
rm_basis_idx <- which(colSums(S_symb) < 100)

## Projectors (parameters, Y0 and Y1)
y_pred <- read.csv(file.path(cache_folder, "y_pred.csv"))
proj_theta <- inla.mesh.projector(mesh_lo,
                               loc = lonlat3D(y_pred$lon, y_pred$lat))
proj_Y0 <- inla.mesh.projector(mesh, 
                               loc = lonlat3D(y_pred$lon,y_pred$lat))
proj_Y1 <- inla.mesh.projector(mesh_fine, 
                               loc = lonlat3D(y_pred$lon,y_pred$lat))

## Compile samples and get samples for the prediction distribution
## at the validation-data locations
for(i in seq_along(samp_num)) {
  
  load(file.path(cache_folder, paste0("etasamp0_", samp_num[i])))
  load(file.path(cache_folder, paste0("etasamp1_", samp_num[i])))
       
  if(i == 1) {
    etasamp0 <- matrix(0, 
                     length(current_samp), 
                     length(samp_num))
    etasamp1 <- matrix(0, 
                     length(current_xsamp1), 
                     length(samp_num))
    y_pred_samps <- matrix(0, 
                           nrow(y_pred), 
                           length(samp_num))
   sd_obs_samps <- matrix(0,
                           nrow(y_pred),
                           length(samp_num))
  }
  
  etasamp0[, i] <- current_samp
  etasamp1[, i] <- current_xsamp1
 
  ## Sample of measurement-error standard-deviation at the validation locations
  sd_obs_samps[, i] <- exp(inla.mesh.project(proj_theta,
                                            field = insert_elements(hyp_samp_log_sigma_eps[samp_num[i], ],
                                                                    rm_basis_idx,
                                                                    rep(0, length(rm_basis_idx)))))
  
  ## Process (not data) sample at validation locations
  y_pred_samps[, i] <- inla.mesh.project(proj_Y0, field = etasamp0[,i]) +
                       inla.mesh.project(proj_Y1, field = etasamp1[,i])
}

## Form prediction and prediction variance at validation-data locations
y_pred$TwoGMRF_pred <- rowMeans(y_pred_samps)
y_pred$TwoGMRF_se_obs <- sqrt(apply(y_pred_samps, 1, var) +
                                rowMeans(sd_obs_samps^2))

## Save data
save(y_pred, file = file.path(cache_folder, "SPDE2_output.rda"))
