######################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Generates predictions for SPDE1 (independent
##          hi-res SPDE models on S2).
##          This code comes with no warranty or guarantee of any kind.
######################################################################

## Load packages and functions
library("dplyr")
library("INLA")
source("utils.R")

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load some details we have on the fine mesh
load(file.path(cache_folder, "data_for_SST2M_no_chunks.rda"))
load(file.path(cache_folder, "C1_matrices.rda"))

## Load trainnig data and validation data
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))
y_pred <- read.csv(file.path(cache_folder, "y_pred.csv"))

## Do SPDE meshes based on first tiling that we use in SPDE2
fnames <- file.path(cache_folder, dir(cache_folder, pattern = "colouring1.rda"))
nchunks <- length(fnames)
true.radius.of.earth = 6371
cutoff <- 100

## Initialise predictions to NA
y_pred$spde_hi_pred <- NA
y_pred$spde_hi_se_obs <- NA

## Prior on range of sd as in SPDE2
meanorig <- 3
paropt <- optim(par = c(log(meanorig), 0.1), fn = lognorm_cost, 
                lx = 0.1, ux = 6)$par
mu_logsigma <-paropt[1]
sd_logsigma <- paropt[2]

## Range of Y1 as in SPDE2
meanorig <- 200/6371
paropt <- optim(par = c(log(meanorig), 0.1), fn = lognorm_cost, 
                lx = 30/6371, ux = 900/6371)$par
mu_logrho1 <-paropt[1]
sd_logrho1 <- paropt[2]

## Convert to tau/kappa
log_kappa1 <- log_rho2log_kappa(mu_logrho1, 1)
log_tau1 <- log_sigma2log_tau(mu_logsigma, log_kappa1, 1)

## We will collect the results from each block in a list
y_pred_list <- list()

## For every tile/chunk
for(i in 1:nchunks) {

    ## Extract numbers from filename -- could either be first digit or 
    ## first two digits
    str <- gsub("[^0-9.]", "",  fnames[i])
    if(nchar(str) == 4) {
        j <- as.numeric(substr(str, 1, 2))
    } else {
        j <- as.numeric(substr(str, 1, 1))
    }
    
    ## Load data for this chunk
    print(paste0("INLA fine res: Doing chunk ", j))
    load(fnames[i])
    ch <- chunk
    
    ## See which training data and validation data we can associate with
    ## thi tile
    sub_mesh <- filter(points_fine, class1 == j)
    C1_sub <- C1[, sub_mesh$id]
    C1_pred_sub <- C1_pred[, sub_mesh$id]
    idx_in <- which(rowSums(C1_sub) >= 0.5)
    idx_pred_in <- which(rowSums(C1_pred_sub) >= 0.5)
    yin <- y_tot[idx_in, ]
    ypredin <- y_pred[idx_pred_in, ]

    ## Create a new mesh for these data
    newmesh = inla.mesh.create(loc = lonlat3D(c(yin$lon, ypredin$lon),
                                              c(yin$lat, ypredin$lat)),
                               cutoff = cutoff/true.radius.of.earth,
                               refine=list(max.edge = 30/true.radius.of.earth,
                                           extend = list(offset = -0.5)))

    ## Create an SPDE on this mesh using our priors
    new_spde <- inla.spde2.matern(mesh = newmesh,
                                  alpha = 2,
                                  B.tau = cbind(log_tau1, -1, 1),
                                  B.kappa = cbind(log_kappa1, 0, -1),
                                  theta.prior.mean = c(0, 0),
                                  theta.prior.prec = c(1/sd_logsigma^2, 1/sd_logrho1^2))

    ## Generate observation matrices for this new mesh
    C1_new <- inla.spde.make.A(newmesh, lonlat3D(yin$lon, yin$lat))
    C1_pred_new <- inla.spde.make.A(newmesh, lonlat3D(ypredin$lon, ypredin$lat))

    ## INLA is much faster if we do not estimate the mean, we therefore 
    ## detrend manually -- OK since we are in a big data scenario
    this_mean <- mean(yin$z)
    
    ## Use INLA to fit SPDE
    proc_time_INLA_hi <- system.time({
          INLA_output_hi <- INLA_spatial_estimation2(new_spde,
                                                     yin$z - this_mean,
                                                     C1_new)})
    ## Save INLA output
    save(INLA_output_hi, newmesh, proc_time_INLA_hi,
         file = paste0(file.path(cache_folder, paste0("SPDE1-indep", i, "_raw.rda"))))
    
    ## Project field onto validation locations
    project <- inla.mesh.projector(newmesh,
                                   loc = lonlat3D(ypredin$lon, ypredin$lat))

    ## Predictions
    y_pred$spde_hi_pred[idx_pred_in] <-
        inla.mesh.project(project, INLA_output_hi$summary.random$spatial.field$mean +
                                   this_mean)
    
    ## Prediction standard errors
    y_pred$spde_hi_se_obs[idx_pred_in] <-
        inla.mesh.project(project, sqrt(INLA_output_hi$summary.random$spatial.field$sd^2 +
                                        1/INLA_output_hi$summary.hyperpar[1,1]))

    ## Return results as a list
    y_pred_list[[i]] <- list(y_pred_sub = y_pred, idx = idx_pred_in)
 }

## Compile results into a data framre
for(i in 1:nchunks) {
    idx <- y_pred_list[[i]]$idx
    ylist <- y_pred_list[[i]]
    y_pred$spde_hi_pred[idx] <- ylist$y_pred_sub$spde_hi_pred[idx]
    y_pred$spde_hi_se_obs[idx] <- ylist$y_pred_sub$spde_hi_se_obs[idx]
}

## Save results to disk
save(y_pred, file = file.path(cache_folder, "SPDE1-indep_output.rda"))
