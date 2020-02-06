######################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Generates predictions for SPDE0 (low-res SPDE model on S2).
##          This code comes with no warranty or guarantee of any kind.
######################################################################

## Load packages
library("INLA")
source("utils.R")

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load meshes
load(file.path(cache_folder, "meshes.rda"))

## Load training and validation data
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))
y_pred <- read.csv(file.path(cache_folder, "y_pred.csv"))

## Prior on range of sd (matches that for SPDE2)
meanorig <- 3
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 0.1, ux = 6)$par
mu_logsigma <-paropt[1]
sd_logsigma <- paropt[2]

## Prior on range of Y0 (matches that for SPDE2)
meanorig <- 2000/6371
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 300/6371, ux = 10000/6371)$par
mu_logrho0 <-paropt[1]
sd_logrho0 <- paropt[2]

## Convert to tau/kappa
log_kappa0 <- log_rho2log_kappa(mu_logrho0, 1)
log_tau0 <- log_sigma2log_tau(mu_logsigma, log_kappa0, 1)

## Form process-data map matrices
C0 <- inla.spde.make.A(mesh = mesh, loc = lonlat3D(y_tot$lon, y_tot$lat))
C0_pred <- inla.spde.make.A(mesh = mesh, loc = lonlat3D(y_pred$lon, y_pred$lat))

## Form SPDE with chosen priors
spde <- inla.spde2.matern(mesh = mesh, alpha = 2,
                          B.tau = cbind(log_tau0, -1, 1),
                          B.kappa = cbind(log_kappa0, 0, -1),
                          theta.prior.mean = c(0, 0),
                          theta.prior.prec = c(1/sd_logsigma^2, 1/sd_logrho0^2))

## Run INLA using this SPDE 
system.time({INLA_output <- INLA_spatial_estimation(spde, y_tot$z, C0, C0_pred)})


## Check what hyperparameters were estimated
output.field <- inla.spde2.result(inla = INLA_output,
                                  name = "spatial.field",
                                  spde = spde,
                                  do.transf = TRUE)

## Nominal variance
exp(output.field$summary.log.variance.nominal$mean + 
      output.field$summary.log.variance.nominal$sd^2/2)
## 3.02

## Nominal range
exp(output.field$summary.log.range.nominal$mean + 
      output.field$summary.log.range.nominal$sd^2/2)
## 0.125

## Measurement-error variance (approx.)
1/INLA_output$summary.hyperpar$mean[1]
## 0.16

## Find mesh projection onto validation data locations
project <- inla.mesh.projector(mesh,
                               loc = lonlat3D(y_pred$lon, y_pred$lat))

## Predict the field at these locations
y_pred$spde_pred <-
    inla.mesh.project(project,INLA_output$summary.random$spatial.field$mean +
                              INLA_output$summary.fixed[1,1])

## Prediction standard error at these locations
y_pred$spde_se_obs <-
    inla.mesh.project(project,sqrt(INLA_output$summary.random$spatial.field$sd^2 +
                                   1/INLA_output$summary.hyperpar[1,1]))

## Save results
save(y_pred, file = file.path(cache_folder, "SPDE0_output.rda"))
save(INLA_output, file = file.path(cache_folder, "SPDE0_output_raw.rda"))
