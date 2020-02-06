######################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Generates predictions for NNGP
##          This code comes with no warranty or guarantee of any kind.
######################################################################

## Load packages
library("dplyr")
library("spNNGP")

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load data
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))
y_tot <- read.csv(file.path(cache_folder, "y_pred.csv"))

## Fix seed
set.seed(1)

## Generate grid for scale (phi) and inverse SNR (alpha)
theta.alpha = expand.grid(phi = seq(0.01, 1, by = 0.1),
                          alpha = seq(0.001, 0.1, by = 0.01),
                          nu = 1) %>% as.matrix()

## Fit conjugage NNGPs on this grid and choose the best one
proc_time_NNGP <- system.time({
    ConjNNGP <- spConjNNGP(z ~ 1, data = y_tot, 
                           coords = as.matrix(y_tot[c("lon", "lat")]), 
                           theta.alpha = theta.alpha,
                           sigma.sq.IG = c(0.001, 0.001),
                           X.0 = matrix(1, nrow(y_pred), 1),
                           coords.0 = as.matrix(y_pred[c("lon", "lat")]),
                           n.neighbors = 15,
                           n.omp.threads = 30,
                           cov.model = "matern")})

## Add to data frame
y_pred$NNGP_pred <- as.numeric(ConjNNGP$y.0.hat)
y_pred$NNGP_se_obs <- as.numeric(sqrt(ConjNNGP$y.0.hat.var))

## Save results
save(ConjNNGP, file = file.path(cache_folder, "NNGP_output_raw.rda"))
save(y_pred, proc_time_NNGP, file = file.path(cache_folder, "NNGP_output.rda"))

## Matern scale parameter
phi <- ConjNNGP$theta.alpha[1]

## Nominal range with nu = 1
sqrt(8)/phi
## 6.89

## Estimated process variance
ConjNNGP$sigma.sq.hat
## 8.25

## alpha = tau^2/sigma^2. We can use this to find the estimated 
## measurement-error variance
alpha <- ConjNNGP$theta.alpha[2]
alpha * sigma2hat
## 0.091
