######################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Generates predictions for LatticeKrig
##          This code comes with no warranty or guarantee of any kind.
######################################################################

## Load packages
library("LatticeKrig")
library("spam64")

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load data
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))
y_pred <- read.csv(file.path(cache_folder, "y_pred.csv"))

## Set seed
set.seed(1)

## Try models with three different a.wght parameters
a.wghtgrid <- c(4.01, 5.01, 6.01)

## For each a.wght
for(k in seq_along(a.wghtgrid)) {
  
    ## Fit LatticeKrig model
    proc_time_LTK1 <- system.time({obj <- LatticeKrig(cbind(y_tot$lon,y_tot$lat), 
                                                      y_tot$z, NC = 120,
                                                      nlevel = 3, nu = 1, 
                                                      a.wght = a.wghtgrid[k],
                                                      LKGeometry = "LKRing",
                                                      verbose = TRUE)})
    
    ## Generate predictions and prediction standard errors
    proc_time_LTK2 <- system.time({
                           out.pred <- predict.LKrig(obj, xnew=cbind(y_pred$lon,y_pred$lat))})
    proc_time_LTK3 <- system.time({
                           J <- LKrig.sim.conditional(obj, x.grid = as.matrix(y_pred[,1:2]),M=30)})

    ## Add to the data frame
    y_pred$LK_pred <- as.numeric(out.pred)
    y_pred$LK_se <- apply(J[[3]],1,sd)
    y_pred$LK_se_obs <- sqrt(y_pred$LK_se^2 + obj$sigma.MLE^2)
    
    ## Save results
    save(obj, file = file.path(cache_folder, paste0("LTK_output_raw", k, ".rda")))
    save(proc_time_LTK1, proc_time_LTK2,
         proc_time_LTK3, y_pred, 
         file = file.path(cache_folder, paste0("LTK_output", k, ".rda")))

}