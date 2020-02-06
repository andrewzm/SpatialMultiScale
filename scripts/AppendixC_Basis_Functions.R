###############################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Shows the effect of grid spacing on predictions in a two-scale
##          model.
##          This code comes with no warranty or guarantee of any kind.
###############################################################################

## Load packages
library("dplyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("Matrix")
library("sparseinv")
library("verification")

## Fix seed
set.seed(1L)

## Model both Y0 and Y1 (2 scale)
modelY0 <- TRUE
modelY1 <- TRUE

## Simulate 1100 data points. We shall use 500 outside the gap for training,
## 500 outside the gap for validation, and 100 inside the gap for validation.
Ntrain <- 500
Nval1 <- 500
Nval2 <- 100
N <- Ntrain + Nval1 + Nval2

## Function to make the precision matrix for an AR1 process
makeQ <- function(N, phi, sigma2v) {
  i <- c(1, 1, rep(2:(N-1), each = 3L), N, N)
  j <- c(1, 2, c(outer(1:3, 0:(N-3), "+")), N - 1, N)
  x <- c(1, -phi, 
         rep(c(-phi, phi^2 + 1, -phi), N - 2),
         -phi, 1)
  Q <- 1/sigma2v * sparseMatrix(i = i, j = j, x = x)
}

## Function to make the incidence matrix mapping grid to observations
makeH <- function(s, grid) {
  i <- 1:length(s)
  j <- apply(outer(s, grid, function(x,y) abs(x - y)), 1, which.min)
  sparseMatrix(i = i, j = j, x = 1, dims = c(length(s), length(grid)))
}

## We will consider several values of Delta
Delta_0 <- seq(0.01, 0.1, by = 0.02)
Delta_1 <- seq(0.001, 0.01, by = 0.002)
Diags <- expand.grid(Delta_0 = Delta_0,
                     Delta_1 = Delta_1) %>%
  mutate(RMSPE1 = 0, CRPS1 = 0,
         RMSPE2 = 0, CRPS2 = 0)

## We will get an RMSPE and CRPS based on 100 simulations
Nsims <- 100

## Initialise the RMSPE and CRPS result matrices
RMSPE1 <- RMSPE2 <- CRPS1 <- CRPS2 <- matrix(0, nrow(Diags), Nsims)
RMSPE_opt1 <- RMSPE_opt2 <- CRPS_opt1 <- CRPS_opt2 <- rep(0, Nsims)

## Fix the process parameters as follows
sigma2e <- 0.0002
tau_0 <- 0.4
sigma2_0 <- 1
tau_1 <- 0.04
sigma2_1 <- 0.05

## Repeat the following Nsims times
for(Sim in 1:Nsims) {
  cat(paste0("Doing Sim ", Sim, "\n"))
  
  ## Randomise a "gap" in our domain of width 0.2
  minbox <- runif(n = 1, 0, 0.8)
  
  ## Put 50% of training data to left of gap, 50% to the right
  strain <- c(runif(Ntrain/2, max = minbox),
              runif(Ntrain/2, min = minbox + 0.2))
  ## Same with validation data 1
  sval1 <- c(runif(Nval1/2, max = minbox),
             runif(Nval1/2, min = minbox + 0.2))
  
  ## These validation data are in the gap
  sval2 <- runif(Nval2, min = minbox, max = minbox + 0.2)
  
  ## Collect ALL data in sobs
  sobs <-c(strain, sval1, sval2)
  
  ## Simulate data for these locations
  Dobs <- as.matrix(dist(sobs))
  Cobs <- sigma2_0*exp(-Dobs/tau_0) + sigma2_1*exp(-Dobs/tau_1)  
  L <- t(chol(Cobs))
  Z <- L %*% rnorm(N) + sqrt(sigma2e) * rnorm(N)

  ## Generate the indices for our data (training and validation)
  idxtrain <- 1:Ntrain
  idxval1 <- (Ntrain + 1) : (Ntrain + Nval1)
  idxval2 <- -(1:(Ntrain + Nval1))
  
  ## Construct the training and validation data sets
  Ztrain <- Z[idxtrain]
  Zval1 <- Z[idxval1]
  Zval2 <- Z[idxval2]

  ## Find the diagnostics using the true covariance function
  Dtrain <- fields::rdist(strain, strain)
  Dpred_train <- fields::rdist(c(sval1, sval2), strain)
  Dpred <- fields::rdist(c(sval1, sval2), c(sval1, sval2))

  Ctrain <- sigma2_0*exp(-Dtrain/tau_0) + sigma2_1*exp(-Dtrain/tau_1)  
  Cpred_train <- sigma2_0*exp(-Dpred_train/tau_0) + sigma2_1*exp(-Dpred_train/tau_1)  
  Cpred <- sigma2_0*exp(-Dpred/tau_0) + sigma2_1*exp(-Dpred/tau_1)
  
  Ctrain_inv <- chol2inv(chol(Ctrain))
  pred_opt <- Cpred_train %*% Ctrain_inv %*% Ztrain
  Sigma_opt <- Cpred - Cpred_train %*% Ctrain_inv %*% t(Cpred_train)
  var_opt <- diag(Sigma_opt)
  
  RMSPE_opt1[Sim] <- sqrt(mean((pred_opt[1:Nval1] - Zval1)^2))
  RMSPE_opt2[Sim] <- sqrt(mean((pred_opt[-(1:Nval1)] - Zval2)^2))
  CRPS_opt1[Sim] <- crps(Zval1, cbind(pred_opt[1:Nval1], sqrt(var_opt)[1:Nval1]))$CRPS
  CRPS_opt2[Sim] <- crps(Zval2, cbind(pred_opt[-(1:Nval1)], sqrt(var_opt)[-(1:Nval1)]))$CRPS
  
  ## Now we compute the diagnostics using GMRF approximations to the 
  ## exponential covariance functions, for various grid-spacings
  for(i in 1:nrow(Diags)) {
    
    ## These spacings
    Delta_0 <- Diags$Delta_0[i]
    Delta_1 <- Diags$Delta_1[i]
    
    ## Construct grid according to these spacings
    grid0 <- seq(Delta_0/2, 1 - Delta_0/2, by = Delta_0)
    grid1 <- seq(Delta_1/2, 1 - Delta_1/2, by = Delta_1)
    
    ## Calculate the equivalent AR1 parameters and variance parameters
    phi_0 <- exp(-Delta_0/tau_0)
    phi_1 <- exp(-Delta_1/tau_1)
    sigma2v_0 <- sigma2_0*(1 - phi_0^2)
    sigma2v_1 <- sigma2_1*(1 - phi_1^2)
    
    ## Construct the precision and incidence matrices on this grid for 
    ## both processes
    Q0 <- makeQ(length(grid0), phi_0, sigma2v_0)
    Q1 <- makeQ(length(grid1), phi_1, sigma2v_1)
    H0train <- makeH(strain, grid0)
    H1train <- makeH(strain, grid1)
    H0val1 <- makeH(sval1, grid0)
    H1val1 <- makeH(sval1, grid1)
    H0val2 <- makeH(sval2, grid0)
    H1val2 <- makeH(sval2, grid1)
    
    ## We can do just Y0 or just Y1 (we do both here, so we concatenate
    ## the matrices appropriately)
    if(!modelY1) {
      Q <- Q0
      Htrain <- H0train
      Hval1 <- H0val1
      Hval2 <- H0val2
    } else if (!modelY0) {
      Q <- Q1
      Htrain <- H1train
      Hval1 <- H1val1
      Hval2 <- H1val2
    } else  {
      Q <- bdiag(Q0, Q1)
      Htrain <- cbind(H0train, H1train)
      Hval1 <- cbind(H0val1, H1val1)
      Hval2 <- cbind(H0val2, H1val2)
    }
    
    ## Measurement-error precision
    Qeps <- 1/sigma2e * Diagonal(n = Ntrain)
    
    ## Conditional precision and mean
    Qstar <- t(Htrain) %*% Qeps %*% Htrain + Q
    mustar <- solve(Qstar, t(Htrain) %*% Qeps %*% Ztrain)

    ## Conditional covariance matrix
    Sigmastar <- chol2inv(chol(Qstar))
    
    ## Validation 1 diagnostics
    predval1 <- Hval1 %*% mustar
    varval1 <- rowSums((Hval1 %*% Sigmastar) * Hval1)
    RMSPE1[i, Sim] <- sqrt(mean((predval1 - Zval1)^2))
    CRPS1[i, Sim] <- crps(Zval1, cbind(predval1, sqrt(varval1)))$CRPS
    
    ## Validation 2 diagnostics
    predval2 <- Hval2 %*% mustar
    varval2 <- rowSums((Hval2 %*% Sigmastar) * Hval2)
    RMSPE2[i, Sim] <- sqrt(mean((predval2 - Zval2)^2))
    CRPS2[i, Sim] <- crps(Zval2, cbind(predval2, sqrt(varval2)))$CRPS
  }
}

## Summarise the diagnostics by finding the mean
Diags$RMSPE1 <- apply(RMSPE1, 1, mean)
Diags$RMSPE2 <- apply(RMSPE2, 1, mean)
Diags$CRPS1 <- apply(CRPS1, 1, mean)
Diags$CRPS2 <- apply(CRPS2, 1, mean)

## Plot the results
g1 <- ggplot(Diags) + geom_tile(aes(Delta_0, Delta_1, fill = RMSPE1)) + 
    theme_bw() + scale_fill_distiller(palette = "Greys") +
    xlab(expression(delta[0])) + ylab(expression(delta[1]))
g2 <- ggplot(Diags) + geom_tile(aes(Delta_0, Delta_1, fill = RMSPE2)) + 
  theme_bw() + scale_fill_distiller(palette = "Greys") +
    xlab(expression(delta[0])) + ylab(expression(delta[1]))
g3 <- ggplot(Diags) + geom_tile(aes(Delta_0, Delta_1, fill = CRPS1)) + 
  theme_bw() + scale_fill_distiller(palette = "Greys") +
    xlab(expression(delta[0])) + ylab(expression(delta[1]))
g4 <- ggplot(Diags) + geom_tile(aes(Delta_0, Delta_1, fill = CRPS2)) + 
  theme_bw() + scale_fill_distiller(palette = "Greys") +
    xlab(expression(delta[0])) + ylab(expression(delta[1]))

gall <- grid.arrange(g1, g2, g3, g4, nrow = 2)
ggsave(gall, file = "../img/RMSPE_CRPS.png", width = 8, height = 6)
