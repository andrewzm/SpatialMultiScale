###############################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Shows the benefit of alternating tilings when Gibbs sampling GMRFs.
##          This code comes with no warranty or guarantee of any kind.
###############################################################################

## Load packages
library("dplyr")
library("Matrix")

## Fix seed
set.seed(1)

## Consider 100 observations
Nobs <- 100

## Function to construct AR1 precision matrix
makeQ <- function(N, phi, sigma2v) {
  i <- c(1, 1, rep(2:(N-1), each = 3L), N, N)
  j <- c(1, 2, c(outer(1:3, 0:(N-3), "+")), N - 1, N)
  x <- c(1, -phi, 
         rep(c(-phi, phi^2 + 1, -phi), N - 2),
         -phi, 1)
  Q <- 1/sigma2v * sparseMatrix(i = i, j = j, x = x)
}

## Function to construct incidence matrix
makeH <- function(s, grid) {
  i <- 1:length(s)
  j <- apply(outer(s, grid, function(x,y) abs(x - y)), 1, which.min)
  sparseMatrix(i = i, j = j, x = 1, dims = c(length(s), length(grid)))
}

## Construct grid
Delta <- 0.01
grid <- seq(Delta/2, 1 - Delta, by = Delta)
n <- length(grid)

## AR1 parameters
sigma2e <- 1
phi <- 0.9
sigma2v <- 0.2
Q <- makeQ(length(grid), phi, sigma2v)
SIGMA <- chol2inv(chol(Q))
L <- t(chol(SIGMA))

## Generate observations and construct incidence matrix
sobs <- runif(Nobs)
Qeps <- 1/sigma2e * Diagonal(n = Nobs)
H <- makeH(sobs, grid)

## Simulate process and data
eta_sim <- (L %*% rnorm(n))
Z <- H %*% eta_sim  + sqrt(sigma2e) * rnorm(Nobs)

## Form tilings ("Groups"), and the Markov blanket of these groups ("Others")
nr <- n + 1e-5
Groups1 <- list(1:(floor(nr / 2)),
               ceiling(nr/2):nr)
Groups2 <- list(1:(floor(nr / 3)),
                ceiling(nr / 3):(floor(2 * nr / 3)),
                ceiling(2 * nr / 3):nr)
Others1 <- list(first(Groups1[[2]]),
                last(Groups1[[1]]))
Others2 <- list(first(Groups2[[2]]),
                c(last(Groups2[[1]]), first(Groups2[[3]])),
                last(Groups2[[2]]))
Groups <- list(Groups1, Groups2)
Others <- list(Others1, Others2)

## Now, for each tiling we construct the incidence matricies, the conditional
## precision matricies, the conditional means, the conditional covariance matrices,
##and the Cholesky factors  of the conditional covariances
Hlist <- Qcondlist <- mucondlist <- Sigmasamplist <- Lsamplist <- list()
for(i in seq_along(Groups)) {
  
  ## For each tiling we have a list, associated with each tile
  Hlist[[i]] <- Qcondlist[[i]] <- mucondlist[[i]] <- 
    Sigmasamplist[[i]] <- Lsamplist[[i]] <- list()
  
  ## For each tile in tiling
  for(j in seq_along(Groups[[i]])) {
    Hlist[[i]][[j]] <- H[, Groups[[i]][[j]]]
    Qcondlist[[i]][[j]] <- Q[Groups[[i]][[j]], Groups[[i]][[j]]]
    mucondlist[[i]][[j]] <- 
      function(i, j, other) -solve(Qcondlist[[i]][[j]]) %*% Q[Groups[[i]][[j]], Others[[i]][[j]]] %*% other
    Q11samp <- Qcondlist[[i]][[j]]
    Sigmasamplist[[i]][[j]] <- solve(Q11samp)
    Lsamplist[[i]][[j]] <- t(chol(Sigmasamplist[[i]][[j]]))
  }
}

## Now we run MCMC. We consider two samplers, and we generate traces for 
## each in eta_samp
N_MCMC <- 10000
eta_samp <- list(matrix(0, n, N_MCMC),
                 matrix(0, n, N_MCMC))

### MCMC scheme (2 possible blocking strategies)
for(grouping in 1:2) {
  for(i in 2:N_MCMC) {
    
    grp <- ifelse(grouping == 2, (i %% 2) + 1, 1)
    eta_current <- eta_samp[[grouping]][, i - 1]
    
    for(j in seq_along(Groups[[grp]])) {
      #mu11samp <- SIGMA11_samp %*% (t(H11) %*% Qeps %*% Z + Q11cond %*% mu11cond(eta_samp[first(grp12), i - 1]))
      musamp <- Sigmasamplist[[grp]][[j]] %*% 
        Qcondlist[[grp]][[j]] %*% 
        mucondlist[[grp]][[j]](grp, j, eta_current[Others[[grp]][[j]]])
      eta_current[Groups[[grp]][[j]]] <- as.numeric(musamp + Lsamplist[[grp]][[j]] %*% 
                                                      rnorm(length(Groups[[grp]][[j]])))
    }
    eta_samp[[grouping]][, i] <- eta_current
  }
}

## Remove first 5000 as burn-in and thin by a factor of two
keep <- seq((i - 5000), i, by = 2)

## Not generate the auto-correlation functions for eta_49 for both samplers
acfseries <- list()
for(grouping in 1:2) {
  MCMCchains <- t(eta_samp[[grouping]][c(1,49), keep])
  par(mfrow = c(2,2)) 
  plot(MCMCchains[,1], type = 'l')
  acf(MCMCchains[,1])
  plot(MCMCchains[,2], type = 'l')
  acf(MCMCchains[,2])
  acfseries[[grouping]] <- sapply(1:n, function(i) acf(eta_samp[[grouping]][i, keep], 
                                                       plot = FALSE)$acf[2])
}

## Plot the results
png("../img/SimACF1.png", width = 1000, height = 1000, res = 300); 
plot(acf(eta_samp[[1]][49, keep]), 
     main = expression('ACF of '* bold(eta)^group("{", 49, "}") * ' from Sampler 1'))
dev.off()

png("../img/SimACF2.png", width = 1000, height = 1000, res = 300); 
plot(acf(eta_samp[[2]][49, keep]), 
     main = expression('ACF of '* bold(eta)^group("{", 49, "}") * ' from Sampler 2'))
dev.off()
