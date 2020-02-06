###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 23 January 2020
## Details: Does Gibbs sampling for the SST data. Requires execution of
##          1_ocean_setup.R.
##          This code comes with no warranty or guarantee of any kind.
###########################################################################

## Load packages
library("dplyr")
library("futile.logger")
library("gtools")
library("Matrix")
library("sparseinv")
library("spam")
library("spam64")
library("parallel")

source("utils.R")

## Set up logger
flog <- flog_console_file()
flog("Running Ocean Sampling via Graph Colouring")
flog("Setting up objects")

## Set local folder where to load data chunks and 
## store intermediates (data and blocks for Gibbs)
cache_folder <- readLines("cache_folder.txt")

## Set seed
set.seed(1L)

## Load the "common" data
load(file = file.path(cache_folder, "data_for_SST2M_no_chunks.rda"))

## Just keep the "id" field of the fine-process data frame, nothing else needed
points_fine <- select(points_fine, id)

## Priors for hyperpars
## ---------------------

## Tune a lognormal prior distribution for sigma (both sigma0 and sigma1)
## such that the lower and upper 0.25 quantiles are 0.1 and 6, respectively
meanorig <- 3
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 0.1, ux = 6)$par
mu_logsigma <- paropt[1]
sd_logsigma <- paropt[2]

## Tune a lognormal prior distribution for rho0 such that the lower 
## sand upper 0.25 quantiles are 300 km and 10000 km, respectively
meanorig <- 2000/6371
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 300/6371, ux = 10000/6371)$par
mu_logrho0 <-paropt[1]
sd_logrho0 <- paropt[2]

## Tune a lognormal prior distribution for rho1 such that the lower 
## sand upper 0.25 quantiles are 30 km and 900 km, respectively
meanorig <- 200/6371
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 30/6371, ux = 900/6371)$par
mu_logrho1 <-paropt[1]
sd_logrho1 <- paropt[2]

## Tune a lognormal prior distribution for sigma_epsilon such that the lower 
## sand upper 0.25 quantiles are 0.5 and 5, respectively
meanorig <- 1
paropt <- optim(par = c(log(meanorig), 0.1), 
                fn = lognorm_cost, lx = 0.5, ux = 5)$par
mu_logsigma_eps <-paropt[1]
sd_logsigma_eps <- paropt[2]

## Initialise parameters and matrices containing samples
## -----------------------------------------------------
N <- 10000         # number of samples
adapt <- 2000      # stop adapting after 2000 samples

## Consolidate the three process scheduling strategies into nchunks
nchunks <- list(nrow(block_wave_map[[1]]),
                nrow(block_wave_map[[2]]),
                nrow(block_wave_map[[3]]))

## Compute the permutation matrix for quick construction of Q0
P0 <- sparseinv:::.amd_Davis(Q0)

## Initialise hyperparameters theta0
hyp_samp0 <- matrix(c(rnorm(1, mean = mu_logsigma, sd = sd_logsigma),
                      rnorm(1, mean = mu_logrho0, sd = sd_logrho0)),
                    2, N)

## Initialise the log_sigma_eps parameters
hyp_samp_log_sigma_eps <- matrix(rnorm(nchunks_kappa, mean = mu_logsigma_eps, sd_logsigma_eps),
                                 N, nchunks_kappa, byrow = TRUE)
hyp_samp_log_sigma_eps_current <- hyp_samp_log_sigma_eps[1, ]

## Initialise the log_rho1 parameters
hyp_samp_log_rho1 <- matrix(rnorm(nchunks_kappa, mean = mu_logrho1, sd_logrho1),
                            N, nchunks_kappa, byrow = TRUE)
hyp_samp_log_rho1_current <- hyp_samp_log_rho1[1, ]

## Initialise the log_sigma1 parameters
hyp_samp_log_sigma1 <- matrix(rnorm(nchunks_kappa, mean = mu_logsigma, sd_logsigma),
                            N, nchunks_kappa, byrow = TRUE)
hyp_samp_log_sigma1_current <- hyp_samp_log_sigma1[1, ]

## Create the interpolated version of the initial parameter fields
load(file.path(cache_folder, "S_matrices.rda"))
points_fine$logsigmainterp <- as.numeric(S %*% hyp_samp_log_sigma1_current)
points_fine$logrhointerp <- as.numeric(S %*% hyp_samp_log_rho1_current)
hyp_samp_log_sigma_eps_current_surface <- as.numeric(S_obs %*% hyp_samp_log_sigma_eps_current)
rm(S, S_obs); gc()

## Randomise the initial sample of Y0
current_samp <- x_samp0_init + rnorm(length(x_samp0_init), mean = 0, sd = 0.2)
xsamp0_mean_all <- 0 * current_samp

## Randomise the initial sample of Y1
current_xsamp1 <- rnorm(nrow(points_fine), mean = 0, sd = 0.1)
xsamp1_mean_all <- xsamp1_sumsq_all <- xsamp1_sum_all <- 0*current_xsamp1

## Initialise proposal distibutions standard deviations
sd_propose_logsigma0 <- sd_propose_logrho0 <- 0.03
sd_propose_logsigma_eps <- sd_propose_logsigma <- sd_propose_logrho1 <-
  rep(0.03, nchunks_kappa)

## contrib 1 is the Y process at the observation locations. Initialise to 0
contrib1 <- 0

## Initialise counter
m = 1

while (m < N) {
  
    m = m + 1
  
    ################################
    ## STEP 1: SAMPLE (theta0, eta0)
    ################################  
    p1 <- proc.time()
    
    ## Update observation error precision matrix
    Qobs <- Diagonal(x = 1/(exp(hyp_samp_log_sigma_eps_current_surface))^2)
    
    ## MH proposal of log(sigma0) and log(rho0)
    propose_lsigma0 <- hyp_samp0[1,(m-1)] + rnorm(1, mean = 0,
                                                  sd = sd_propose_logsigma0)
    propose_lrho0 <- hyp_samp0[2,(m-1)] + rnorm(1, mean = 0,
                                                sd = sd_propose_logrho0)
    propose0 <- c(propose_lsigma0, propose_lrho0)

    ## For old sample and proposed theta0 sample, sample eta0
    LL <- list(hyp_samp0[, (m-1)], propose0)
    eta0_theta0_samp <- mclapply(LL,
                                 logf_marg0,
                                 mu_logsigma = mu_logsigma,
                                 sd_logsigma = sd_logsigma,
                                 mu_logrho = mu_logrho0,
                                 sd_logrho = sd_logrho0,
                                 spde = spde,
                                 P = P0,
                                 z = z_all - contrib1,
                                 Qobs = Qobs,
                                 C_obs = C0,
                                 mc.cores = 2L)
    
    ## Accept/reject based on posterior density evaluation
    logf0_old <- eta0_theta0_samp[[1]]$logf
    logf0_new <- eta0_theta0_samp[[2]]$logf
    alpha <- exp(logf0_new - logf0_old)
    if (runif(1) < alpha) {
        hyp_samp0[, m] <- propose0
        current_samp <- eta0_theta0_samp[[2]]$xsamp
    } else {
        hyp_samp0[, m] <- hyp_samp0[, m-1]
        current_samp <- eta0_theta0_samp[[1]]$xsamp  # Always update x0
    }

    ## Every 30 samples, if still in adaptation phase, update the
    ## proposal distribution standard deviation
    if ((m %% 30) == 0 & m < adapt) {
        sd_propose_logsigma0 <- MHadapt(sd_propose_logsigma0,
                                        hyp_samp0[1, (m-29):m])
        sd_propose_logrho0 <- MHadapt(sd_propose_logrho0,
                                      hyp_samp0[2, (m-29):m])
    }                                         

    ###########################################
    ## STEP 2: SAMPLE (theta1, theta_eps, eta1)
    ########################################### 
    p2 <- proc.time()
    
    ## For each colour
    for (j in 1:4) {
      
      ## Get out the tiles associated with this colour
      blocks_kappa <- block_wave_map_kappa[[j]]
      
      ## For each tile, in parallel, do the following:
      all_hyp_samps <- mclapply(blocks_kappa, function(i) { 
        
        ## Load the data and extract the variables we need
        load(file = file.path(cache_folder, paste0("chunks_kappa", i, ".rda")))
        inidx <- chunk$inidx
        M0 <- chunk$M0
        M1 <- chunk$M1
        M2 <- chunk$M2
        M0_in <- chunk$M0_in
        M1_in <- chunk$M1_in
        M2_in <- chunk$M2_in
        P <- chunk$P
        neighbs <- chunk$neighbs
        xidx <- chunk$affected_states
        S_sub <- chunk$S_sub
        S_own <- chunk$S_own
        zidx_kappa <- chunk$affected_obs_kappa
        zidx_omega <- chunk$affected_obs_omega
        z_kappa <- chunk$z_kappa
        z_omega <- chunk$z_omega
        C_obs_kappa <- chunk$C_obs_kappa
        C_obs_omega <- chunk$C_obs_omega
        
        lsigma_neighb <- hyp_samp_log_sigma1_current[-i]
        lrho_neighb <- hyp_samp_log_rho1_current[-i]
        
        ## Establish what is the old (sigma1, rho1) and the old surfaces
        old_lsigma <- hyp_samp_log_sigma1_current[i]
        old_lrho <- hyp_samp_log_rho1_current[i]
        old_lsigma_surface <- S_sub %*% lsigma_neighb + S_own*old_lsigma
        old_lrho_surface <- S_sub %*% lrho_neighb + S_own*old_lrho
        old <- c(old_lsigma, old_lrho)
        
        ## Propose a new sigma1 and rho1 and update the surfaces
        propose_lsigma <- old_lsigma + rnorm(1, mean = 0,
                                             sd = sd_propose_logsigma[i])
        propose_lrho <- old_lrho + rnorm(1, mean = 0,
                                         sd = sd_propose_logrho1[i])
        propose_lsigma_surface <- S_sub %*% lsigma_neighb + S_own*propose_lsigma
        propose_lrho_surface <- S_sub %*% lrho_neighb + S_own*propose_lrho
        propose <- c(propose_lsigma, propose_lrho)
        this_xsamp <- current_xsamp1[xidx]
        
        ## Compute log density under the old parameters 
        eta1_theta1_samp_old <- logf_marg1(theta = old,
                                           theta_surface = cbind(old_lsigma_surface,
                                                                 old_lrho_surface),
                                           mu_logsigma = mu_logsigma,
                                           sd_logsigma = sd_logsigma,
                                           mu_logrho = mu_logrho1,
                                           sd_logrho = sd_logrho1,
                                           M0 = M0, M1 = M1, M2 = M2,
                                           M0_in = M0_in, M1_in = M1_in, M2_in = M2_in,
                                           Yx = this_xsamp,
                                           P = P,
                                           inidx = inidx,
                                           z = z_kappa - C0[zidx_kappa, ] %*% current_samp,
                                           Qobs = Diagonal(x = 1/(exp(hyp_samp_log_sigma_eps_current_surface[zidx_kappa]))^2),
                                           C_obs = C_obs_kappa)
          
        ## Compute log density under the new proposed parameters
        eta1_theta1_samp_new <- logf_marg1(theta = propose,
                                           theta_surface = cbind(propose_lsigma_surface,
                                                                 propose_lrho_surface),
                                           mu_logsigma = mu_logsigma,
                                           sd_logsigma = sd_logsigma,
                                           mu_logrho = mu_logrho1,
                                           sd_logrho = sd_logrho1,
                                           M0 = M0, M1 = M1, M2 = M2,
                                           M0_in = M0_in, M1_in = M1_in, M2_in = M2_in,
                                           Yx = this_xsamp,
                                           P = P,
                                           inidx = inidx,
                                           z = z_kappa - C0[zidx_kappa, ] %*% current_samp,
                                           Qobs = Diagonal(x = 1/(exp(hyp_samp_log_sigma_eps_current_surface[zidx_kappa]))^2),
                                           C_obs = C_obs_kappa)

        ## Accept/reject eta1 and theta jointly
        alpha <- exp(eta1_theta1_samp_new$logf - eta1_theta1_samp_old$logf)
        if (runif(1) < alpha) {
            hypsamp1 <- propose
            xsamp1_part <- eta1_theta1_samp_new$xsamp
            this_xsamp[inidx] <- xsamp1_part
        } else {
            hypsamp1 <- old
            xsamp1_part <- eta1_theta1_samp_old$xsamp
            this_xsamp[inidx] <- xsamp1_part                
        }

        ## Update measurement-error standard deviation
        ## First extract quantities we need from the chunk
        S_obs_sub <- chunk$S_obs_sub
        S_obs_own <- chunk$S_obs_own
        lsigma_eps_neighb <- hyp_samp_log_sigma_eps_current[-i]
        
        ## Old log_sigma_eps and the old surface
        old_lsigma_eps <- hyp_samp_log_sigma_eps_current[i]
        old_lsigma_eps_surface <- S_obs_sub %*% lsigma_eps_neighb + S_obs_own*old_lsigma_eps
          
        ## Propose new log_sgma_eps and generate new surface
        propose_lsigma_eps <- old_lsigma_eps + rnorm(1, mean = 0,
                                                     sd = sd_propose_logsigma_eps[i])
        propose_lsigma_eps_surface <- S_obs_sub %*% lsigma_eps_neighb + S_obs_own*propose_lsigma_eps
        
        ## Use new eta generated above when updating; generate a corresponding
        ## data/eta mapping matrix
        this_C <- cbind(C0[zidx_omega, ], C_obs_omega)
        this_xsamp_all <- c(current_samp, this_xsamp)
          
        ## Compute log density at old log_sigma_eps
        logf0 <- logf_sigma_eps(theta = old_lsigma_eps,
                                theta_surface = as.numeric(old_lsigma_eps_surface),
                                mu_lsigma_eps = mu_logsigma_eps,
                                sd_lsigma_eps = sd_logsigma_eps,
                                C_obs = this_C,
                                Yx = this_xsamp_all,
                                z = z_omega)
        
        ## Compute log density at new log_sigma_eps  
        logf1 <- logf_sigma_eps(theta = propose_lsigma_eps,
                                  theta_surface = as.numeric(propose_lsigma_eps_surface),
                                  mu_lsigma_eps = mu_logsigma_eps,
                                  sd_lsigma_eps = sd_logsigma_eps,
                                  C_obs = this_C,
                                  Yx = this_xsamp_all,
                                  z = z_omega)
          
        ## Accept/Reject
        alpha <- exp(logf1 - logf0)
        if (runif(1) < alpha) {
          lsigma_eps <- propose_lsigma_eps
        } else {
          lsigma_eps <- old_lsigma_eps
        }
        
        ## Return the sampled eta and parameters
        list(hypsamp1 = hypsamp1,
             lsigma_eps = lsigma_eps,
             xsamp1_part = xsamp1_part,
             xidx = xidx,
             inidx = inidx)
      }, mc.cores = min(length(blocks_kappa), detectCores() - 2L))  
      
      ## Once a colour is completed, we need to update the current chain
      ## We do this for each block in sequence
      for (i in seq_along(blocks_kappa)) {
        block <- blocks_kappa[i]
        hyp_samp_log_sigma1_current[block] <- all_hyp_samps[[i]]$hypsamp1[1]
        hyp_samp_log_rho1_current[block] <- all_hyp_samps[[i]]$hypsamp1[2]
        hyp_samp_log_sigma_eps_current[block] <- all_hyp_samps[[i]]$lsigma_eps
        
        hyp_samp_log_sigma1[m, block] <- all_hyp_samps[[i]]$hypsamp1[1]
        hyp_samp_log_rho1[m, block] <- all_hyp_samps[[i]]$hypsamp1[2]
        hyp_samp_log_sigma_eps[m, block] <- all_hyp_samps[[i]]$lsigma_eps

        ## Update eta samples 
        ## There is no need to update contrib1 since we don't use it for Step 3;
        ## we update it after sampling all of eta1
        inidx <- all_hyp_samps[[i]]$inidx
        xidx <- all_hyp_samps[[i]]$xidx
        current_xsamp1[xidx][inidx] <- all_hyp_samps[[i]]$xsamp1_part     
        
        ## Every 30 samples, if still in adaptation phase, update the
        ## proposal distributions standard deviations
        if ((m %% 30) == 0 & m < adapt) {
          sd_propose_logrho1[block] <- MHadapt(sd_propose_logrho1[block],
                                               hyp_samp_log_rho1[(m-29):m, block])
          sd_propose_logsigma[block] <- MHadapt(sd_propose_logsigma[block],
                                                hyp_samp_log_sigma1[(m-29):m, block])
          sd_propose_logsigma_eps[block] <- MHadapt(sd_propose_logsigma_eps[block],
                                                    hyp_samp_log_sigma_eps[(m-29):m, block])
        }                                         
      }
    }
    
    ## Interpolate to obtain the updated parameter surfaces
    load(file.path(cache_folder, "S_matrices.rda"))
    points_fine$logsigmainterp <- as.numeric(S %*% hyp_samp_log_sigma1_current)
    points_fine$logrhointerp <- as.numeric(S %*% hyp_samp_log_rho1_current)
    hyp_samp_log_sigma_eps_current_surface <- as.numeric(S_obs %*% hyp_samp_log_sigma_eps_current)
    rm(S, S_obs); gc()

    ############################
    ## STEP 3: SAMPLE eta1 again
    ############################
    p3 <- proc.time()

    ## We rotate the tilings (colourings) every m, so that no states are 
    ## near the boundaries for each m
    colouring <- (m %% 3) + 1
    
    ## Extract the number of waves (colours) we need to go through for this
    ## tiling
    wave_nums <- sort(unique(block_wave_map[[colouring]]$wave))
    
    ## For each colour
    for (j in wave_nums) {
      
      ## Extract the tiles
      block_nums <- subset(block_wave_map[[colouring]], wave == j)$block
      
      ## In parallel, for each tile, do the following:
      all_samps <- mclapply(block_nums, function(i) { 
        
          ## Load the required chunk and extract the objects we need
          load(file = file.path(cache_folder, paste0("chunks",i,"_colouring", colouring, ".rda")))
          C0_in <- chunk$C0_in
          C0_out <- chunk$C0_out
          in_idx0 <- chunk$in_idx0
          out_idx0 <- chunk$out_idx0
          obs_idx <- chunk$obs_idx
          Qobs_in <- Diagonal(x = 1/(exp(hyp_samp_log_sigma_eps_current_surface[obs_idx]))^2)
          C1_in <- chunk$C1_in
          C1_out <- chunk$C1_out
          M0 <- chunk$M0
          M1 <- chunk$M1
          M2 <- chunk$M2
          Mall_df1 <- chunk$Mall_df
          QC1 <- chunk$QC
          ext_idx1 <- chunk$ext_idx1
          out_idx1 <- chunk$out_idx1
          all_idx1 <- chunk$all_idx1
          in_idx1 <- chunk$in_idx1
          
          ## Compute the fine-scale contributions from neighbouring vertices
          neigh_contrib1 <- C1_out %*% current_xsamp1[out_idx1]
          
          ## Compute kappa1 and tau1 from rho1 and sigma1 (the surfaces)
          this_log_kappa1 <- log_rho2log_kappa(points_fine$logrhointerp[all_idx1], 1)
          this_log_tau1 <- log_sigma2log_tau(points_fine$logsigmainterp[all_idx1], this_log_kappa1, 1)
          
          ## Compute the precision matrix associated with these values of
          ## rho and sigma
          Q1_all <- inla.spde.precision.fast2(
              M0, M1, M2,
              theta = cbind(this_log_tau1,
                            this_log_kappa1))
          Q1_in <- Q1_all[-ext_idx1, -ext_idx1]
          Q1_out <- Q1_all[-ext_idx1, ext_idx1]
       
          ## Update eta1 in this tile conditional on Y0 and neighbouring vertices
          contrib0 <- C0_in %*% current_samp[in_idx0] +
              C0_out %*% current_samp[out_idx0]
          Qtildejj1 <- t(C1_in) %*% Qobs_in %*% C1_in + Q1_in
          ybar1 <- t(C1_in) %*% Qobs_in %*% (chunk$z - neigh_contrib1 - contrib0) -
              Q1_out %*% current_xsamp1[chunk$out_idx1]
          Qtildejj1_spam <- as.spam.dgCMatrix(Qtildejj1)
          Ppost <- sparseinv:::.amd_Davis(Qtildejj1)
          R1 <- spam::chol(Qtildejj1_spam, pivot = Ppost)
          this_mean <- backsolve(R1, forwardsolve(R1, ybar1, transpose = T, upper.tri= T))
          xsamp1 <- as.vector(this_mean + backsolve(R1, rnorm(length(this_mean))))
          contrib1 <- C1_in %*% xsamp1
          
          ## Same output but slower:
          ## X <- cholPermute(Qtildejj1)
          ## this_mean <- cholsolve(Qtildejj1, ybar1, cholQp = X$Qpermchol, P = X$P)
          ## this_samp <- this_mean + X$P %*% solve(t(X$Qpermchol),z)
      
          ## Collect into a list and return
          list(i = i,
               xsamp1 = xsamp1,
               contrib1 = contrib1,
               obs_idx = chunk$obs_idx)
    }, mc.cores = length(block_nums))
    
    
    ## Wave completed, update the current state of information
    for (i in 1:length(all_samps)) {
      current_xsamp1[in_idx_lists[[colouring]][[block_nums[i]]]] <- all_samps[[i]]$xsamp1
      contrib1[all_samps[[i]]$obs_idx] <- all_samps[[i]]$contrib1
    }
  }
  
  ## Keep a running mean, sum of squares, and variance of the traces
  xsamp0_mean_all <- (current_samp + (m-1)*xsamp0_mean_all)/m          
  xsamp1_mean_all <- (current_xsamp1 + (m-1)*xsamp1_mean_all)/m          
  xsamp1_sumsq_all <- xsamp1_sumsq_all + current_xsamp1^2
  xsamp1_sum_all <- xsamp1_sum_all + current_xsamp1
  xsamp1_var_all <- (xsamp1_sumsq_all - (xsamp1_sum_all)^2/m)/(m - 1)
  p4 <- proc.time()
  
  ## Save sample to disk
  save(current_samp, file = file.path(cache_folder, paste0("etasamp0_", m)))
  save(current_xsamp1, file = file.path(cache_folder, paste0("etasamp1_", m)))


  ## Every 20th sample also save the hyperparameters
  if ((m %% 20) == 0) {
    save(hyp_samp0, hyp_samp_log_rho1, hyp_samp_log_sigma1, hyp_samp_log_sigma_eps, m,
         file = file.path(cache_folder, paste0("hypsamps_MCMC.rda")))
  }
  
  ## Log progress
  flog(paste0("Sample ", m, " completed. Time for sampling (theta0, eta0): ",
              formatC((p2 - p1)[3], digits = 3)," s. ",
              "Time for sampling (theta1, thetaeps, eta1): ",
              formatC((p3 - p2)[3],digits = 3)," s. ",
              "Time for sampling eta1: ",
              formatC((p4 - p3)[3], digits = 3)," s. "))
  flog(paste0("...Coarse scale hyperparameters: ", hyp_samp0[1, m],"  ", hyp_samp0[2, m]))
  flog(paste0("...Range coarse states: ", range(current_samp)))
  flog(paste0("...Range fine states: ", range(current_xsamp1)))
  flog(paste0("----------------------------------------------------"))
}
