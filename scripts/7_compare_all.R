###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Compares the results from all the models on the SST data.
##          This code comes with no warranty or guarantee of any kind.
###########################################################################

## Load packages
library("dplyr")
library("fields")
library("ggplot2")
library("parallel")
library("tidyr")
library("verification")
library("xtable")
source("utils.R")

## Define the functions we will use for comparing across methods

## 90% coverage
coverage90 <- function(z, mu, se) {
    lower <- mu - 1.64*se
    upper <- mu + 1.64*se
    sum((z < upper) & (z > lower)) / length(z)
}

## 90% interval score
IS90 <- function(true, mu, se) {
  alpha = 0.1
  pred90l <- mu - 1.64*se
  pred90u <- mu + 1.64*se
  ISs <- (pred90u - pred90l) + 2/alpha * (pred90l - true) * (true < pred90l) +
                        2/alpha * (true - pred90u) * (true > pred90u)
  mean(ISs)
}

## Root-mean-squared prediction error
RMSPE <- function(z,pred) {
    Y <- (z - pred)^2
    sqrt(mean(Y))
}

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load training data
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))

## SPDE1-indep
load(file = file.path(cache_folder, "SPDE1-indep_output.rda"))
y_pred_all <- y_pred %>% dplyr::select(lon, lat, sst, z, spde_hi_pred, spde_hi_se_obs)

## SPDE0
load(file = file.path(cache_folder, "SPDE0_output.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, spde_pred, spde_se_obs))

## LTK a = 6.01
load(file = file.path(cache_folder, "LTK_output3.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, LK6_pred = LK_pred, LK6_se_obs = LK_se_obs))

## LTK a = 5.01
load(file = file.path(cache_folder, "LTK_output2.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, LK5_pred = LK_pred, LK5_se_obs = LK_se_obs))

## LTK a = 4.01
load(file = file.path(cache_folder, "LTK_output1.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, LK4_pred = LK_pred, LK4_se_obs = LK_se_obs))

## NNGP
load(file = file.path(cache_folder, "NNGP_output.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, NNGP_pred, NNGP_se_obs))

## FSA
FSA_pred <- read.csv(file.path(cache_folder, "FSA_results.csv"), header = FALSE)
names(FSA_pred) <- c("lon", "lat", "FSA_pred", "FSA_se")
y_pred_all <- left_join(y_pred_all, FSA_pred)
## sigma2:5.2529,phi:3.0625,nu: 1.0000, nugget:0.1058,MSPE:0.12584

## Two-scale GMRF
load(file = file.path(cache_folder, "SPDE2_output.rda"))
y_pred_all <- cbind(y_pred_all, dplyr::select(y_pred, TwoGMRF_pred, TwoGMRF_se_obs))

## There are some few (5 in all) NAs from SPDE1-indep, remove these
rmidx <- which(is.na(y_pred_all$spde_hi_pred))
y_pred_all <- y_pred_all[-rmidx, ]

## Now we find the indices for Z_v^(2)
## Note: This takes a long time to compute
if(!file.exists(file.path(cache_folder, "idx_0_75_box_val.rda"))) {
    
  ## Do this in batches of 1000
  batches <- cut(1 : nrow(y_pred_all), 
                 seq(0, nrow(y_pred_all) + 1000, by = 1000), 
                 labels = FALSE)
  
  idx_logical <- rep(NA, nrow(y_pred_all))
  for(j in unique(batches)) {
    this_idx <- which(batches == j)
  
    idx_logical[this_idx] <- (mclapply(this_idx,#nrow(y_pred_all),
                               function(i) {
                                 thisy <- y_pred_all[i, ]
                                 ysub <- filter(y_tot, lon < thisy$lon + 0.75 & lon > thisy$lon - 0.75 &
                                                  lat < thisy$lat + 0.75 & lat > thisy$lat - 0.75)
                                 
                                 ## If there are no data points then add this to Z_v^(2)
                                 ifelse(nrow(ysub) == 0, TRUE, FALSE)
                               }, mc.cores = 10L) %>% unlist())
    cat(paste0("Finished Batch ", j, " out of ", length(unique(batches)), "\n"))
  }
  
  idx <- which(idx_logical)

  ## Save the indices of the data points for Z_v^(2)
  save(idx, file = file.path(cache_folder, "idx_0_75_box_val.rda"))
} else {
  load(file = file.path(cache_folder, "idx_0_75_box_val.rda"))
}

## Now find those validation data inside the 8 x 8 box in the Pacific Ocean for Z_v^(3)
idx_box <- which(y_pred_all$lat > -4 & y_pred_all$lat < 4 & 
                y_pred_all$lon > -159 & y_pred_all$lon < -151)

## Remove indices from Z_v^(2) that are in this box
idx_outside_box_sparse <- setdiff(idx, idx_box)

## Remove indices from Z_v^(1) that are in the box and in Z_v^(2)
idx_outside_box_dense <- setdiff(1:nrow(y_pred_all), 
                                 c(idx_box, idx_outside_box_sparse))

## Construct Z_v^(1)
y_pred_dense <- y_pred_all[idx_outside_box_dense, ]

## Construct Z_v^(2)
y_pred_sparse <- y_pred_all[idx_outside_box_sparse, ]

## Construct Z_v^(3)
y_pred_box <- y_pred_all[idx_box, ]

## Do some distances check. How far away are Z_v^(2) from training data?
dists2 <- fields::rdist(cbind(y_pred_all$lon[idx_outside_box_sparse],
                              y_pred_all$lat[idx_outside_box_sparse]),
                        cbind(y_tot$lon, y_tot$lat))
min_dists2 <- apply(dists2, 1, min)

## How far away are Z_v^(3) from training data?
## Do in batches of 1000
batches <- cut(1 : length(idx_box), 
               seq(0, length(idx_box) + 1000, by = 1000), 
               labels = FALSE)
min_dists3 <- c()
for (i in unique(batches)) {
  dists3 <- fields::rdist(cbind(y_pred_all$lon[idx_box[batches == i]], 
                                y_pred_all$lat[idx_box[batches == i]]),
                          cbind(y_tot$lon, y_tot$lat))
  min_dists3 <- c(min_dists3, apply(dists3, 1, min))
  cat(paste0("Finished Batch ", i, " out of ", max(batches), "\n"))
}

quantile(min_dists2, c(0,0.95))
quantile(min_dists3, c(0,0.6))

## Create a function that generates the results
generate_results <- function(df) {
  data.frame(Model = c("SPDE2", "FSA",  "SPDE1-indep", "LTK4", "LTK5", 
                       "LTK", "NNGP", "SPDE0"),
             RMSPE = c(RMSPE(df$z, df$TwoGMRF_pred),
                      RMSPE(df$z, df$FSA_pred),
                      RMSPE(df$z, df$spde_hi_pred),
                      RMSPE(df$z, df$LK4_pred),
                      RMSPE(df$z, df$LK5_pred),
                      RMSPE(df$z, df$LK6_pred),
                      RMSPE(df$z, df$NNGP_pred),
                      RMSPE(df$z, df$spde_pred)),
             CRPS = c(crps(df$z, cbind(df$TwoGMRF_pred, df$TwoGMRF_se_obs))$CRPS,
                      crps(df$z, cbind(df$FSA_pred, df$FSA_se))$CRPS,
                      crps(df$z, cbind(df$spde_hi_pred, df$spde_hi_se_obs))$CRPS,
                      crps(df$z, cbind(df$LK4_pred, df$LK4_se_obs))$CRPS,
                      crps(df$z, cbind(df$LK5_pred,df$LK5_se_obs))$CRPS,
                      crps(df$z, cbind(df$LK6_pred ,df$LK6_se_obs))$CRPS,                      
                      crps(df$z, cbind(df$NNGP_pred, df$NNGP_se_obs))$CRPS,
                      crps(df$z, cbind(df$spde_pred, df$spde_se_obs))$CRPS),
             IS90 = c(IS90(df$z,df$TwoGMRF_pred,df$TwoGMRF_se_obs),
                      IS90(df$z,df$FSA_pred,df$FSA_se),
                      IS90(df$z,df$spde_hi_pred,df$spde_hi_se_obs),
                      IS90(df$z,df$LK4_pred,df$LK4_se_obs),
                      IS90(df$z,df$LK5_pred,df$LK5_se_obs),
                      IS90(df$z,df$LK6_pred,df$LK6_se_obs),                      
                      IS90(df$z,df$NNGP_pred,df$NNGP_se_obs),
                      IS90(df$z,df$spde_pred,df$spde_se_obs)),
             Cov90 = c(coverage90(df$z, df$TwoGMRF_pred, df$TwoGMRF_se_obs),
                       coverage90(df$z, df$FSA_pred, df$FSA_se),
                       coverage90(df$z, df$spde_hi_pred, df$spde_hi_se_obs),
                       coverage90(df$z, df$LK4_pred, df$LK4_se_obs),
                       coverage90(df$z, df$LK5_pred, df$LK5_se_obs),
                       coverage90(df$z, df$LK6_pred, df$LK6_se_obs),                       
                       coverage90(df$z, df$NNGP_pred, df$NNGP_se_obs),
                       coverage90(df$z, df$spde_pred, df$spde_se_obs)))
}

## Generate the results for Z_v^(1), Z_v^(2), and Z_v^(3)
df_results_dense <- generate_results(y_pred_dense)
df_results_sparse <- generate_results(y_pred_sparse)
df_results_box <- generate_results(y_pred_box)

## Print for latex table
print(xtable(df_results_dense %>% arrange(Model), 
             digits = 3, align = c("llcccc")), 
      include.rownames = FALSE)
print(xtable(df_results_sparse %>% arrange(Model), 
             digits = 3, align = c("llcccc")), 
      include.rownames = FALSE)
print(xtable(df_results_box %>% arrange(Model), 
             digits = 3, align = c("llcccc")), 
      include.rownames = FALSE)

##############################################################
## Generate Fig.7 of the paper, visualising predictions in box
##############################################################

## First, extract data in a region enclosing the box of interest
idx_bigger_box <- which(y_pred_all$lat > -6 & y_pred_all$lat < 6 & 
                          y_pred_all$lon > -161 & y_pred_all$lon < -149)
y_pred_bigger_box <- y_pred_all[idx_bigger_box, ]

## Rename, put into long format, and remove the standard error fields
box_diag_data <- rename(y_pred_bigger_box, LTK = LK6_pred, FSA = FSA_pred,
                        NNGP = NNGP_pred, SPDE1indep = spde_hi_pred,
                        SPDE0 = spde_pred, SPDE2 = TwoGMRF_pred) %>%
    gather(method, value, -lon, -lat, -sst, -z) %>%
    filter(!grepl("_se", method)) %>%
    filter(!method %in% c("LK4_pred", "LK5_pred"))

## Extract the actual validation data, and label as "Observed"
truth_df <- filter(box_diag_data, method == "SPDE1indep") %>%
    mutate(method = "Observed", value = z)
box_diag_data <- rbind(box_diag_data, truth_df)

## Make the method field a factor for faceting
box_diag_data$method <- factor(box_diag_data$method, 
                               levels = c("FSA", "LTK", "NNGP", 
                                          "SPDE0", "SPDE1indep", "SPDE2", "Observed"))

## Plot the facets
g1 <- ggplot(box_diag_data) + geom_point(aes(lon, lat, colour = pmax(pmin(value, 2), -2)), 
                                        size = 0.4) +
    scale_colour_gradientn(colours = nasa_palette, name = "degC") +
    facet_wrap(~method) + theme_bw(base_size = 18) + coord_fixed() +
    geom_rect(data = data.frame(ymin = -4, ymax = 4, xmin = -151, xmax = -159),
               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                            colour = "black", fill = NA)

## Save the plot
ggsave(g1, file = "../img/Z3_pred.png", width = 7, height = 7)

## Plot squared error maps
g2 <- ggplot(box_diag_data) + geom_point(aes(lon, lat, colour = pmin((z - value)^2,1)), 
                                   size = 0.4) +
    facet_wrap(~method) + theme_bw() + coord_fixed() +
    scale_colour_distiller(palette = "Spectral") +
    geom_rect(data = data.frame(ymin = -4, ymax = 4, xmin = -151, xmax = -159),
               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                            colour = "black", fill = NA)
## Save the plot
g2 <- ggsave(g2, file = "../img/Z3_spe.png", width = 7, height = 7)

