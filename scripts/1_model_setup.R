###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 23 January 2020
## Details: Sets up the blocks for Gibbs sampling in the SST application.
##          This code comes with no warranty or guarantee of any kind.
###########################################################################

##############################################
### Part 1: SETUP
##############################################
print("Part 1: Setting up...")

## Load packages and functions
library("data.table")
library("digest")
library("dggrids") # Install from github.com/andrewzm/dggrids
library("dplyr")
library("INLA")
library("fields")
library("ggplot2")
library("igraph")
library("Matrix")
library("reshape2")
library("sparseinv")
source("utils.R")

## Set local folder where to store intermediates (data and blocks for Gibbs)
cache_folder <- readLines("cache_folder.txt")

## Setup an MD5 check so that if something is already computed it is not
## re-computed
md5_wrapper <- md5_cache(path = cache_folder)

## Plot figures when running this script?
plot_figs <- 0

## Fix seed
set.seed(1)

## Earth radius in km
true_radius_of_earth <- 6371

## Just some values with which to construct precision matrices with INLA for
## k = 0 and k = 1 scales
sigma0 <- 10
range0 <- 0.3
kappa0 <- sqrt(8 / 2) / range0
tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma0)

sigma1 <- 5
range1 <- 0.03
kappa1 <- sqrt(8 / 2) / range1
tau1 <- 1 / (sqrt(4 * pi) * kappa1 * sigma1)

## Generate mesh for theta_1
earth_points <- lonlat3D(runif(1000, -180, 180),
                         runif(1000, -90, 90))
mesh_lo <- md5_wrapper(inla.mesh.create,
                      loc = earth_points,
                      cutoff = 2000 / true_radius_of_earth,
                      refine = list(max.edge = 2400 / true_radius_of_earth))

## Generate mesh for Y_0
mesh <- gen_ocean_mesh_lo_res("../data/meshes/australia.txt",
                              "../data/meshes/antarctica.txt",
                              "../data/meshes/eurasia_africa_americas.txt",
                              max.edge = 150)

## Generate mesh for Y_1
mesh_fine <- md5_wrapper(fn = gen_ocean_mesh_lo_res,
                       aus_path = "../data/meshes/australia.txt",
                       ant_path = "../data/meshes/antarctica.txt",
                       eur_path = "../data/meshes/eurasia_africa_americas.txt",
                       max.edge = 30)

## Save meshes to local folder
save(mesh_lo, mesh, mesh_fine, file = file.path(cache_folder, "meshes.rda"))

## Plot process meshes
if (plot_figs) {
  plot_3d_mesh(mesh, col = "white")
  plot_3d_mesh(mesh_fine, new_window = FALSE, col = "yellow")
}

## Load the data, clean, and keep only unique values
load("../data/SST_sub_1000000.rda")
y_tot <- SST_sub_1000000 %>%
         filter(lat < 90) %>%
         distinct(lon, lat, .keep_all = TRUE) # remove repetitions

## Remove intercept and latitude and let's just focus on the residuals
y_tot$lat2 <- y_tot$lat^2
LM <- lm(sst ~ 1 + lat + lat2, data = y_tot)
y_tot$z <- LM$residuals

## Do the same with the validation data
load("../data/SST_sub_1000000_val.rda")
y_pred <- SST_sub_1000000_val %>%
          filter(lat < 90) %>%
          distinct(lon, lat, .keep_all = TRUE)
y_pred$lat2 <- y_pred$lat^2
coeff <- LM$coefficients
y_pred <- mutate(y_pred, z = sst - coeff[1] - coeff[2] * lat - coeff[3] * lat2)

## Generate the SPDEs and observation matrices
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
Q0 <- inla.spde.precision(spde, theta = c(log(tau0), log(kappa0)))
C0 <- inla.spde.make.A(mesh = mesh, loc = lonlat3D(y_tot$lon, y_tot$lat))
C0_pred <- inla.spde.make.A(mesh = mesh, loc = lonlat3D(y_pred$lon, y_pred$lat))

## Re,pve observations that are outside the coarse mesh and inside an
## 8x8 box in the Pacific Ocean
rm_idx <- which(rowSums(C0) == 0)
rm_idx <- c(rm_idx, which(y_tot$lon > -159 & y_tot$lon < -151 &
                            y_tot$lat > -4 & y_tot$lat < 4))
if (length(rm_idx) > 0) {
  y_tot <- y_tot[-rm_idx, ]
  C0 <- C0[-rm_idx, ]
}

## Create measurement-error precision matrix
Qobs <- Diagonal(x = pmin(1 / y_tot$error^2, 10))

## Save the processed training data to cache folder
write.csv(y_tot, file = file.path(cache_folder, "y_tot.csv"),
          row.names = FALSE)

## Remove validation data outside the main mesh
rm_idx <- which(rowSums(C0_pred) == 0)
if (length(rm_idx) > 0) {
  y_pred <- y_pred[-rm_idx, ]
  C0_pred <- C0_pred[-rm_idx, ]
}

## Save the processed validation data to cache folder
write.csv(y_pred, file = file.path(cache_folder, "y_pred.csv"),
          row.names = FALSE)

## Visualise data
if (plot_figs) {
  lo_res_world <- lo_res_mask("../data/meshes/australia.txt",
                              "../data/meshes/antarctica.txt",
                              "../data/meshes/eurasia_africa_americas.txt")
  mapWorld <- borders("world", colour = "black", fill = "light gray")
  g <- ggplot(y_tot) + geom_point(aes(lon, lat, colour = z), size = 0.01) + 
    nasa_scale("degC") + theme_bw() +
    geom_polygon(data = lo_res_world, aes(lon, lat, group = id), 
                 fill="light gray") +
    mapWorld + coord_fixed(xlim = c(-180, 180), ratio = 1.2)
}

#####################################
### Part 2: Setting up the partitions
#####################################
print("Part 2: Setting up the partitions...")

## Get lan-lot coordinates of mesh associated with Y0
lat_mesh = asin(mesh$loc[,3])*360/(2*pi)
lon_mesh = atan2(mesh$loc[,2], mesh$loc[,1])*360/(2*pi)
df_mesh <- data.frame(lat = lat_mesh, lon = lon_mesh)
points <- df_mesh[,c("lon","lat")]
names(points) <- c("x","y")
points$id <- 1:nrow(points)

## Now find neighbous from this triangulation
mesh_nei2 <- neighb_from_prec2(mesh$graph$vv)

## Partition the domain according to an ISEA3H grid at resolution 3
data("isea3h")
isea3h_centroids <- filter(isea3h, centroid == 1 & res == 3)

## Create two rotated versions to have different blocking schemes
isea3h_centroids_rotated1 <- rotate_3d(df = isea3h_centroids, east_deg = 0, 
                                       north_deg = 6)
isea3h_centroids_rotated2 <- rotate_3d(df = isea3h_centroids, east_deg = 7, 
                                       north_deg = 0) 

## Set up tilings with two rules:
## 1. Don't let blocks have too few triangles and
## 2. Don't let tiles have too few data points
## Accommodate these rules by merging tiles when needed

min_obs <- 300 # minimum observations per tile
min_tri <- 200 # minimum triangles per tile

## Our three sets of tile centroids
centroids <- list(isea3h_centroids, 
                  isea3h_centroids_rotated1,
                  isea3h_centroids_rotated2)

## For each tiling
for(i in 1:3) { 
  
  ## Allocate vertices to tile according to distance to centroids
  ## and delete centroids when conditions are not met
  ## Repeat until both condiions are satisfied
  
  OK_triangles <- OK_data <-  0
  
  while(!(OK_triangles & OK_data)) {
    D <- fields::rdist.earth(df_mesh[,c("lon","lat")],
                             centroids[[i]][,c("lon","lat")])
    points$class <- this_class <- apply(D,1,function(x) which.min(x))
    
    ## Too few triangles?
    if(all(table(this_class) >= min_tri)) {
      OK_triangles <- 1
      
      ## Too few data points?
      obs_class <- points %>% group_by(class) %>%
        summarise(obsability = sum(C0[,id]))
      if(all(obs_class$obsability >= 100)) {OK_data <- 1}
      else {
        OK_data <- 0
        rm_basis <- which(obs_class$obsability < 100)
        centroids[[i]] <- centroids[[i]][-rm_basis,]
      }
    } else {
      OK_triangles <- 0
      rm_basis <- which(table(this_class) < min_tri)
      centroids[[i]] <- centroids[[i]][-rm_basis,]
    }
  }
  points[paste0("class", i)] <- this_class
}

## For each tiling we now colour
for(i in 1:3) { 
  this_class <- points[[paste0("class", i)]] # vector of vertex tiling association
  nchunks <- length(unique(this_class))      # number of unique tilings
  
  # Relabel tiles to start from 1
  relabel <- data.frame(newclass = 1 : nchunks) 
  relabel[paste0("class", i)] <- unique(this_class)
  points <- left_join(points, relabel)
  points[paste0("class", i)] <- points$newclass
  points <- points %>% dplyr::select(-newclass)
  
  ## Plot tilings for this class
  if(plot_figs) plot_map(mesh, points$class)

  ## Now set the tiling to colour in the field class
  points$class <- points[[paste0("class", i)]]
  
  ## Find the neighbouhood list OF THE TILES
  partition_nei_list <- partition_nei_from_points(points, mesh_nei2)
  
  ## Adjacency matrix of the tiles
  adj.matrix <- adj_matrix_from_nei_list(partition_nei_list)
  
  ## Colour the tiles
  partition_colour <- colour_graph(adj.matrix,
                                   method = "backtrack",
                                   numcolours = 30,
                                   startnode = 2, obs=NULL)
  
  ## Associate the colouring with the tiling
  points <- left_join(points, partition_colour, by = "class")
  points <- points[order(points$id),]
  points[paste0("class", i)] <- points$class
  points[paste0("colour", i)] <- points$colour
  points$colour <-NULL
  points$class <- NULL
}

### Stop if we do not have four colours (we can always get 4 colours)
stopifnot(length(unique(points$colour1)) == 4)
stopifnot(length(unique(points$colour2)) == 4)
stopifnot(length(unique(points$colour3)) == 4)

## head(points)
##   x        y         id   class1 class2 class3 colour1 colour2 colour3
## 1 131.777 -1.172090  1      1      1      1       3       2       2
## 2 133.786 -0.837209  2      1      1      1       3       2       2
## 3 134.791 -2.846510  3      1      1      1       3       2       2

## Assign the blocks/colours to the fine mesh
## First find the lon-lat coordinates of the fine mesh (Y1)
## and put into a new data frame
lat_mesh_fine <- asin(mesh_fine$loc[,3])*360/(2*pi)
lon_mesh_fine <- atan2(mesh_fine$loc[,2], mesh_fine$loc[,1])*360/(2*pi)
df_mesh_fine <- data.frame(lat = lat_mesh_fine, lon = lon_mesh_fine)
points_fine <- df_mesh_fine[,c("lon","lat")]
names(points_fine) <- c("x","y")
points_fine$id <- 1:nrow(points_fine)

## Find the mapping matrix from the fine mesh to the coarse mesh
coarse_to_fine_map <- inla.spde.make.A(mesh = mesh,
                                       loc = mesh_fine$loc)
CTF_dgT <- as(coarse_to_fine_map, "dgTMatrix")
CTF_long <- data.frame(i = CTF_dgT@i + 1, j = CTF_dgT@j + 1, x = CTF_dgT@x)

## Associate the vertices of the fine triangulation with those in the coarse
## by simply finding which basis function in the coarse triangulation is the largest
CTF_alloc <- group_by(CTF_long, i) %>% dplyr::summarise(lxsens = j[which.max(x)])

## For each tiling, associate the fine-scale vertices with a tile an a colour
for(i in 1:3) {
  points_fine[paste0("class", i)] <- points[[paste0("class", i)]][CTF_alloc$lxsens]
  points_fine[paste0("colour", i)] <- points[[paste0("colour", i)]][CTF_alloc$lxsens]
}

## Find the observation matrix for the training data on the fine triangulation
C1 <- md5_wrapper(fn = inla.spde.make.A,
                  mesh = mesh_fine, 
                  loc = lonlat3D(y_tot$lon, y_tot$lat),
                  print_output = TRUE)

## Find the observation matrix for the validation data on the fine triangulation
C1_pred <- md5_wrapper(fn = inla.spde.make.A,
                       mesh = mesh_fine,
                       loc = lonlat3D(y_pred$lon, y_pred$lat),
                       print_output = TRUE)

## Save as we use them for SPDE1-indep
save(C1, C1_pred, file = file.path(cache_folder, "C1_matrices.rda"))

## Find the neighbourhood list for the fine mesh
mesh_fine_nei <- md5_wrapper(fn = neighb_from_prec2,
                             Q = mesh_fine$graph$vv)

## Construct the SPDE for the fine mesh
spde_fine <- md5_wrapper(fn = inla.spde2.matern,
                         mesh = mesh_fine,
                         alpha = 2)

## Initial precision matrix for the fine mesh
Q1 <- md5_wrapper(fn = inla.spde.precision,
                  spde = spde_fine, 
                  theta = c(log(tau0), log(kappa0)))



if(plot_figs) {
  ## Plot the three colourings to file
  g1 <- plot_map(mesh, points$colour1) + theme_bw() + 
    xlab("Longitude (deg)") + ylab("Latitude (deg)") +
    coord_fixed(ratio = 1.3, expand = FALSE, ylim = c(-85, 90)) +
    theme(text = element_text(size = 13))
  ggsave(g1, file="SSTcolouring1.png", width=5, height=3)
  
  g2 <- plot_map(mesh, points$colour2) + theme_bw() + 
    xlab("Longitude (deg)") + ylab("Latitude (deg)") +
    coord_fixed(ratio = 1.3, expand = FALSE, ylim = c(-85, 90)) +
    theme(text = element_text(size = 13))
  ggsave(g2, file="SSTcolouring2.png", width = 5, height=3)
  
  g3 <- plot_map(mesh, points$colour3) + theme_bw() + 
    xlab("Longitude (deg)") + ylab("Latitude (deg)")  +
    coord_fixed(ratio = 1.3, expand = FALSE, ylim = c(-85, 90)) +
    theme(text = element_text(size = 13))
  ggsave(g3, file="SSTcolouring3.png", width = 5, height = 3)
}

## Now we set the parallelisation schedule. Associate each wave with a colour.
ib <- block_wave_map <- list()
for(i in 1:3) {
  
  ## Block is class and wave is color (i.e., wave goes from 1 to 4)
  bl_list <- sort(unique(points[[paste0("class", i)]]))
  numblocks <- length(bl_list)
  ib[[i]] <- vector("list", numblocks)
  
  for (j in 1:numblocks) 
    ib[[i]][[j]] <-  which(points[[paste0("class", i)]] == bl_list[j])
  block_wave_map[[i]] <- as.data.frame(unique(cbind(points[[paste0("class", i)]],
                                                    points[[paste0("colour", i)]])))
  names(block_wave_map[[i]]) <- c("block","wave")
}
## Example: In the following tile 1 is updated in the third wave (colour)
##          Tile 2 in the second wave (colour) etc.
## head(block_wave_map[[1]])
##     block wave
##1     1    3
##2     2    2

## Finally we allocate observations to the tiles in the three tilings

## Do this for the training data...
print("...Seeing in which blocks observations fall...")
C0_dgT <- as(C0,"dgTMatrix")
C0_long <- data.frame(i = C0_dgT@i + 1,j= C0_dgT@j+1,x = C0_dgT@x)
obs_alloc <- group_by(C0_long,i) %>% dplyr::summarise(lxsens = j[which.max(x)])
for(i in 1:3) {
  y_tot[paste0("class", i)] <- points[[paste0("class", i)]][obs_alloc$lxsens]
}

## And the validation data...
rm_idx <- which(rowSums(C0_pred) == 0)
C0_pred <- C0_pred[-rm_idx,]
y_pred <- y_pred[-rm_idx,]
C0_pred_dgT <- as(C0_pred,"dgTMatrix")
C0_pred_long <- data.frame(i = C0_pred_dgT@i+1,j= C0_pred_dgT@j+1,x = C0_pred_dgT@x)
val_alloc <- group_by(C0_pred_long,i) %>% dplyr::summarise(lxsens = j[which.max(x)])
for(i in 1:3) {
  y_pred[paste0("class", i)] <- points[[paste0("class", i)]][val_alloc$lxsens]
}


########################################################################
### Part 3: Distributing components for sampling the PROCESS into chunks
########################################################################

## Here we distribute the data and matrices into "computational blocks"
## which are saved to local directory on disk for loading during MCMC
## It is suggested that the local folder is actually a RAM disk for quick
## loading when Gibbs sampling, but this is not necessary.

print("Part 3: Distributing components for sampling the process into chunks")
chunks_list <- list()
for(i in 1:3) {
  cat(paste0("### PROCESSING COLOURING ", i, " ###\n"))
  chunks_list[[i]] <- lapply(unique(points[[paste0("class", i)]]), function(cl) {
    cat(paste0(cl, " "))
    
    points$class <- points[[paste0("class", i)]]
    points_fine$class <- points_fine[[paste0("class", i)]]
    
    ### Start off by chunking up the data
    ### When we do a chunk, we of course need all the states in that chunk
    states_coarse <- filter(points, class == cl)
    states_coarse_id <- states_coarse$id
    
    ## We need those states in X, but also those states in other partitions that are affected
    ## by the border observations. First, find the observations affecting states in our chunk
    obs_in_chunk <- which(rowSums(C0[, states_coarse_id]) > 0)
    
    ## Do the same for the fine mesh
    states_fine <- filter(points_fine, class == cl)
    states_fine_id <- states_fine$id
    obs_in_chunk <- union(obs_in_chunk,
                          which(rowSums(C1[, states_fine_id]) > 0))
    obs_in_chunk <- sort(obs_in_chunk)
    
    ## Then see all the states that these observations affect, possible in other partitions
    affected_states_coarse <- which(colSums(C0[obs_in_chunk, ]) > 0)
    
    ## and then finding the union
    considered_states_coarse <- sort(union(affected_states_coarse, states_coarse_id))
    
    ## the externally affected states are the difference
    external_affected_states_coarse <- setdiff(considered_states_coarse, states_coarse_id)
    
    ## Uncomment the below to check the border states
    # g <- ggplot(points[states_coarse_id, ])  +
    #   geom_point(aes(x,y),colour="black") +
    #   geom_point(data = points[external_affected_states_coarse,],
    #              aes(x,y), colour="red")
    
    C0_in <- C0[obs_in_chunk, states_coarse_id]
    C0_out <- C0[obs_in_chunk, external_affected_states_coarse]
    
    ### Find all second-order neighbours (these should include all the complementary states)
    ### We need second-order because alpha = 2
    first_order_coarse_neighb <- setdiff(unique(unlist(mesh_nei2[states_coarse_id])), 
                                         states_coarse_id)
    second_order_coarse_neighb <- setdiff(unique(unlist(mesh_nei2[first_order_coarse_neighb])), 
                                          states_coarse_id)
    
    ## Do the same thing for the fine mesh
    affected_states_fine <- which(colSums(C1[obs_in_chunk, ]) > 0)
    considered_states_fine <- sort(union(affected_states_fine, states_fine_id))
    external_affected_states_fine <- setdiff(considered_states_fine, states_fine_id)
    C1_in <- C1[obs_in_chunk, states_fine_id]
    C1_out <- C1[obs_in_chunk, external_affected_states_fine]
    first_order_fine_neighb <- setdiff(unique(unlist(mesh_fine_nei[states_fine_id])), states_fine_id)
    second_order_fine_neighb <- setdiff(unique(unlist(mesh_fine_nei[first_order_fine_neighb])), 
                                        states_fine_id)
    
    ## Data in this chunk
    z01 <- y_tot$z[obs_in_chunk]
    
    ## Extract neighbouring tile indices
    all_neighb <- unique(points[external_affected_states_coarse, ]$class)
    
    ## Now subset the Q matrix
    Q1_all <- Q1[considered_states_fine, considered_states_fine]
    ext_idx <- match(external_affected_states_fine, considered_states_fine)
    Q1_in <- Q1_all[-ext_idx, -ext_idx]
    Q1_out <- Q1_all[-ext_idx, ext_idx]
    
    ## For quick construction of Q when sampling we need the following:
    P <- sparseinv:::.amd_Davis(Q1_in)
    Mall_df <- elements_in_df(spde_fine, idx = considered_states_fine)
    QC <- as(spde_fine$param.inla$M2[considered_states_fine, 
                                     considered_states_fine],
             "dgCMatrix") ## symbolic dgCMatrix
    C_allj <- cbind(C0_in, C1_in)
    Qj_in <- bdiag(Q0[states_coarse_id, states_coarse_id], Q1_in)
    Qtildejj2 <- t(C_allj) %*%  Qobs[obs_in_chunk,obs_in_chunk] %*% C_allj + Qj_in
    Ppost <- sparseinv:::.amd_Davis(Qtildejj2)
    B0 <- spde_fine$param.inla$B0[considered_states_fine, ]
    B1 <- spde_fine$param.inla$B1[considered_states_fine, ]
    B2 <- spde_fine$param.inla$B2[considered_states_fine, ]
    M0 <- spde_fine$param.inla$M0[considered_states_fine, considered_states_fine]
    M1 <- spde_fine$param.inla$M1[considered_states_fine, considered_states_fine]
    M2 <- spde_fine$param.inla$M2[considered_states_fine, considered_states_fine]

    ## Undo the above
    points$class <- NULL
    
    ## Save all items into list (one chunk)
    list(C0_in = C0_in, C0_out = C0_out,
         B0 = B0, B1 = B1, B2 = B2,
         M0 = M0, M1 = M1, M2 = M2,
         obs_idx = obs_in_chunk,
         z = z01, 
         nei_idx = all_neighb,
         Qobs = Qobs[obs_in_chunk,obs_in_chunk],
         in_idx0 = states_coarse_id, 
         out_idx0 = external_affected_states_coarse,
         all_idx0 = considered_states_coarse,
         in_idx1 = states_fine_id, 
         out_idx1 = external_affected_states_fine,
         all_idx1 = considered_states_fine,
         ext_idx1 = ext_idx,
         Q1 = Q1_all, 
         C1_in = C1_in,
         C1_out = C1_out,
         P = P, 
         Ppost = Ppost,
         Mall_df = Mall_df, 
         QC = QC)
  })
}


## Get out the indices for the fine-scale nodes by block for each tiling
## I.e., in_idx_lists[[1]][[3]] contains the index numbers of the vertices
## of the fine triangulation associated with tiling 1 and tile 3.
in_idx_lists <- list()
for(i in 1:3) {
  in_idx_lists[[i]] <- list()
  for(cl in unique(points[[paste0("class", i)]])) {
    in_idx_lists[[i]][[cl]] <- chunks_list[[i]][[cl]]$in_idx1
  }
}

## Save the individual chunks to disk and not compressing for fast load
print("...Saving individual chunks to disk...")
for(i in 1:3) {
  sapply(1:length(chunks_list[[i]]), function(j) {
    chunk <- chunks_list[[i]][[j]]
    save(chunk, file = file.path(cache_folder, paste0("chunks", j, "_colouring", i, ".rda")),
         compress = FALSE)})
}


###########################################################################
### Part 4: Distributing components for sampling the PARAMETERS into chunks
###########################################################################
print("Part 4: Distributing components for sampling the parameters into chunks")

## First, let's check how many fine scale PROCESS vertices we have associated 
## with each parameter basis function
spobj <- points_fine[c("x", "y")] %>%
  SpatialPoints(proj4string = CRS("+proj=longlat +ellps=sphere"))
S <- inla.spde.make.A(mesh_lo, 
                      lonlat3D(points_fine$x,
                               points_fine$y)) 


## There are probably some parameter basis functions that do not affect many states
## or any at all because they are on land. These can be removed.
## Note: Every state is still affected by a parameter basis function in some way
## otherwise we would get errors.
S_symb <- S
S_symb@x <- S_symb@x / S_symb@x
rm_basis_idx <- which(colSums(S_symb) < 100)
S <- S[, -rm_basis_idx]
S_symb <- S_symb[, -rm_basis_idx]

centroids_kappa_df <- lonlatfr3D(mesh_lo$loc[,1], 
                                 mesh_lo$loc[,2],
                                 mesh_lo$loc[,3])
centroids_kappa_df <- centroids_kappa_df[-rm_basis_idx, ]
centroids_kappa_df <- as.data.frame(centroids_kappa_df)

## For each parameter basis function find the neighbours and the
## associated process states (vertices)
names(centroids_kappa_df) <- c("x", "y")
kappa_affected_states <- kappa_states_neis <-
     kappa_states_incl_neis <- in_idx <- list()
nei_kappa <- list()
for(i in 1:nrow(centroids_kappa_df)) {
  tempidx1 <- which(S[, i] > 0)                       # affected states
  tempidx2 <- which(abs(colSums(Q1[tempidx1,])) > 0)  # incl. neis
  tempidx3 <- setdiff(tempidx2, tempidx1)             # neis
  tempidx4 <- intersect(which(abs(colSums(Q1[tempidx3,])) > 0), tempidx1) # boundary of "affected states"
  tempidx5 <- setdiff(tempidx1, tempidx4)  # interior states

  kappa_affected_states[[i]] <- tempidx5
  kappa_states_neis[[i]] <- tempidx4
  kappa_states_incl_neis[[i]] <- tempidx1
  in_idx[[i]] <- match(tempidx5, tempidx1)
  S_sub <- S[kappa_affected_states[[i]], ]
  nei_kappa[[i]] <- setdiff(which(colSums(S_sub) > 0), i)
}

## Add an ID to the parameter basis functions and a tile number (class)
centroids_kappa_df$id <- 1:nrow(centroids_kappa_df)
centroids_kappa_df$class <- 1:nrow(centroids_kappa_df)

## Form adjacency list and colour the tiles
adj.matrix_kappa <- adj_matrix_from_nei_list(nei_kappa)
partition_nei_list <- nei_kappa
partition_colour <- colour_graph(adj.matrix_kappa,
                                 method = "backtrack",
                                 numcolours = 30,
                                 startnode = 20, obs = NULL)
centroids_kappa_df <- left_join(centroids_kappa_df, partition_colour)
nchunks_kappa <- nrow(centroids_kappa_df)

## Uncomment the below to see an example realisation of a log parameter field
# points_fine$logkappaidw <- as.numeric(S %*% rnorm(ncol(S)))
# ggplot(points_fine[seq(1,nrow(points_fine), length.out = 100000),]) + 
#  geom_point(aes(x, y, colour = logkappaidw))

## As with the process chunks, now create a scheduling for sampling the
## parameter basis function weights
block_wave_map_kappa <- list()
for(i in 1:4)
  block_wave_map_kappa[[i]] <- 
     filter(centroids_kappa_df, colour == i)$class

## For the measurement-error standard deviation we also need to know
## which observations are associated with which basis function. This is
## found below
S_obs <- inla.spde.make.A(mesh_lo, 
                          lonlat3D(y_tot$lon,
                                   y_tot$lat))
S_obs <- S_obs[, -rm_basis_idx]
omega_affected_obs <- list()
for(i in 1:nrow(centroids_kappa_df)) {
  omega_affected_obs[[i]] <- which(S_obs[, i] > 0)
}

## Now we can put everything in chunks to be loaded when Gibbs sampling
chunks_kappa_list <- lapply(centroids_kappa_df$class, function(i) {

  ## Indices of process vertices affected by this parameter basis function
  ## including neighbours
  idx <- kappa_states_incl_neis[[i]]
  
  ## Some matrices for constructing Q, and let's put in an initial
  ## Q matrix for the states and their neighbours while we're at it
  M0 = spde_fine$param.inla$M0[idx, idx]
  M1 = spde_fine$param.inla$M1[idx, idx]
  M2 = spde_fine$param.inla$M2[idx, idx]
  Q <- inla.spde.precision.fast2(M0, M1, M2,
                                 theta = c(1, 1))

  ## Do the same as above but now excluding the neighbours
  inidx <- in_idx[[i]]
  M0_in <- M0[inidx, inidx]
  M1_in <- M1[inidx, inidx]
  M2_in <- M2[inidx, inidx]  
  Qin <- Q[inidx, inidx]  
  P = sparseinv:::.amd_Davis(Qin) # permutation matrix
  
  ## Put in the neighbours of this parameter coefficient
  neighbs = nei_kappa[[i]]

  ## Define the effective data footprint as those data points that are
  ## affected by process states that fall within the parameter basis function
  internal_data_idx <- which(rowSums(C1[, kappa_affected_states[[i]]]) > 0)
  z_kappa <- y_tot$z[internal_data_idx]
  
  ## Record those data that are affected tbe the parameter basis function 
  ## (for the measurement-error standard deviation)
  z_omega <- y_tot$z[omega_affected_obs[[i]]]

  ## Subset the corresponding matrices
  C_obs_kappa <- C1[internal_data_idx, kappa_affected_states[[i]]]
  C_obs_omega <- C1[omega_affected_obs[[i]], kappa_states_incl_neis[[i]]]
    
  ## Save into a list
  list(M0 = M0,
       M1 = M1,
       M2 = M2,
       M0_in = M0_in,
       M1_in = M1_in,
       M2_in = M2_in,
       P = P,
       inidx = inidx,
       S_sub = S[idx, -i],
       S_own = S[idx, i],
       neighbs = neighbs,
       affected_states = idx,
       affected_obs_kappa = internal_data_idx,
       affected_obs_omega = omega_affected_obs[[i]],
       z_kappa = z_kappa,
       z_omega = z_omega,
       S_obs_sub = S_obs[omega_affected_obs[[i]], -i],
       S_obs_own = S_obs[omega_affected_obs[[i]], i],
       C_obs_kappa = C_obs_kappa,
       C_obs_omega = C_obs_omega)  })


## Save to disk without compressing for fast loading
sapply(1:nchunks_kappa, function(j) {
  chunk <- chunks_kappa_list[[j]]
  save(chunk, file = file.path(cache_folder, paste0("chunks_kappa", j, ".rda")),
       compress = FALSE)})

## Save also the S matrices
save(S, S_obs, file = file.path(cache_folder, "S_matrices.rda"))


#######################################################
### Part 5: Final steps and savng of "global" variables
#######################################################
print("Part 5: Final steps and savng of 'global' variables")

## As initial condition for Y0 we can take the empirical average of the observations
## in a tile, and assign that value to each vertex in that tile :
nchunks <- length(unique(points$class1))
zmeans <- rep(0, nchunks)
for(j in 1:nchunks) {
  load(file.path(cache_folder, paste0("chunks", j, "_colouring", 1, ".rda")))
  zmeans[j] <- mean(chunk$z)
}
zmean_df <- data.frame(class1 = 1:nchunks, zmean = zmeans)
x_samp0_init <- left_join(points, zmean_df)$zmean

## As initial condition (first sample) of Y1 we can just use 0
## Actually both these initial conditions will be modulated when using
## multiple chains.
current_xsamp1 <- rep(0, mesh_fine$n)

## Now we save the matrices and vectors needed for running the Gibbs sampler 
## that are not associated with a specific tiling
print("...Saving chunks and information for HPC sampling...")
B0 = spde$param.inla$B0
B1 = spde$param.inla$B1
B2 = spde$param.inla$B2
Mall_df0 <- elements_in_df(spde)
QC0 <- as(spde$param.inla$M2, "dgCMatrix") ## symbolic dgCMatrix
nobs <- nrow(y_tot)
z_all <- y_tot$z
save(points, mesh, spde, Q0, Mall_df0, QC0, C0, z_all,
     block_wave_map, current_xsamp1, x_samp0_init,
     tau0, kappa0, tau1, kappa1, nobs, S_obs,
     B0, B1, B2, in_idx_lists, points_fine, C0_pred,
     nchunks_kappa, block_wave_map_kappa, S,
     file = file.path(cache_folder, "data_for_SST2M_no_chunks.rda"))

##############################################
### Part 6: Sampling
##############################################
## From here on run on the HPC
## Can be run independently from all of the above
print("Data ready for MCMC. Please run 2_ocean_sampling.R")


