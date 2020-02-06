###########################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Plots several of the images in the paper,
##          including those of meshes and tilings.
##          This code comes with no warranty or guarantee of any kind.
###########################################################################

library("data.table")
library("dggrids")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("INLA")
library("FRK")
library("reshape2")
source("utils.R")

## Set local folder where everything is saved
cache_folder <- readLines("cache_folder.txt")

## Load meshes
load(file.path(cache_folder, "meshes.rda"))

##########################################################################
## Plot  the base diagram for illustrating the 
## effective process footprint (left panel)
##########################################################################
set.seed(1)
X <- cbind(runif(100), runif(100))
temp_mesh <- inla.mesh.2d(X, max.edge = 0.01)
png("../img/effective_footprint_base.png", 
    width = 3000, height = 3000, res = 300)
plot(temp_mesh, edge.color = rgb(0.5,0.5,0.5),
     xlim = c(0.4, 0.6), ylim = c(0.4, 0.6), main = "")
dev.off()

###################################################################
## Plot the different tilings on the sphere (Figure 1, right panel)
###################################################################
data("isea3h")
isea3h_polys <- filter(isea3h, centroid == 0 & res == 3)
isea3h_polys_rotated1 <- rotate_3d(isea3h_polys, 0, 6)
isea3h_polys_rotated2 <- rotate_3d(isea3h_polys, 7, 0)

lonlim <- 180
latlim <- 90

g <- ggplot(filter(isea3h_polys, abs(lon) < lonlim & abs(lat) < latlim)) +
  geom_path(aes(lon,lat,group=id)) +
  geom_path(data = filter(isea3h_polys_rotated1,
                          abs(lon) < lonlim & abs(lat) < latlim),
            aes(lon,lat,group=id), col = "green") +
  geom_path(data = filter(isea3h_polys_rotated2,
                          abs(lon) < lonlim & abs(lat) < latlim),
            aes(lon,lat,group=id), col = "red") +
  coord_map("ortho", orientation = c(0, 125, 65)) + theme_void()
ggsave(g, file = "../img/partitions.png", width = 5, height = 5)

###################################
## Plot the different tilings on R2
###################################
xygrid <- expand.grid(x = seq(0, 1, by = 0.1),
                      y = seq(0, 1, by = 0.1))
polys <- NULL
for(i in 1:nrow(xygrid)) {
  x1 = xygrid$x[i] - 0.05; x2 = xygrid$x[i] + 0.05;
  y1 = xygrid$y[i] - 0.05; y2 = xygrid$y[i] + 0.05;
  polys[[i]] <- data.frame(x = c(x1, x2, x2, x1, x1),
                           y = c(y1, y1, y2, y2, y1),
                           id = i)
}
polys_df <- rbindlist(polys)
g <- ggplot(polys_df) + geom_path(aes(x, y, group = id), col = "red") +
  geom_path(aes(x + 0.05, y + 0.025, group = id), col = "green") +
  geom_path(aes(x + 0.025, y + 0.05, group = id), col = "black") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
ggsave(g, file = "../img/partitionsR2.png", width = 5, height = 5)

#################################################################
## Plot the raw data globally and zoomed in around PNG (Figure 2)
#################################################################
y_tot <- read.csv(file.path(cache_folder, "y_tot.csv"))
g1 <- (ggplot(sample_n(y_tot, 100000)) +
         geom_point(aes(lon, lat, colour = pmin(pmax(z,-8),8)), pch = 46) +
         scale_colour_gradientn(colours = nasa_palette, limits = c(-8, 8),
                                name = "degC") +
         xlab("Longitude (deg)") + ylab("Latitude (deg)") +
         xlim(c(-180, 180)) + ylim(c(-90, 90)) +
         geom_rect(aes(xmin = -159, xmax = -151, ymin = -4, ymax = 4),
                   fill = NA, colour = "black", lwd = 1) +
         theme_bw()) +
  geom_map(data = map_data("world"), map = map_data("world"),
           aes(group = group, map_id=region),
           fill = "white", colour = "black", size= 0.1) +
  geom_rect(aes(xmin = 120, xmax = 160, ymin = -14, ymax = 18),
            fill = NA, colour = "red", lwd = 1) +
  coord_fixed(expand = FALSE, ratio = 1.3, ylim = c(-85, 90))

g1_zoom <- (ggplot(filter(y_tot, lon > 120 & lon < 160 &
                            lat > -14 & lat < 18)) +
              geom_point(aes(lon, lat, colour = pmin(pmax(z,-4),4)), pch = 46) +
              scale_colour_gradientn(colours = nasa_palette, limits = c(-4, 4),
                                     name = "degC")  +
              xlab("Longitude (deg)") + ylab("Latitude (deg)") +
              coord_fixed(xlim = c(120, 160),
                          ylim = c(-14, 18), expand = FALSE) +
              theme_bw()) +
  geom_map(data = map_data("world"), map = map_data("world"),
           aes(group = group, map_id=region),
           fill = "white", colour = "black", size= 0.1)

ggsave(g1, file = "../img/data.png", width = 5, height = 3)
ggsave(g1_zoom, file = "../img/data_zoom.png", width = 5, height = 3)


##############################
## Visualise meshes (Figure 3)
##############################
df_lo <- mesh3d_to_df(mesh = mesh_lo, xlim = c(0, 60), ylim = c(-45, 45))
df <- mesh3d_to_df(mesh = mesh, xlim = c(10, 35), ylim = c(-10, 15))
df_fine <- mesh3d_to_df(mesh = mesh_fine, xlim = c(14, 31), ylim = c(-6, 11))

g1 <- (ggplot() + geom_rect(data = data.frame(xmin = 15- 180, xmax= 30- 180, ymin = -5, ymax = 10),
                            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                            colour = "red", fill = NA) +
         coord_fixed(xlim = c(-180, -50)) + theme_void() + xlab("") + ylab("") +
         theme(panel.background = element_rect(fill = "white"),
               panel.border = element_blank())) %>%
  FRK::draw_world(inc_border = FALSE)

g2 <- ggplot(df_lo) + geom_path(aes(x, y, group = id)) +
  geom_path(aes(x, y, group = id)) +
  geom_path(data = df, aes(x, y, group = id), colour = "blue", alpha = 0.2) +
  geom_path(data = df_fine, aes(x, y, group = id), colour = "red", alpha = 0.15) +
  coord_fixed(xlim = c(15, 30), ylim = c(-5, 10)) + xlab("Longitude (deg)") +
  ylab("Latitude (deg)") + theme_bw()

g3 = g2 + annotation_custom(grob = ggplotGrob(g1), 
                            xmin = 15, xmax = 20, 
                            ymin = 1.5, ymax = 11.5)
ggsave(g3, file = "../img/Three_meshes_world.png", width = 5, height = 5)

########################################################
## Plot posterior inferences globally and PNG (Figure 5)
########################################################

## Extract max number of samples
samplelist <- dir(cache_folder, pattern = "*etasamp0_*") %>%
              strsplit("etasamp0_")
m <- sapply(samplelist, function(l) l[[2]]) %>% as.integer() %>% max()

## Initialise (we are thinning by a factor of 50)
etasamp1 <- etasamp0 <- NULL
samp_num <- seq(5000, m, by = 50)
samps_tot <- length(samp_num)

for(i in seq_along(samp_num)) {
  
        ## Load the samples
        load(file.path(cache_folder, paste0("etasamp0_", samp_num[i])))
        load(file.path(cache_folder, paste0("etasamp1_", samp_num[i])))
  
        ## If this is the first sample, initialise the matrix
        if(i == 1) {
            etasamp0 <- matrix(0, length(current_samp), samps_tot)
            etasamp1 <- matrix(0, length(current_xsamp1), samps_tot)
        }
  
        ## Update the matrix with the thinned samples
        etasamp0[, i] <- current_samp
        etasamp1[, i] <- current_xsamp1
}   

## Below is the function which plots the predictions, both global and
## over PNG.
do_pred_plots <- function(plotregion = c("Global", "PNG", "Malvinas")) {

  plotregion <- match.arg(plotregion)
  
  ## The base map we will use depends on whether we are global or over PNG
  if(plotregion == "Global") {
    region_xlim <- NULL
    region_ylim <- c(-85, 90)
    DPI <- 500
    mumax <- 8  
    mumin <- -8  
    semax <- 2
    ratio <- 1.3
    width <- 5
    height <- 3
  } else if (plotregion == "PNG") {
    region_xlim <- c(120, 160)
    region_ylim <- c(-14, 18)
    DPI <- 300
    mumax <- 4 
    mumin <- -4  
    semax <- 1
    ratio <- 1
    width <- 5
    height <- 3
  } else if (plotregion == "Malvinas") {
    region_xlim = c(-60, -48)
    region_ylim = c(-50, -35)
    DPI <- 500
    mumax <- 8  
    mumin <- -8  
    semax <- 2
    ratio <- 1
    width <- 3
    height <- 3
  }

  base_map0 <- plot_map(mesh, etasamp0[,1], plot = FALSE,
                        xlim = region_xlim, ylim = region_ylim)
  base_map1 <- plot_map(mesh_fine, etasamp1[,1], plot = FALSE,
                        xlim = region_xlim, ylim = region_ylim)
  
  ## Extract projections
  proj0 <- base_map0$proj
  proj1 <- base_map1$proj
  
  ## Find global field for each sample
  zsamp0 <- zsamp1 <- zsamp <- 
        matrix(0, nrow(base_map1$grid_data), ncol(etasamp1))
  
  for (i in 1:ncol(etasamp1)) {
      zsamp_mat0 <- inla.mesh.project(proj0, field = etasamp0[,i])
      zsamp_mat1 <- inla.mesh.project(proj1, field = etasamp1[,i])
      zsamp_mat <- zsamp_mat0 + zsamp_mat1
  
      zsamp0[, i] <- reshape2::melt(zsamp_mat0)[,3]
      zsamp1[, i] <- reshape2::melt(zsamp_mat1)[,3]
      zsamp[, i] <- reshape2::melt(zsamp_mat)[,3]
  }
  
  ## Find global mean and s.e.
  base_map0$grid_data$mu0 <- rowMeans(zsamp0)
  base_map0$grid_data$se0 <- sqrt(rowSums((zsamp0 - rowMeans(zsamp0))^2)/(dim(zsamp0)[2] - 1))
  base_map0$grid_data$mu1 <- rowMeans(zsamp1)
  base_map0$grid_data$se1 <- sqrt(rowSums((zsamp1 - rowMeans(zsamp1))^2)/(dim(zsamp1)[2] - 1))
  base_map0$grid_data$mu <- rowMeans(zsamp)
  base_map0$grid_data$se <- sqrt(rowSums((zsamp - rowMeans(zsamp))^2)/(dim(zsamp)[2] - 1))
  
  ## The grid spacing is not exactly EQUAL. We need to correct this to not get
  ## vertical lines when plotting
  lonaxis <- unique(base_map0$grid_data$lon)
  lataxis <- unique(base_map0$grid_data$lat)
  dlon <- round(mean(diff(lonaxis)), 4)
  dlat <- round(mean(diff(lataxis)), 4)
  coord_fix_lon <- data.frame(lon = lonaxis, lon2 = lonaxis[1] + (0:(length(lonaxis)-1))*dlon)
  coord_fix_lat <- data.frame(lat = lataxis, lat2 = lataxis[1] + (0:(length(lataxis)-1))*dlat)
  base_map0$grid_data2 <- left_join(base_map0$grid_data, coord_fix_lon) %>%
      left_join(coord_fix_lat) %>%
      mutate(lon = lon2, lat = lat2)

  grad_scale <- scale_fill_gradientn(colours = nasa_palette, limits = c(mumin, mumax),
                                     name = "degC", na.value = "white") 
  lab1 <- xlab("Longitude (deg)")
  lab2 <- ylab("Latitude (deg)")
  
  ggworld <- geom_map(data = map_data("world"), map = map_data("world"),
                      aes(group = group, map_id=region),
                      fill = "white", colour = "black", size = 0.1)
  
  g0 <- ggplot(base_map0$grid_data2) +
          geom_tile(aes(lon, lat, fill = pmin(pmax(mu0, mumin), mumax))) +
        grad_scale + lab1 + lab2 + theme_bw() + ggworld  +  
        coord_fixed(xlim = region_xlim, ylim = region_ylim,
                    ratio = ratio, expand = FALSE)
  ggsave(g0, file = paste0("../img/Ypred0_", plotregion, ".png"),
         width = width, height = height, dpi = DPI)
  
  g1 <- ggplot(base_map0$grid_data2) +
            geom_tile(aes(lon, lat, fill = pmin(pmax(mu1, mumin), mumax))) +
            grad_scale + lab1 + lab2 + theme_bw() + ggworld +
         coord_fixed(xlim = region_xlim, ylim = region_ylim,
                     ratio = ratio, expand = FALSE)
  ggsave(g1, file = paste0("../img/Ypred1_", plotregion, ".png"),
         width = width, height = height, dpi = DPI)
  
  g01 <- ggplot(base_map0$grid_data2) +
            geom_tile(aes(lon, lat, fill = pmin(pmax(mu, mumin), mumax))) +
            grad_scale + lab1 + lab2 + theme_bw() + ggworld + 
         coord_fixed(xlim = region_xlim, ylim = region_ylim, 
                     ratio = ratio, expand = FALSE)
  ggsave(g01, file = paste0("../img/Ypred01_", plotregion, ".png"),
         width = width, height = height, dpi = DPI)
  
  g01err <- ggplot(base_map0$grid_data2) +
          geom_tile(aes(lon, lat, fill = pmin(pmax(se, 0), semax))) +
        scale_fill_distiller(palette = "BrBG", limits = c(0, semax),
                             name="degC", na.value="white") +
        lab1 + lab2 + theme_bw() + ggworld +
        coord_fixed(xlim = region_xlim, ylim = region_ylim,
                    ratio = ratio, expand = FALSE)
  ggsave(g01err, file = paste0("../img/Ypred01err_", plotregion, ".png"), 
         width = width, height = height, dpi = DPI)
}

do_pred_plots("Global")
do_pred_plots("PNG")
do_pred_plots("Malvinas")
