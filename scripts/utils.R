###############################################################################
## Author: Andrew Zammit-Mangion
## Date: 24 January 2020
## Details: Internal functions for the SPDE2 model
##          This code comes with no warranty or guarantee of any kind.
###############################################################################

plot_map <- function(mesh,field,plot=TRUE,bigwidth = 3*361,bigheight =  3*181, leg_title ="", xlim = NULL, ylim = NULL) {

    lonlatpts <- lonlatfr3D(mesh$loc[,1],mesh$loc[,2],mesh$loc[,3])
   if(is.null(xlim) | is.null(ylim)) {
        proj = inla.mesh.projector(mesh, dims=c(bigwidth, bigheight))
        lon_lat_grid <- expand.grid(lon = seq(-180, 180, length = bigwidth),
                                    lat = seq(-90, 90, length = bigheight))
   } else {
       proj = inla.mesh.projector(mesh, dims=c(bigwidth, bigheight), xlim = xlim, ylim = ylim)
        lon_lat_grid <- expand.grid(lon = seq(xlim[1], xlim[2], length = bigwidth),
                                    lat = seq(ylim[1], ylim[2], length = bigheight))
    }
    zproj <- inla.mesh.project(proj, field = field)
    grid_data <- reshape2::melt(zproj)
    names(grid_data) <- c("lon","lat","z")
    grid_data[,1:2] <- lon_lat_grid

    XX <- filter(grid_data,!is.na(z))
    minx <- min(XX[,1])
    maxx <- max(XX[,1])
    miny <- min(XX[,2])
    maxy <- max(XX[,2])


    if(plot) {
        mapWorld <- borders("world", colour="black", fill="light gray")
        return(ggplot(grid_data) +
        geom_raster(aes(lon,lat,fill=z)) +
        scale_fill_distiller(palette="Spectral",name=leg_title,na.value="light gray") +
        coord_fixed(xlim=c(minx,maxx), ylim=c(miny,maxy),ratio=1.2) + mapWorld + theme_bw())
    } else {
        return(list(grid_data = grid_data,proj = proj))
    }
}

MHadapt <- function(sd, samps) {
    acceptrate <- (length(unique(samps)) - 1) / length(samps)
    if(acceptrate > 0.8) sd <- sd*2
    if(acceptrate < 0.1) sd <- sd/2
    sd        
}


plot_3d_mesh <- function(mesh,new_window=TRUE,col="white",vertex.color="red",size=2) {
    if(new_window) {
        require(rgl)
        windowRect.globe = c(50,50,50+840,50+400)
        open3d(windowRect=windowRect.globe)
    }
    plot(mesh, rgl=TRUE, col=col,
         draw.vertices=TRUE, draw.edges=TRUE, add=TRUE,
         vertex.color=vertex.color,size=size)
}


mesh3d_to_df <- function(mesh, xlim = c(0, 360), ylim = c(-90, 90)) {
  c2d <- lonlatfr3D(mesh$loc[,1],mesh$loc[,2],mesh$loc[,3])
  tv <- mesh$graph$tv
  
  rmidx <- which(c2d[,1] < xlim[1] | c2d[,1] > xlim[2] |
                   c2d[,2] < ylim[1] | c2d[,2] > ylim[2])
  
  rmtv1 <- which(tv[,1] %in% rmidx)
  rmtv2 <- which(tv[,2] %in% rmidx)
  rmtv3 <- which(tv[,3] %in% rmidx)
  rmtv <- union(union(rmtv1, rmtv2), rmtv3)
  tv <- tv[-rmtv, ]
  
  l <- lapply(1:nrow(tv), function(i) {
    data.frame(x = c(c2d[tv[i, ], 1],
                     c2d[tv[i, 1], 1]),
               y = c(c2d[tv[i, ], 2],
                     c2d[tv[i, 1], 2]),
               id = i)}
  )
  rbindlist(l)
}

mesh2d_to_df <- function(mesh, xlim = NULL, ylim = NULL) {
  c2d <- mesh$loc
  tv <- mesh$graph$tv
  
  if(!is.null(xlim)) {
    rmidx <- which(c2d[,1] < xlim[1] | c2d[,1] > xlim[2] |
                     c2d[,2] < ylim[1] | c2d[,2] > ylim[2])
    
    rmtv1 <- which(tv[,1] %in% rmidx)
    rmtv2 <- which(tv[,2] %in% rmidx)
    rmtv3 <- which(tv[,3] %in% rmidx)
    rmtv <- union(union(rmtv1, rmtv2), rmtv3)
    tv <- tv[-rmtv, ]
  }
  l <- lapply(1:nrow(tv), function(i) {
    data.frame(x = c(c2d[tv[i, ], 1],
                     c2d[tv[i, 1], 1]),
               y = c(c2d[tv[i, ], 2],
                     c2d[tv[i, 1], 2]),
               id = i)}
  )
  rbindlist(l)
}

lonlat3D=function(lon,lat){
    cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
          sin((lon/180)*pi)*cos((lat/180)*pi),
          sin((lat/180)*pi))
}

insert_elements <- function(x, idx, vals) {
  lorig <- length(x)
  x <- c(x, rep(0, length(vals)))
  for(i in seq_along(idx)) {
    l <- lorig + i - 1
    if(idx[i] < length(x)) { # Shift
      x[(idx[i]+1):(l + 1)] <- x[idx[i]:l]
    }
    x[idx[i]] <- vals[i]
  }
  x
}

lonlatfr3D <- function(x,y,z) {
    # lat = asin(z)*360/(2*pi)
    # lon = atan2(y, z)*360/(2*pi)
    lat = acos(z)*360/(2*pi)
    lon = atan2(y,x)*360/(2*pi)+180
    cbind(lon,(90-lat))
}

neighb_from_prec2 <- function(Q) {
    Q <- as(Q,"dgTMatrix")
    df <- data.frame(i=Q@i+1,j=Q@j+1) %>% group_by(i) %>% dplyr::summarise(j=list(j))
    i_full <- data.frame(i=1:nrow(Q))
    left_join(i_full,df)$j
}





logf_marg0 <- function(theta, mu_logsigma, sd_logsigma,
                       mu_logrho, sd_logrho, spde,
                       P,  z, Qobs, C_obs) {

  ## Here theta = (log_sigma,log_rho)
  log_sigma <- theta[1]
  log_rho <- theta[2]
  log_kappa <- log_rho2log_kappa(log_rho, 1)
  log_tau <- log_sigma2log_tau(log_sigma, log_kappa, 1)
  
  if(log_sigma < (mu_logsigma - 50*sd_logsigma) |
     log_sigma > (mu_logsigma + 50*sd_logsigma) |
     log_rho < (mu_logrho - 50*sd_logrho) |
     log_rho > (mu_logrho + 50*sd_logrho))
  {
    -Inf
  } else {
      
      LQobs <- chol(Qobs)
      logdetQobs <- logdet(LQobs)
      M0 <- spde$param.inla$M0
      M1 <- spde$param.inla$M1
      M2 <- spde$param.inla$M2

      Qx = inla.spde.precision.fast2(M0, M1, M2, 
                                     theta = cbind(log_tau,
                                                   log_kappa))
      LQx_perm <- t(Matrix::chol(Qx[P, P]))
      Pmat <- sparseMatrix(i = P, j = 1:length(P), x = 1)
      logdetQx <- logdet(LQx_perm)
      
      Qstar <- crossprod(t(LQobs) %*% C_obs) + Qx
      LQstar_perm <- t(Matrix::chol(Qstar[P, P]))
      logdetQstar <- logdet(LQstar_perm)
      mustar <- t(C_obs) %*% (Qobs %*% z)
      Qstarinv_mustar <- cholsolve(Q = Qstar, y = mustar, 
                                   perm = TRUE, cholQp = LQstar_perm, P = Pmat)
      xsamp <- as.numeric(Qstarinv_mustar + Pmat %*% solve(t(LQstar_perm),t(Pmat) %*% rnorm(length(mustar))))
      
      list(logf = 0.5 * as.numeric(
                            logdetQobs + logdetQx - logdetQstar - 
                            crossprod(t(LQobs) %*% z) +
                            t(mustar) %*% Qstarinv_mustar -
                            (log_sigma - mu_logsigma)^2/(sd_logsigma)^2 -
                            (log_rho - mu_logrho)^2/(sd_logrho)^2),
           xsamp = xsamp)
  }
}


logf_marg1 <- function(theta, theta_surface, mu_logsigma,
                      sd_logsigma, mu_logrho,sd_logrho, 
                      M0, M1, M2, M0_in, M1_in, M2_in, Yx, P, inidx,
                      z, Qobs, C_obs) {
  ## Here theta = (log_sigma,log_rho)
  log_sigma <- theta[1]
  log_rho <- theta[2]
  log_kappa <- log_rho2log_kappa(log_rho, 1)
  log_tau <- log_sigma2log_tau(log_sigma, log_kappa, 1)
  
  log_kappa_surface <- log_rho2log_kappa(theta_surface[, 2], 1)
  log_tau_surface <- log_sigma2log_tau(theta_surface[, 1], log_kappa_surface, 1)
  
  if(log_sigma < (mu_logsigma - 50*sd_logsigma) |
     log_sigma > (mu_logsigma + 50*sd_logsigma) |
     log_rho < (mu_logrho - 50*sd_logrho) |
     log_rho > (mu_logrho + 50*sd_logrho)) {
    -Inf
  } else {

      LQobs <- chol(Qobs)
      logdetQobs <- logdet(LQobs)
      
      Qx = inla.spde.precision.fast2(M0, M1, M2, 
                                     theta = cbind(log_tau_surface,
                                                   log_kappa_surface))
      Qx_in = inla.spde.precision.fast2(M0_in, M1_in, M2_in, 
                                        theta = cbind(log_tau_surface[inidx],
                                                      log_kappa_surface[inidx]))
      LQx_perm <- t(Matrix::chol(Qx_in[P, P]))
      Pmat <- sparseMatrix(i = P, j = 1:length(inidx), x = 1)
      logdetQx <- logdet(LQx_perm)
      nei_contrib <- cholsolve(Q = Qx_in, y = Qx[inidx, -inidx] %*% Yx[-inidx], 
                               perm = TRUE, cholQp = LQx_perm, P = Pmat)
      
      Qstar <- crossprod(t(LQobs) %*% C_obs) + Qx_in
      LQstar_perm <- t(Matrix::chol(Qstar[P, P]))
      logdetQstar <- logdet(LQstar_perm)
      mustar <- t(C_obs) %*% (Qobs %*% z) - Qx_in %*% nei_contrib
      Qstarinv_mustar <- cholsolve(Q = Qstar, y = mustar, 
                                   perm = TRUE, cholQp = LQstar_perm, P = Pmat)
      xsamp <- as.numeric(Qstarinv_mustar + Pmat %*% solve(t(LQstar_perm),
                                                           t(Pmat) %*% rnorm(length(mustar))))

      logf <- 0.5 * as.numeric(
                            logdetQobs + logdetQx - logdetQstar - 
                            t(nei_contrib) %*% Qx_in %*% nei_contrib +
                            t(mustar) %*% Qstarinv_mustar -
                            (log_sigma - mu_logsigma)^2/(sd_logsigma)^2 -
                            (log_rho - mu_logrho)^2/(sd_logrho)^2)

      if(length(z) > 0) logf <- logf -  as.numeric(crossprod(t(LQobs) %*% z))
                           
      
      list(logf = logf,
           xsamp = xsamp)

  }
}

logf_sigma_eps <- function(theta, 
                         theta_surface,
                         mu_lsigma_eps,
                         sd_lsigma_eps,
                         C_obs,
                         Yx, 
                         z) {
  Qobs <- Diagonal(x = 1/(exp(theta_surface))^2)
  GAMMA <- as.numeric(t(z - C_obs %*% Yx) %*% 
                        Qobs %*% 
                        (z - C_obs %*% Yx))
  0.5*logdet(chol(Qobs)) - 0.5*GAMMA -
    0.5*(theta - mu_lsigma_eps)^2/(sd_lsigma_eps^2)
}


logdet <- function (L) 
{
  diagL <- diag(L)
  return(2 * sum(log(diagL)))
}

## extract elements from spde M matrices into data frame. Note: M2 non zeroes are a superset of the M0 and M1 non zeroes
elements_in_df <- function(spde, idx = NULL) {
  if(is.null(idx)) {
    M0 <- spde$param.inla$M0
    M1 <- spde$param.inla$M1
    M2 <- spde$param.inla$M2
  } else {
    M0 <- spde$param.inla$M0[idx, idx]
    M1 <- spde$param.inla$M1[idx, idx]
    M2 <- spde$param.inla$M2[idx, idx]
  }

    M0_df <- data.frame(i = M0@i+1, j = M0@j + 1, x0 = M0@x)
    M1_df <- data.frame(i = M1@i+1, j = M1@j + 1, x1 = M1@x)
    M2_df <- data.frame(i = M2@i+1, j = M2@j + 1, x2 = M2@x)
    M_alldf <- left_join(left_join(M2_df, M1_df,
                                   by=c("i","j")),
                         M0_df, by = c("i","j"))
    M_alldf$x0[which(is.na(M_alldf$x0))] <- 0
    M_alldf$x1[which(is.na(M_alldf$x1))] <- 0
    M_alldf
}

inla.spde.precision.fast2 <- function (M0, M1, M2, theta)
{
  if(length(theta) > 2) {
    tau <- exp(theta[, 1])
    kappa2 <- exp(theta[, 2])^2
    D0 <- Diagonal(x = tau)
    D1 <- Diagonal(x = kappa2)
    Q = (D0 %*% (D1 %*% M0 %*% D1 +
                   D1 %*% M1 + t(M1) %*% D1 + M2) %*% D0)
  } else {
    tau <- exp(theta[1])
    kappa2 <- exp(theta[2])^2
    Q = tau^2 * (kappa2^2 * M0 + 2 * kappa2 * M1 + M2)
  }
  return(as(Q, "dgCMatrix"))
}



## Generates an ocean mesh -- requires INLA
gen_ocean_mesh_hi_res <- function(coast_path,antarctica_path,
                           cutoff=100,max.edge=150) {

    true.radius.of.earth = 6371
    radius.of.earth = 1

    coastline <- readShapeSpatial(coast_path)
    antarctica <- readShapeSpatial(antarctica_path)
    poly_2_segments <- function(polys) {
        dfs <- lapply(polys@polygons,
                      function(l) l@Polygons[[1]]@coords)
        lapply(dfs,function(l) {
            coords3d <- lonlat3D(l[,1],l[,2])
            inla.mesh.segment(loc=coords3d)
        })
    }
    segm_coastlines <- poly_2_segments(coastline)
    segm_antarctica <- poly_2_segments(antarctica)

    inla.mesh.create(boundary=c(segm_coastlines[c(1:500)],
                                segm_antarctica),
                     cutoff=cutoff/true.radius.of.earth,
                     refine=list(max.edge=max.edge/true.radius.of.earth))
}


## Returns data frame
lo_res_mask <- function(aus_path,ant_path,eur_path) {

  aus <- read.table(aus_path)
  aus = rbind(aus,aus[1,])
  ant <- read.table(ant_path)
  ant <- rbind(ant,ant[1,])
  eur <-  read.table(eur_path)
  eur <- rbind(eur,eur[1,])
  colnames(aus) <- colnames(ant) <- colnames(eur) <- c("lon","lat")
  aus$id <- 1
  ant$id <- 2
  america <- eur[45:126,]
  eurasia <- eur[-c(45:126),]
  america$id <- 3
  eurasia$id <- 4
  df <-   rbind(aus,ant,america,eurasia)
  df$lon[df$lon > 180] <- df$lon[df$lon > 180] - 360
  df
}

gen_ocean_mesh_lo_res <- function(aus_path,ant_path,eur_path, 
                                  max.edge = 150, cutoff = 100) {

    true.radius.of.earth = 6371
    
    aus <- read.table(aus_path)
    aus = rbind(aus,aus[1,])
    ant <- read.table(ant_path)
    ant <- rbind(ant,ant[1,])
    eur <-  read.table(eur_path)
    eur <- rbind(eur,eur[1,])
    #Make everything 3D
    aus3 <- lonlat3D(aus[,1],aus[,2])
    ant3 <- lonlat3D(ant[,1],ant[,2])
    eur3 <- lonlat3D(eur[,1],eur[,2])
    segm_aus = inla.mesh.segment(loc=aus3)
    segm_ant = inla.mesh.segment(loc=ant3)
    segm_eur = inla.mesh.segment(loc=eur3)

    segm_world = list(segm_aus,segm_ant,segm_eur)
    mesh=inla.mesh.create(boundary=segm_world,
                          cutoff=cutoff/true.radius.of.earth,
                          refine=list(max.edge=max.edge/true.radius.of.earth,
                                      extend=list(offset=-0.5)))
}


INLA_spatial_estimation <- function(spde, z, C, C_pred) {
    
    s_index <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = spde$n.spde)
    ## First stack: the model for estimation
    stack_est <- inla.stack(data = list(z = z),
                            A = list(C),
                            effects = list(c(s_index,list(Intercept=1))),
                            tag="est")
    ## Second stack: the model for prediction
    stack_pred <- inla.stack(data = list(z = NA),
                             A = list(C_pred),
                             effects = list(c(s_index,list(Intercept=1))),
                             tag="pred")
    ## Now combine
    stack <- inla.stack(stack_est,stack_pred)

    ## Formula
    formula <- z ~ -1 + Intercept +
        f(spatial.field,
          model=spde)

    output <- inla(formula,
                   data=inla.stack.data(stack,spde=spde),
                   family="gaussian",
                   control.predictor = list(A=inla.stack.A(stack),
                                            compute=FALSE),
                   control.results=list(return.marginals.random=FALSE,
                                        return.marginals.predictor=FALSE))
}


INLA_spatial_estimation2 <- function(spde,z,C) {
  s_index <- inla.spde.make.index(name = "spatial.field",
                                  n.spde = spde$n.spde)
  ## First stack: the model for estimation
  stack_est <- inla.stack(data = list(z = z),
                          A = list(C),
                          effects = list(c(s_index)),
                          tag="est")

  stack <- stack_est

  ## Formula
  formula <- z ~ -1 + f(spatial.field,model=spde)

  output <- inla(formula,
                 data=inla.stack.data(stack,spde=spde),
                 family="gaussian",
                 control.predictor = list(A=inla.stack.A(stack),
                                          compute=FALSE),
                 control.results=list(return.marginals.random=FALSE,
                                      return.marginals.predictor=FALSE))
}


if(require(ggplot2)) {
  nasa_palette <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                  "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d","#f6da0c","#f5a009",
                  "#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")
  nasa_scale <- function(leg_title = " ") scale_color_gradientn(colours = nasa_palette, name = leg_title)
}


rotate_3d <- function(df, east_deg = 0, north_deg = 0) {
  if(east_deg != 0) {
    df$lon <- df$lon + east_deg/cos(df$lat * 2 * pi / 360)
    
    idx <- which(df$lon > 180)
    if(length(idx) > 0)
      df$lon[idx] <- df$lon[idx] - 360
    
    idx <- which(df$lon < -180)
    if(length(idx) > 0)
      df$lon[idx] <- -df$lon[idx] + 360
  }
  
  if(north_deg != 0) {
    coords3d <- lonlat3D(df$lon, df$lat)
    
    ## Now flip
    temp <- coords3d[, 1]
    coords3d[, 1] <- coords3d[, 3]
    coords3d[, 3] <- temp
    
    ## Convert back to lat-lon and shift
    coords2d <- lonlatfr3D(coords3d[, 1],coords3d[, 2], coords3d[, 3])
    coords2d[,1] <- coords2d[,1] - 180
    coords2d[, 1]  <- coords2d[, 1] + north_deg / cos( coords2d[, 2] * 2 * pi / 360)
    
    idx <- which(coords2d[, 1] > 180)
    if(length(idx) > 0)
      coords2d[idx, 1] <- coords2d[idx, 1] - 360
    
    idx <- which(coords2d[, 1] < -180)
    if(length(idx) > 0)
      coords2d[idx, 1] <- -coords2d[idx, 1] + 360
    
    ## Convert to 3d and flip back
    coords3d <- lonlat3D(coords2d[, 1], coords2d[, 2])
    
    ## Now flip
    temp <- coords3d[, 1]
    coords3d[, 1] <- coords3d[, 3]
    coords3d[, 3] <- temp
    
    ## Convert back to 2D
    coords2d <- lonlatfr3D(coords3d[, 1], coords3d[, 2], coords3d[, 3])
    df$lon <- coords2d[, 1] - 180
    df$lat <- coords2d[, 2]
  }
  
  df
}


flog_console_file <- function(fname = "run.log", loggername = "logger") {
    flog.appender(appender.file(fname), loggername)
    function(s) {
        flog.info(s)
        flog.info(s, name = loggername)
    }
}

log_rho2log_kappa <- function(log_rho, nu) {
  log(8 * nu)/2 - log_rho
}

log_sigma2log_tau <- function(log_sigma, log_kappa, nu, d = 2) {
  alpha <- nu + d/2
  0.5*log(gamma(nu) / (gamma(alpha) * (4*pi)^(d/2))) - log_sigma - nu*log_kappa
}

lognorm_cost <- function(par, lx, ux) {
    mu <- par[1]
    sigma <- par[2]
    l <- qlnorm(0.025, meanlog = mu, sdlog = sigma)
    u <- qlnorm(0.975, meanlog = mu, sdlog = sigma)
    sum(((l - lx)/lx)^2 + ((u - ux)/ux)^2)
    }


circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#' Check md5 sum in cache directory
#'
#' @description This function takes an md5 and checks whether it is in the default cache directory
#' @param md5 the md5 checksum
#' @param path the directory where the cache is stored
#' @return True or False, depending on whether the file exists in the cache.
#' @export
#' @examples
#' require(digest)
#' md5 <- digest(2)
#' check_md5(md5,".")
check_md5 <- function(md5,path) {
  return(file.exists(file.path(path,paste0(md5,".rda"))))
}

#' md5 checksum function
#'
#' Creates a wrapper for generating the md5 checksum for the calling function and its arguments. If a file with the md5 as its filename exists in the specified folder then this is loaded. Otherwise the function is evaluated and the results are stored in the specified folder.
#' @param path the path to where the files will be cached. 
#' @return An md5 wrapper for intensive computations. This has arguments
#' \itemize{
#' \item \code{fn}: the function being called
#' \item \code{...}: the arguments to be passed to the function
#' \item \code{print_output}: if T, details of md5 digest are outputted to screen.
#' }
#' @details md5 checksums on the called function and any function-arguments are generated by creating a text file in the specified folder and digesting 
#' the text file before deleting it. md5 checksums on character arguments are carried out on the file (if the file exists) or on the character as appropriate.
#' @keywords md5, digest
#' @export
#' @examples
#' myfun1 <- function(x) x^2
#' myfun2 <- function(y,f) y^3 + f(y)
#' md5_wrapper <- md5_cache(".")
#' x <- md5_wrapper(myfun2,2,myfun1)
md5_cache <- function(path) {
  stopifnot(is.character(path))
  dir.create(file.path(path), showWarnings = FALSE)
  md5 <- NA
  
  function(fn,...,print_output=F) {
    args <- list(...)
    if (length(args)==0) stop("Nothing to check md5 against")
    md5 <<- md5_fun(fn,path)
    lapply(args,function(x) {  
      if(is.character(x) & length(x) == 1) {
        if(file.exists(x)) {
          md5_temp <- digest(file=x)
        } else {
          md5_temp <- digest(x)
        }
      } else {
        if(class(x) == "function") {
          md5_temp <- md5_fun(x,path)
        } else {
          md5_temp <- digest(x)
        }
      }
      md5 <<- digest(c(md5,md5_temp))
    })
    
    if(check_md5(md5,path)) {
      if(print_output) cat("Found existing MD5 for data. Loading...",sep="\n")
      load(file=file.path(path,paste(md5,".rda",sep=""))) 
    } else {
      cat("Data not found. Evaluating expression...",sep="\n")
      flush.console()
      X <- fn(...)
      save(X,file=file.path(path,paste0(md5,".rda")))
    }
    if(print_output) cat(paste("md5 for this result is ",md5,sep=""),sep="\n")
    return(X)  
  }
}


#' Generate md5 checksum from function
#'
#' Creates an md5 checksum of a function by dumping it and then digesting the file.
#' @param fn the function to digest
#' @param path the directory where dumping will be carried out
#' @details Note that a cache folder needs to be present for this to work.
#' @keywords md5, digest
#' @export
#' @examples
#' myfun1 <- function(x) x^2
#' x <- md5_fun(myfun1,".")
md5_fun <- function(fn,path) {
  path <- file.path(path,"fun_check.R")
  dump(ls(pattern = 'fn'),path)
  #save(fn,file=path)
  fun_digest <- digest::digest(file=path)
  unlink(path)
  return(fun_digest)
}


#' @title Find cpartition neighbourhood from underlying graph
#' @author Andrew Zammit Mangion
#' @description This function takes a set of points with a field \code{class} and the neighbourhood list of those points,
#' and computes the neighbourhood list of the super-graph composed by treating the \code{class} as a vertex in its own right. If even
#' one underlying vertex has a nerghbour in another partition, then that partition is the present partition's neighbour.
#' @param points a data frame with fields \code{x}, \code{y} and \code{class}
#' @param nei_list a list where each element is a vector containing the indices of the neighbours of the corresponding node
#' @export
#' @return a neighbourhood list of the superset
partition_nei_from_points <- function(points,nei_list) {
  partition_nei_list <- plyr::dlply(points,"class", function(df) {
    classes <- c()
    for (i in df$id) {
      points_to_examine <- nei_list[[i]]
      classes <- c(classes,points$class[points_to_examine])
    }
    classes <- classes[-which(classes == df$class[1])]
    return(unique(classes))
  })
  return(partition_nei_list)
}


#' @title Return adjacency matrix from list of neighbours
#' @author Andrew Zammit Mangion
#' @description Creates a sparse adjacency matrix from list of vertex neighbours.
#' @param nei_list a list of vectors, where each vector contains the indices of the adjacent neighbours.
#' @return a sparse adjacency matrix (with zero on the diagonal)
#' @export
#' @examples
#' nei_list <- list(c(2,3),1,1)
#' adj.matrix <- adj_matrix_from_nei_list(nei_list)
adj_matrix_from_nei_list <- function(nei_list) {
  
  adj.matrix <- sparseMatrix(i= unlist(sapply(1:length(nei_list),function(i) {rep(i,length(nei_list[[i]]))})),
                             j = as.vector(unlist(nei_list)),
                             x=1)
  return(adj.matrix)
}


#' @title Greedy graph colouring algorithm
#' @author Andrew Zammit Mangion
#' @description This algorithm colours the nodes in sequential order, where the order is given through a breadth first search algorithm started on a
#' pre-specified start-node.
#' @param adj.matrix the graphs adjacency matrix (usually sparse)
#' @param numcolours the maximum number of colours to try out with the \code{BFS} and \code{DSATUR} methods. The function gives an error if this is exceeded.
#' @param method takes values \code{"BFS"}, \code{"DSATUR"} or \code{"backtrack"} for the breadth-first-search, the maximum-degree-of-saturation and the backtracking respectively. 
#' In the latter method, a four-colour configuration is attempted using brute force, where the order is obtained by first running a \code{DSATUR} run.
#' @param obs an \code{m x n} matrix identifiying which partitions each observations affect. If present, the algorithm
#' will not produce a colouring where an observation influences two or more partitions of the same colour
#' (as such this increases the chromatic number of the problem). This option is only available with \code{BFS} and 
#' \code{DSATUR}
#' @param startnode the starting node for the BFS algorithm
#' @return a data frame with the associated colour for each vertex (also \code{class})
#' @references \url{http://community.topcoder.com/longcontest/?module=Static&d1=match_editorials&d2=intel_mtcs_10}
#' @export
colour_graph <- function(adj.matrix,numcolours,method="BFS",startnode=1,obs=NULL) {
  partition_colour <- data.frame(class = 1:nrow(adj.matrix),colour=0)
  G <- graph.adjacency(adj.matrix)
  
  if(method == "BFS") {
    G.bfs <- graph.bfs(G,startnode)
    for(i in G.bfs$order) {
      # check neighbour colours
      # use minimum of intersection of 1:numcolours and neighbour colours as this colour
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      if (suppressWarnings(any(nei_colours))) {
        if(!(is.null(obs))) {
          offending_obs <- which(obs[,i] == 1)
          affected_partitions <- Reduce("union",apply(D[offending_obs,],1,function(x) which(x==1)))
          prohibited_colours <- partition_colour$colour[affected_partitions]
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),union(nei_colours,prohibited_colours)))
        } else {
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),nei_colours)) 
        }
      } else {
        partition_colour$colour[i] <- 1 
      }  
    }
  }
  
  
  if(method %in% c("DSATUR","backtrack")) {
    degrees <- rowSums(adj.matrix)
    i <- which.max(degrees)
    col.matrix <- adj.matrix
    done_list <- NULL
    while(any(partition_colour$colour == 0)) {
      done_list <- c(done_list,i)
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      ## Color node
      if (suppressWarnings(any(nei_colours))) {
        if(!(is.null(obs))) {
          offending_obs <- which(D[,i] == 1)
          affected_partitions <- Reduce("union",apply(D[offending_obs,],1,function(x) which(x==1)))
          prohibited_colours <- partition_colour$colour[affected_partitions]
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),union(nei_colours,prohibited_colours)))
        } else {
          partition_colour$colour[i] <- min(setdiff(c(1:numcolours),nei_colours))
        }
      } else {
        partition_colour$colour[i] <- 1 
      }  
      
      if(length(done_list) < nrow(adj.matrix)) {
        
        col.matrix[which(col.matrix[,i]>0),i] <-  partition_colour$colour[i]+1
        ## Find next node
        sat.level <- apply(col.matrix,1,function(x) length(unique(x)))
        maxsat <- max(sat.level)
        found_node <- 0
        while(!found_node) {
          vv <- setdiff(which(sat.level == maxsat),done_list)
          if(length(vv) == 0) {
            maxsat <- maxsat - 1
          } else {
            i <- vv[which.max(degrees[vv])]
            found_node <- 1
          }
        }
      }  
    }
  }
  
  if(method=="backtrack") { 
    numcolours=4
    partition_colour <- data.frame(class = 1:nrow(adj.matrix),colour=0)
    order <- done_list
    done <- 0
    i <- order[1]
    ilist <- NULL
    tried_colours <- list(); length(tried_colours) <- nrow(adj.matrix)
    program_state <- list()
    global_count <- 0
    while (!done) {
      backtrack = 0
      nei_colours <- partition_colour$colour[partition_nei_list[[i]]]
      if (suppressWarnings(any(nei_colours))) {
        colours_to_try <- setdiff(1:numcolours,tried_colours[[i]])
        if (length(setdiff(colours_to_try,nei_colours)) ==0) {
          ## BACKTRACK
          X <- match(partition_nei_list[[i]],ilist) # return first positions of neighbours in ilist
          go.back.to <- max(X,na.rm=T)
          i <- ilist[go.back.to]
          #restore state
          tried_colours <-  program_state[[go.back.to]]$tried_colours
          partition_colour <-  program_state[[go.back.to]]$partition_colour
          global_count <- go.back.to - 1
          backtrack = 1
          ilist <- ilist[1:go.back.to-1] ##### TEST THIS LINE
        } else {
          colour <- min(setdiff(colours_to_try,nei_colours))  
          if(global_count + 1 == nrow(adj.matrix)) done = 1
        }
        
      } else {
        colour <- min(setdiff(c(1 :numcolours),tried_colours[[i]]))
      }  
      
      if(!backtrack) {
        global_count <- global_count + 1
        ilist <- c(ilist,i)
        partition_colour$colour[i] <- colour
        tried_colours[[i]] <- c(tried_colours[[i]],colour)
        i <- order[which(order == i) + 1]
        program_state[[global_count]] <- list(tried_colours = tried_colours, partition_colour = partition_colour) 
      }
      print(i)
    }
  }
  return(partition_colour)
}



#' @title Find convex hull from a set of points
#' @author Andrew Zammit Mangion
#' @description Returns the convex hull from a set of points with fields \code{x} and \code{y}
#' @param points a data frame with fields \code{x} and \code{y}
#' @export
#' @return subset of \code{points} corresponding to the convex hull
#' @examples
#' N = 99
#' points <- data.frame(x = runif(n=N), y = runif(n=N),id=1:N)
#' hull <- find_hull(points)
#' plot(points$x,points$y)
#' lines(hull$x,hull$y,col="red")
find_hull <- function(points) {
  return(points[chull(points$x, points$y), ])
}


compute_neff <- function(chain) {
  
  stopifnot(is.matrix(chain))
  nvars <- ncol(chain)
  N <- nrow(chain)
  ess <- rep(0, nvars)
  
  ess <- mclapply(1 : nvars, function(i) {
    if (all(chain[, i] == chain[1, i])) {
      ess[i] <- 0
    } else {
      acf_all <- acf(chain[, i], plot = FALSE, lag.max = N/2)$acf
      l <- length(acf_all)
      pair_adds <- acf_all[2:l] + acf_all[1:(l - 1)]
      if(any(pair_adds < 0)) stopping_pt <- min(which(pair_adds < 0)) else
        stopping_pt <- l
      acf_sum <- sum(acf_all[2 : max(stopping_pt - 1, 2)])
      ess[i] <- N / (1 + 2 * acf_sum)
    }
  }, mc.cores = min(detectCores() - 2L, 16L)) %>% unlist()
  
  ess
}


###############################################################################
###############################################################################
### NOTES ON INLA PARAMETERISATION
### THE PRIOR CONSTRUCTION IS IN INLA::param2.matern.orig
### RANGE PRIOR IS A FUNCTION OF MESH
### PRIOR.VARIANCE.NOMINAL IS ALWAYS 1 unless changed
### PRECISIONS OF PARAMETERS ARE ALWAYS 0.1 unless changed
### Below is sample code from the function:
# if (is.null(prior.range.nominal)) {
#   mesh.range = inla.ifelse(d == 2, (max(c(diff(range(mesh$loc[,1])), 
#                                           diff(range(mesh$loc[,2])), 
#                                           diff(range(mesh$loc[,3]))))), 
#                               diff(mesh$interval))
#   prior.range.nominal = mesh.range * 0.2
# }
# if (is.null(prior.kappa)) {
#   prior.kappa = sqrt(8 * nu.nominal)/prior.range.nominal
# }
# if (is.null(prior.tau)) {
#   prior.tau = sqrt(gamma(nu.nominal)/gamma(alpha.nominal)/
#                      (4 *pi * prior.kappa^(2 * nu.nominal) * 
#                         prior.variance.nominal))
# }
### NOTE ALSO WEIRD CONSTRUCTION. 
### B.kappa = (log(kappa0),0,1) => that
### theta0 = log(kappa0) (fixed), and theta2 weights an offset from log(kappa0)
### B.kappa = (log(tau0),1,0) => that
### theta0 = log(tau0) (fixed), and theta1 weights an offset from log(tau0)
### If left undefined, then then the mean if captured in the prior of theta
### Note that the equation is log(kappa) = log(kappa0) + theta2 where the prior
### is on theta2, similarly for tau
