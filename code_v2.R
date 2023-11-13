# Load required packages
require(pacman) # Package manager for R
p_load(sf, tidyverse, tmap, ggplot2, geoR, mvtnorm, msm, terra, raster, stars, PrevMap)

# Define the empirical logit transform function
elogit <- function(num, den) {
  log((num + 0.5) / (den - num + 0.5))
}

# Function to convert longitude and latitude to Cartesian coordinates in a specified CRS
longlat2xy <- function(long, lat, crs){
  data <- data.frame(long = long, lat= lat) %>% st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
    st_transform(., crs= crs) %>%
    st_coordinates() %>% data.frame # %>% rename(., utmx = x, utmy = y)
  return(data.frame(data))
}

# Function to create labels for plot axes
create_labels <- function(x, greater = F, smaller = F) {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = " - ")
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }
  
  return(labs)
}

# Function to convert EPSG code to kilometers
epsgKM <- function(x) {
  crs <- st_crs(x)
  proj4KM <- gsub(pattern = "+.units=m", replacement = "+units=km", 
                  crs$proj4string)
  return(proj4KM)
}

# Function to estimate parameters for the non-stationary model
estimate_parameters <- function(model_number, n, beta, scale_s1, scale_s2, scale_s3, sigma_s, sigma_z = NULL,
                                simulated_data, kappa1, kappa2, kappa3) {
  
  true_par <- c(beta, sigma_s,  scale_s1, scale_s2, scale_s3,  sigma_z)
  # Simulate data
  # simulated_data <- simulate_data_fun(model_number, n, beta, scale_s1, scale_s2, scale_s3, sigma_s, sigma_z)
  
  # Extract simulated data
  X <- simulated_data$X
  y <- simulated_data$y
  p  <- length(beta)
  # Define log-likelihood function
  log_likelihood <- function(params, X, y, model_number) {
    beta <- params[1:p]
    sigma_s <- exp(params[p+1])
    scale_s1 <- exp(params[p+2])
    scale_s2 <- exp(params[p+3])
    scale_s3 <- exp(params[p+4])
    if (!is.null(sigma_z)) {
      sigma_z <- exp(params[p+5])
    } else {
      sigma_z = NULL
    }
    
    # Compute spatial random effect at observed locations
    D <- as.matrix(dist(simulated_data$coords))
    D_e <- as.matrix(dist(simulated_data$e))
    D_t <- as.matrix(dist(simulated_data$t))
    R1 <- matern(u = D, phi = scale_s1, kappa = kappa1)
    R2 <- matern(u = D_e, phi = scale_s2, kappa = kappa2)
    R3 <- matern(u = D_t, phi = scale_s3, kappa = kappa3)
    
    
    # Simulate spatial random effects based on the selected model
    if (model_number == 1) {
      if (!is.null(sigma_z)) {
        Z <- rmvnorm(1, sigma = sigma_z * diag(n))
        W <- sigma_s * (R1 * R2 * R3) + sigma_z * diag(n)
      } else {
        Z <- matrix(0, nrow = n, ncol = 1)
        W <- sigma_s * (R1 * R2 * R3) 
      }
    } else if (model_number == 2) {
      if (!is.null(sigma_z)) {
        Z <- rmvnorm(1, sigma = sigma_z * diag(n))
        W <- sigma_s*((R1*R2) + (R1*R3)) + sigma_z * diag(n)
      } else {
        Z <- matrix(0, nrow = n, ncol = 1)
        W <- sigma_s*((R1*R2) + (R1*R3)) 
      }
    } else if (model_number == 3) {
      if (!is.null(sigma_z)) {
        Z <- rmvnorm(1, sigma = sigma_z * diag(n))
        W <- sigma_s*(R1+R2+R3) + sigma_z * diag(n)
      } else {
        Z <- matrix(0, nrow = n, ncol = 1)
        W <- sigma_s*(R1+R2+R3) 
      }
    } else {
      stop("Invalid model number. Choose model_number = 1, 2, or 3.")
    }
    
    mu <- as.numeric(X %*% beta) 
    
    # Compute log-likelihood
    loglik <- dmvnorm(x = y, mean = mu, sigma = W, log = TRUE) 
    return(loglik)
  }
  
  # Initial parameter values
  # init_params <- c(rnorm(p), log(runif(4)))
  init_params <- c(beta, log(runif(4)))
  if (!is.null(sigma_z)) {
    init_params <- c(init_params, log(runif(1)))
  }
  
  # Optimize log-likelihood function to estimate parameters
  fit <- optim(init_params, log_likelihood, X = X, y = y, model_number=model_number,  control = list(fnscale = -1))
  
    
    # Create a data frame to store the results
  results <- fit$par
  
  return(list(results = results, data = simulated_data))
}

#### prediction function

# Define a function named 'prediction_fun' with parameters 'par', 'model_number', 'sigma_z', 'pred_loc', and 'pred_predictor'
prediction_fun <- function(par, model_number, sigma_z, pred_loc = NULL, pred_predictor = NULL){
  
  # Load required packages
  library(proxy)
  library(Matrix)
  
  # Extract beta coefficients from the results
  res <- par
  beta <- res$results$Estimated[1:5]
  n <- length(res$simulated_data$y)
  
  # Generate prediction locations if not provided
  if(is.null(pred_loc)){
    n.pred <- 50
    pred.ind <- sample(x= n, size = n.pred)
    predictors <- res$simulated_data$X[pred.ind, ]
    mu.pred <- as.numeric(predictors %*% beta)
    grid.pred <- res$simulated_data$coords[pred.ind, ]
  } else {
    n.pred <- nrow(pred_loc)
    pred.ind <- 1:n.pred
    predictors <- pred_predictor
    mu.pred <- as.numeric(cbind(1, predictors) %*% beta)
    grid.pred <- pred_loc
  }
  
  mu <- as.numeric(res$simulated_data$X %*% beta)
  
  # Compute distance matrices for data
  D <- as.matrix(dist(res$simulated_data$coords))
  D_e <- as.matrix(dist(res$simulated_data$e))
  D_t <- as.matrix(dist(res$simulated_data$t))
  
  # Compute correlation matrices for data
  R1 <- matern(u = D, phi = res$results$Estimated[5], kappa = 0.5)
  R2 <- matern(u = D_e, phi = res$results$Estimated[6], kappa = 0.5)
  R3 <- matern(u = D_t, phi = res$results$Estimated[7], kappa = 0.5)
  
  # Compute distance matrices for predictions
  D.pred <- as.matrix(dist(grid.pred))
  D_e.pred <- as.matrix(dist(predictors[,1]))
  D_t.pred <- as.matrix(dist(predictors[,2]))
  
  # Compute correlation matrices for predictions
  R1.pred <- matern(u = D.pred, phi = res$results$Estimated[5], kappa = 0.5)
  R2.pred <- matern(u = D_e.pred, phi = res$results$Estimated[6], kappa = 0.5)
  R3.pred <- matern(u = D_t.pred, phi = res$results$Estimated[7], kappa = 0.5)
  
  # Compute distance matrices for cross-correlation
  D.cross <- as.matrix(pdist::pdist(grid.pred, res$simulated_data$coords))
  D_e.cross <- as.matrix(dist(predictors[,1], res$simulated_data$e))
  D_t.cross <- as.matrix(dist(predictors[,2], res$simulated_data$t))
  
  # Compute correlation matrices for cross-correlation
  R1.cross <- matern(u = D.cross, phi = res$results$Estimated[5], kappa = 0.5)
  R2.cross <- matern(u = D_e.cross, phi = res$results$Estimated[6], kappa = 0.5)
  R3.cross <- matern(u = D_t.cross, phi = res$results$Estimated[7], kappa = 0.5)
  
  # Simulate spatial random effects based on the selected model
  if (model_number == 1) {
    if (!is.null(sigma_z)) {
      Z <- rmvnorm(1, sigma = res$results$Estimated[8] * diag(n))
      W <- res$results$Estimated[4] * (R1 * R2 * R3) + res$results$Estimated[8] * diag(n)
      W.pred <- res$results$Estimated[4] * (R1.pred * R2.pred * R3.pred) + res$results$Estimated[8] * diag(n.pred)
      W.cross <- res$results$Estimated[4] * (R1.cross * R2.cross * R3.cross) + 
        res$results$Estimated[8] * diag(x = 1, nrow = n.pred, ncol = n)
    } else {
      Z <- matrix(0, nrow = n, ncol = 1)
      W <- res$results$Estimated[4] * (R1 * R2 * R3) 
      W.pred <- res$results$Estimated[4] * (R1.pred * R2.pred * R3.pred)
      W.cross <- res$results$Estimated[4] * (R1.cross * R2.cross * R3.cross)
    }
  } else if (model_number == 2) {
    if (!is.null(sigma_z)) {
      Z <- rmvnorm(1, sigma = res$results$Estimated[8] * diag(n))
      W <- res$results$Estimated[4]*((R1*R2) + (R1*R3)) + res$results$Estimated[8] * diag(n)
      W.pred <- res$results$Estimated[4]*((R1.pred*R2.pred) + (R1.pred*R3.pred)) + res$results$Estimated[8] * diag(n.pred)
      W.cross <- res$results$Estimated[4]*((R1.cross*R2.cross) + (R1.cross*R3.cross)) + 
        res$results$Estimated[8] * diag(x = 1, nrow = n.pred, ncol = n)
    } else {
      Z <- matrix(0, nrow = n, ncol = 1)
      W <- res$results$Estimated[4]*((R1*R2) + (R1*R3)) 
      W.pred <- res$results$Estimated[4]*((R1.pred*R2.pred) + (R1.pred*R3.pred)) 
      W.cross <- res$results$Estimated[4]*((R1.cross*R2.cross) + (R1.cross*R3.cross))
    }
  } else if (model_number == 3) {
    if (!is.null(sigma_z)) {
      Z <- rmvnorm(1, sigma = res$results$Estimated[8] * diag(n))
      W <- res$results$Estimated[4]*(R1+R2+R3) + res$results$Estimated[8] * diag(n)
      W.pred <- res$results$Estimated[4]*(R1.pred+R2.pred+R3.pred) + res$results$Estimated[8] * diag(n.pred)
      W.cross <- res$results$Estimated[4]*(R1.cross+R2.cross+R3.cross) + 
        res$results$Estimated[8] * diag(x = 1, nrow = n.pred, ncol = n)
    } else {
      Z <- matrix(0, nrow = n, ncol = 1)
      W <- res$results$Estimated[4]*(R1+R2+R3) 
      W.pred <- res$results$Estimated[4]*(R1.pred+R2.pred+R3.pred) 
      W.cross <- res$results$Estimated[4]*(R1.cross+R2.cross+R3.cross)
    }
  } else {
    stop("Invalid model number. Choose model_number = 1, 2, or 3.")
  }
  
  # Compute inverse of Sigma.W matrix
  Sigma.W.inv <- solve(W)
  
  # Compute A matrix for conditional mean
  A <- W.cross %*% Sigma.W.inv
  
  # Compute conditional mean and variance
  mu.cond <- as.numeric(mu.pred + A %*% (res$simulated_data$y -  mu))
  sigma.cond <- W.pred - A %*%t(W.cross)
  
  # Simulate 100 spatial predictions
  Y.hat.all <- try(rmvnorm(100, mean = mu.cond, sigma = sigma.cond))
  
  # Handle errors and use nearPD if needed
  if(inherits(Y.hat.all, "try-error")){
    Y.hat.all <- try(rmvnorm(100, mean = mu.cond, sigma = nearPD(sigma.cond)$mat))
    if (inherits(Y.hat.all, "try-error")){
      sigma.cond <- (sigma.cond + t(sigma.cond))/2
      Y.hat.all <- rmvnorm(100, mean = mu.cond, sigma = nearPD(sigma.cond)$mat)
    }
  }
  
  # Compute mean, lower, and upper quantiles of the simulated predictions
  Y.hat <- apply(Y.hat.all, 2, mean)
  Y.hat.lower <- apply(Y.hat.all, 2, quantile, prob=c(0.025))
  Y.hat.upper <- apply(Y.hat.all, 2, quantile, prob=c(0.975))
  
  # Return a list of results
  to_return <- list(Y=res$simulated_data$y, Y.hat = Y.hat, Y.hat.lower = Y.hat.lower, 
                    Y.hat.upper = Y.hat.upper)
  return(to_return)
}



############################### Analysis #######################################

# Load malaria data from the specified RDS file
data <- readRDS("../../data/malaria_data_ethiopia_edited.rds")

############ fit a non-spatial model ##########

# Fit a non-spatial linear model (lm) using logit-transformed response and predictor variables
fit0 <- lm(logitp ~ alt + temp  + hum + dist_aqua, data = data)

# Display summary statistics of the fitted model
summary(fit0)

# Perform stepwise variable selection on the fitted model
summary(step(fit0))

######## Create data for the analysis #########

# Create a list containing design matrices, response, and coordinates for spatial analysis
data_dummy <- list(X = model.matrix(fit0), y = model.response(model.frame(fit0)), 
                   coords = as.matrix(data[, c("x", "y")]), e = data$alt, t = data$temp)


# Estimate parameters for the spatial model (model_number = 2) using the new non-spatial data
fit1 <- estimate_parameters(model_number = 2, n = nrow(data), beta = coef(fit0), 
                            scale_s1 = 1, scale_s2 = 1, scale_s3 = 1, 
                            sigma_s = 1, sigma_z = 0.05, simulated_data = data_dummy,
                            kappa1= 1.5, kappa2 = 2.5, kappa3 = 1.5)

# Display results of the new fitted spatial model
fit1$results

######### prediction ######

# Read prediction data from the specified CSV file
pred_data <- read_csv("../../data/dp.csv")

# Extract coordinates and predictors for prediction
pred_grid <- pred_data[, 1:2]  %>% st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = epsgKM(2736)) %>%
  st_coordinates()
predictors <- as.matrix(pred_data[, c(3, 4, 6, 7)])

# Generate predictions using the previously fitted spatial model (model_number = 2)
predictions <- prediction_fun(par = fit1, model_number = 2, sigma_z = 0.05, 
                              pred_loc = pred_grid, pred_predictor = predictors)



# Load required libraries
require(tmap)

# Read and transform shapefile data
shp <- st_read("../../data/gadm41_MOZ_shp/gadm41_MOZ_0.shp") %>%
  st_transform(., crs = epsgKM(2736))

# Define breaks, labels, and color palette for the map
brks <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.91, 0.92, 0.95, 1)
labs <- create_labels(brks)
pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)

# Create the first map (m1) with observed data
m1 <- data %>% st_as_sf(., coords = c("x", "y"), crs = epsgKM(2736)) %>%
  tm_shape() + tm_symbols(col = "prev", size =0.1, palette = pal, 
                          breaks = brks, style = "fixed", legend.col.show = F) + 
  tm_shape(shp, is.master = TRUE) + tm_borders(lwd = 1.5, col="black") +
  tm_compass(position = c("right", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 1, legend.title.fontface = 2.35, 
            legend.position = c("left","top"),  panel.labels = "Empirical prevalence", panel.label.size=1)

# Display the first map (m1)
m1

# Create stars object and rasterize predicted data
r0 <- st_as_stars(st_bbox(shp), nx = 130, ny = 120)
p <- st_as_sf(data.frame(cbind(pred_grid, prev = plogis(predictions$Y.hat), lprev = plogis(predictions$Y.hat.lower),
                               uprev = plogis(predictions$Y.hat.upper))), coords = c("X", "Y"), crs = epsgKM(2736))
ras_shp <- st_rasterize(p, r0)

# Crop the raster to the shape of the country
ras_shp <- st_crop(ras_shp, shp)

# Create three additional maps (m2, m3, m4) for predicted mean and confidence intervals
m2 <- tm_shape(ras_shp[1]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Predicted mean \nprevalence (%)")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"), panel.labels = "Predicted Mean",  panel.label.size=1)

m3 <- tm_shape(ras_shp[2]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Lower limit of \n95% CI of prevalence")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"),  panel.labels = "Lower 95%", panel.label.size=1)

m4 <- tm_shape(ras_shp[3]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Upper limit of \n95% CI of prevalence")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"), panel.labels = "Upper 95%", panel.label.size=1)

# Arrange the maps in a grid layout
tmap_arrange(m1, m2, m3, m4)


################## fit a stationary model to the data ########################

######### Linear model ################

# Fit a linear model using Maximum Likelihood Estimation (MLE)
fit <- linear.model.MLE(formula = logitp ~ alt + temp  + hum + dist_aqua, 
                        coords = ~ x + y, 
                        data = data,
                        kappa = 1.5, fixed.rel.nugget = NULL, 
                        start.cov.pars = c(30, 1)) 

# Display summary of the linear model
summary(fit, log = F)

#### calculate AIC and BIC

# Calculate AIC and BIC for model comparison
AIC <- 2 * (-(fit$log.lik-(nrow(data)/2)*log(2*pi))) + 2 * length(fit$estimate)
BIC <- 2 * (-(fit$log.lik-(nrow(data)/2)*log(2*pi))) + length(fit$estimate) * log(nrow(data))
AIC
BIC

######## confidence interval

# Get the parameter estimate 
require(msm)
ress <- fit$estimate
ress[6:8] <- exp(ress[6:8])
mn <- ress
se <- deltamethod(g = list(~x1, ~x2, ~x3, ~x4, ~x5, ~exp(x6), ~exp(x7), ~exp(x8)), mean = fit$estimate, cov = fit$covariance)

# Function to calculate confidence intervals
myconinter <- function(mn, se, ci=c(0.025, 0.975)){
  a <- matrix(rbind(mn), nrow = 8, ncol = 2) + cbind(se) %*%  rbind(qnorm(ci))
  b <- cbind(mn, a)
  round(b, digits = 4)
}

# Display confidence intervals
myconinter(mn, se)

# Make predictions using the linear model
predictions <- spatial.pred.linear.MLE(object = fit, 
                                       grid.pred = pred_grid, 
                                       predictors = data.frame(predictors), 
                                       scale.predictions = c("logit", "prevalence", "odds"), 
                                       quantiles = c(0.025, 0.975), 
                                       standard.errors = T, n.sim.prev = 100)

# Create stars object and rasterize predicted data
r0 <- st_as_stars(st_bbox(shp), nx = 130, ny = 120) 
p <- st_as_sf(data.frame(cbind(pred_grid, prev = plogis(predictions$logit$predictions), 
                               lprev = plogis(predictions$logit$quantiles[,1]),
                               uprev = plogis(predictions$logit$quantiles[,2]))), coords = c("X", "Y"), crs = epsgKM(2736))
ras_shp <- st_rasterize(p, r0)

# Crop the raster to the shape of the country
pred_mean_sd <- st_crop(ras_shp, shp)

# Plot the results
plot(pred_mean_sd)

# Set up color palette for mapping
brks <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.91, 0.92, 0.95, 1)
labs <- create_labels(brks)
pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)

# Create and display map for predicted mean prevalence
tm_shape(pred_mean_sd) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = "Predicted mean\nlead concentration\n(log-scale)") +
  tm_shape(shp) +
  tm_borders(col = "black") +
  tm_compass() +
  tm_scale_bar(position = c("left", "bottom")) 

# Create three additional maps for lower and upper limits of the confidence interval
m5 <- tm_shape(pred_mean_sd[1]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Predicted mean \nprevalence (%)")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"), panel.labels = "Predicted Mean",  panel.label.size=1)

m6 <- tm_shape(pred_mean_sd[2]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Lower limit of \n95% CI of prevalence")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"),  panel.labels = "Lower 95%", panel.label.size=1)

m7 <- tm_shape(pred_mean_sd[3]) +
  tm_raster(palette = pal, style = "fixed", legend.show = F, breaks = brks, 
            title = c("Upper limit of \n95% CI of prevalence")) +
  tm_shape(shp, is.master=T) +
  tm_borders("black", lwd = 1.5) +
  tm_compass(position = c("LEFT", "TOP")) +
  tm_scale_bar(position = c("right", "BOTTOM")) +
  tm_layout(frame.lwd = 3, 
            legend.title.size = 0.95, legend.title.fontface = 1.35, 
            legend.position = c("left","top"), panel.labels = "Upper 95%", panel.label.size=1)

# Arrange the maps in a grid layout
tmap_arrange(m2, m3, m4, m5, m6, m7)
