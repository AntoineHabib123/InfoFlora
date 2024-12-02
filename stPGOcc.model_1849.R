#This script create the mhttp://127.0.0.1:18963/graphics/plot_zoom_png?width=1443&height=1018odel on all sites other than those indicated as na in wsprop (sites that are adjacent to WS with no observations but have atlas information)

library(pacman)
p_load(httr, httpuv,brms, StanHeaders,  spOccupancy, jsonlite, dplyr, data.table, readr, sf, ggplot2, sp, terra, mapview, MCMCvis, raster, viridisLite, viridis, tidyr, tibble, tidyverse, boot)
set.seed(500)

####building the model specific data####
wsprop_na_indicies <- which(is.na(taxa_all_info$wsprop_model))

#subseting again to create data only for model run
taxa_all_info_wsmodel <- subset(taxa_all_info, !is.na(wsprop_model))

y_matrix <- unlist(st_drop_geometry(taxa_all_info_wsmodel[35:64]))

y <- array(y_matrix,dim=c(nrow(taxa_all_info_wsmodel),30,1))

years <- matrix(rep(c(1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023),each = nrow(taxa_all_info_wsmodel)),nrow = nrow(taxa_all_info_wsmodel))

years <- matrix(rep(c(1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023),each = nrow(taxa_all_info_wsmodel)),nrow = nrow(taxa_all_info_wsmodel))

Occcov <- list(as.vector(array(st_drop_geometry(taxa_all_info_wsmodel[[1]]))),as.vector(array(st_drop_geometry(taxa_all_info_wsmodel[[3]]))),as.vector(array(st_drop_geometry(taxa_all_info_wsmodel[[4]]))),years)
names(Occcov) <- c("id_grid","wsprop","suitability", "years")

effort <- unlist((st_drop_geometry(taxa_all_info_wsmodel[5:34])))
effort <- array(effort,dim=c(nrow(taxa_all_info_wsmodel),30))

Detcov <- list(effort,st_drop_geometry(taxa_all_info_wsmodel[[1]]), years) 
names(Detcov) <- c("effort","id_grid","years")

coords <- st_coordinates(st_centroid(taxa_all_info_wsmodel))

data_list <- list(
  y = y,
  occ.covs = Occcov,
  det.covs = Detcov,
  coords = coords
)

####want to pass into model the relevant information####
n.sites = nrow(taxa_all_info)
z.inits <- y[,,1]
# Starting data.frame()# Starting values for the MCMC chains
inits <- list(
  beta = 0, #occurence coefficient
  alpha = 0,  #detection coefficient
  sigma.sq.psi = 1, # occurence random offect variance
  z = z.inits
)
#alpha beta priors calculated specifically, see markdown (7.priors) for more info
priors <- list( 
  beta.normal = list(mean = 0, var = 100),   # Vague prior for occupancy
  alpha.normal = list(mean = 0, var = 100),  # Vague prior for detection
  sigma.sq.psi.ig = list(a = 0.01, b = 0.01))

n.chains <- 5
n.batch <- 400
batch.length <- 50
n.samples <- n.batch * batch.length 
n.burn <- 2000
n.thin <- 12
ar1 <- TRUE

occ_formula <- ~ scale(suitability)+I(scale(suitability)^2)+scale(wsprop)+I(scale(wsprop)^2)
#occ_formula <- ~ scale(wsprop)+I(scale(wsprop)^2)

det_formula <- ~scale(effort)
####the model + functions to save####
out <- stPGOcc(occ.formula = occ_formula, #environ 13 mins taxa dependent
               det.formula = det_formula, 
               data = data_list, 
               n.batch = n.batch, 
               batch.length = batch.length,
               inits = inits,
               priors = priors,
               ar1 = ar1,
               n.burn = n.burn, 
               n.thin = n.thin, 
               n.chains = n.chains, 
               n.report = 50,
               n.omp.threads = 11)
summary(out)
saveRDS(out, file=paste0(your.path,"DynamicModels30/stPGOcc_5ts_6chains_fes.",taxon_names[l],".Model.Rdata"))
####Predictions####

# Number of prediction sites.
J.pred <- n.sites
# Number on.sites# Number of prediction years.
n.years.pred <- 30
# Number of predictors (including intercept)
p.occ <- ncol(out$beta.samples)
# Get covariates and standardize them using values used to fit the model
predict.data <- taxa_all_info
suitability.pred <- (predict.data$suitability - mean(data_list$occ.covs$suitability)) / sd(data_list$occ.covs$suitability)
year.pred <- matrix(rep((c(1994, 2023) - mean(data_list$occ.covs$years)) / 
                          sd(data_list$occ.covs$years), 
                        length(suitability.pred)), J.pred, n.years.pred, byrow = TRUE) # a priori pas besoin pour notre modèle
wsprop.pred <- (predict.data$wsprop_pred - mean(data_list$occ.covs$wsprop)) / sd(data_list$occ.covs$wsprop)
# Create three-dimensional array
X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ))
# Fill in the array
# suitability
X.0[, , 2] <- suitability.pred
# suitability^2
X.0[, , 3] <- suitability.pred^2
# wsprop
X.0[, , 4] <- wsprop.pred
# wsprop
X.0[, , 5] <- wsprop.pred^2
# Check out the structure
str(X.0)

# Indicate which primary time periods (years) we are predicting for
t.cols <- c(1:30)


coords.full <- st_coordinates(st_centroid(taxa_all_info))

# Approx. run time: < 30 sec
out.pred <- predict(out, X.0, t.cols = t.cols, coords.full, ignore.RE = TRUE, type = 'occupancy')
# Check out the structure
str(out.pred)
saveRDS(out.pred, file=paste0(your.path,"DynamicModels30/stPGOcc_",taxon_names[l],"2_outpred.Rdata"))

####view predictions at time n####
n=30 #n=1 is 1994, n=30 is 2023
mean.psi.n <- apply(out.pred$psi.0.samples[,,n],2,mean)
taxa_all_info$psi.pred_time.n <- mean.psi.n
mapview(taxa_all_info, zcol="psi.pred_time.n")

####graphing z.0.samples (pred_presence) with actual (DB)####
#sum to find predicted presences (z-samples)
total_presence <- vector("numeric",length=30)
for(i in 1:length(total_presence)){
  total_presence[i]<- sum(st_drop_geometry(taxa_all_info[,i+34]))
}
z_0_sum_all_times <- apply(out.pred$psi.0.samples, c(1,3), sum) #summing to get #presence
mean_each_time <- apply(z_0_sum_all_times, 2, mean) #mean across all mcmc samples
sum_zsampels <-  data.frame(Time =years[1,], sum_zsamples = (mean_each_time), presence=total_presence)

#Presence and Zsamples
taxa_zplot <- ggplot(sum_zsampels[2:30,], aes(x = Time)) +
  geom_line(aes(y = sum_zsamples, color = "Zsample sum")) +
  geom_line(aes(y = presence, color = "Presence")) +
  labs(x = "Time", y = "Values", title = paste0("Line Graph of Zsamples & Presence over Time (Full Model)","Saxifraga")) +
  scale_color_manual(values = c("Zsample sum" = "blue", "Presence" = "red"),
                     name = "Legend") +
  theme_minimal()

saveRDS(taxa_zplot, file=paste0(your.path,"DynamicModels30/zsample_graph_2x2",taxon_names[l],".Rdata"))
#note:we don't think it is of value to include the first index of sum_zsamples since it is data based on years 1800-1994, unlike all subsequent time frames of 1 year. 


####caculating difference of presences for each time
sum(abs(sum_zsampels$sum_zsamples-sum_zsampels$presence))

readRDS(paste0(your.path,"DynamicModels30/zsample_graph_","festuca",".Rdata"))
