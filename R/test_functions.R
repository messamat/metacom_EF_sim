source("R/packages.R")
source("R/sim_setup.R")
source("R/MC_simulate.R")

rootdir = rprojroot::find_root(rprojroot::has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

timesteps = 100
out_SSNdir = file.path(resdir,
                       paste0('simSSN_', 
                              format(Sys.Date(), "%Y%m%d")))

landscape <- generate_OCNigraph(patches=100, 
                                cellsize = 0.5, 
                                dimX = 25, 
                                dimY = 25, 
                                plot = T,
                                out_format = 'SSN',
                                out_SSNdir = out_SSNdir)


if (spatial_autocor) {
  ##############################################################################
  #To compute random fields over stream distances, tried generating synthetic
  #x and y coordinates with distances in euclidean space (because
  #providing "distances" to RFSimulate didn't work) but that's really hard and 
  #no options yield exact results due to non-euclidean space
  dist_mat <- compute_distmat(landscape = landscape)
}

if (inherits(landscape, 'igraph')) {
  patches <- gorder(landscape)
  landscape <- as.data.table(get.vertex.attribute(landscape))
  
} else if (inherits(landscape, 'data.frame')) {
  patches <- nrow(landscape)
}

Sys.setenv(USE_CXX14 = 1)

library('Rcpp')
library('ggplot2')
library('dplyr')
library('RColorBrewer')
library('rstan')
library('bayesplot')
library('nlme')
library('SSN')
library('SSNbayes')
library('viridis')

## Set some useful options for modelling
rstan_options(auto_write = TRUE) # avoid recompilation
options(mc.cores = parallel::detectCores())
RNGkind(sample.kind = "Rounding")

## Set the seed for reproducibility
seed <- 202008
set.seed(seed)

ssn <- landscape

plot(ssn, lwdLineCol = "addfunccol",  lwdLineEx = 8,
     lineCol = 4,  col = 1,  pch = 16,  xlab = "x-coordinate",  ylab = "y-coordinate")

## Create stream distance matrices
SSN::createDistMat(ssn, o.write=TRUE)


## Extract the data.frames for the observed and prediction location data
rawDFobs <- SSN::getSSNdata.frame(ssn, Name = "Obs")

## Extract the geographic coordinates from the SpatialStreamNetwork
## object and add to data.frames
obs_data_coord <- data.frame(ssn@obspoints@SSNPoints[[1]]@point.coords)
obs_data_coord$pid<- as.numeric(rownames(obs_data_coord))
rawDFobs<- rawDFobs %>% left_join(obs_data_coord, by = c("pid"),
                                  keep = FALSE)
rawDFobs$point <- "Obs" ## Create label for observed points


## Generate 3 continous covariates at observed and prediction locations
set.seed(seed)
rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X3"] <- rnorm(length(rawDFobs[,1]))

## Ensure the rownames still match the pid values used in the
## SpatialStreamNetwork object
rownames(rawDFobs)<- as.character(rawDFobs$pid)

## Generate 3 continous covariates at observed and prediction locations
set.seed(seed)
rawDFpred[,"X1"] <- rnorm(length(rawDFpred[,1]))
rawDFpred[,"X2"] <- rnorm(length(rawDFpred[,1]))
rawDFpred[,"X3"] <- rnorm(length(rawDFpred[,1]))

rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFobs[,"X3"] <- rnorm(length(rawDFobs[,1]))

## Ensure the rownames still match the pid values used in the
## SpatialStreamNetwork object
rownames(rawDFobs)<- as.character(rawDFobs$pid)
rownames(rawDFpred)<- as.character(rawDFpred$pid)

## Put the new covariates back in the SpatialStreamNetwork object
ssn <- SSN::putSSNdata.frame(rawDFobs,ssn, Name = 'Obs')

## Simulate the response variable at observed and prediction locations
set.seed(seed)
sim.out <- SSN::SimulateOnSSN(ssn.object = ssn,
                              ObsSimDF = rawDFobs, ## observed data.frame
                              #PredSimDF = rawDFpred, ## prediction data.frame
                              #PredID = "preds", ## name of prediction dataset
                              formula = ~ X1 + X2 + X3,
                              coefficients = c(10, 1, 0, -1), ## regression coefficients
                              CorModels = c("Exponential.taildown"), ## covariance model
                              use.nugget = TRUE, ## include nugget effect
                              CorParms = c(3, 10, .1)) ## covariance parameters

## Extract the SpatialStreamNetwork object from the list returned by
## SimulateOnSSN and extract the observed and prediction site
## data.frames. Notice the new column Sim_Values in the data.frames
sim.ssn <- sim.out$ssn.object
simDFobs <- SSN::getSSNdata.frame(sim.ssn,"Obs")

## Create a data.frame containing training and test data.
df_obs <- SSN::getSSNdata.frame(sim.ssn, "Obs") ## Extract observed dataset
df_obs$dataset <- 'train' ## Create new column 'dataset' and set to 'train'

## Expand data.frames to include t=timesteps days per location
df_obs <- do.call("rbind", replicate(timesteps, df_obs, simplify = FALSE))# replicating the df
df_obs$timestep <- rep(1:timesteps, each = (nrow(df_obs)/timesteps)) # Set date variable

## Create a copy of the pid value used in the SpatialStreamNetwork
## object and create a new pid value for use in SSNbayes
## package. Values must be consequtively ordered from 1 to the number
## of rows in the data.frame
df_obs <- df_obs %>% mutate(pid.ssn = pid,
                            pid = rep(1:nrow(.)))

set.seed(seed)
phi <- 0.8 ## lag 1 autocorrelation value
ar1_sim <- nlme::corAR1(form = ~ timestep, value = phi) # can also use corExp function
AR1 <- nlme::Initialize(ar1_sim, 
                        data = data.frame(timestep = unique(df_obs$timestep)))

## Create a vector of AR1 errors for each date and expand to all all locations
#NB AR1 error
epsilon <- t(chol(corMatrix(AR1))) %*% rnorm(length(unique(df_obs$timestep)), 0, 3) 
epsilon <- rep(epsilon, each = length(unique(df_obs$locID)) ) +
  rnorm(length(epsilon)*length(unique(df_obs$locID)), 0, 0.25) # for all the locations
epsilon_df <- data.frame(timestep = rep(unique(df_obs$timestep), each = length(unique(df_obs$locID))),
                         locID = rep(unique(df_obs$locID), times = length(unique(df_obs$timestep))),
                         epsilon = epsilon)
df_obs <- df_obs %>% left_join(epsilon_df, by = c('timestep' = 'timestep', 'locID' = 'locID'))

## Create a new simulated response variable, y, with errors added
df_obs$y <- df_obs$Sim_Values + df_obs$epsilon 

## Create line plots of the response over time for training and test datasets
ggplot(df_obs) +
  geom_line(aes(x = timestep, y = y, group = locID), alpha = 0.4) +
  ylab("Simulated Temperature (\u00B0C)")+
  facet_wrap(~dataset)+
  theme_bw()

## Extract stream (edge) network structure, including the additive function value
nets <- SSNbayes::collapse(ssn, par = 'addfunccol')


## Create additive function value categories for plotting
nets$afv_cat <- cut(nets$addfunccol,
                    breaks = seq(min(nets$addfunccol),
                                 max(nets$addfunccol),
                                 length.out=6),
                    labels = 1:5,
                    include.lowest = T)




## Plot simulated temperature, by date, with line width proportional to afv_cat
ggplot(nets) +
  geom_path(aes(X1, X2, group = slot, size = afv_cat), lineend = 'round',
            linejoin = 'round', col = 'lightblue')+
  geom_point(data = dplyr::filter(df_obs, timestep %in% 1:10),
             aes(x = coords.x1, y = coords.x2, col = y, shape = point),
             size = 1)+
  scale_size_manual(values = seq(0.2,2,length.out = 5))+
  facet_wrap(~timestep, nrow = 2)+
  scale_color_viridis(option = 'C')+
  scale_shape_manual(values = c(16,15))+
  xlab("x-coordinate") +
  ylab("y-coordinate")+
  theme_bw()





