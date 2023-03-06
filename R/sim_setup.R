### ~~~~~~~~~~~~~~~~~~~~~~~~~~ convert_OCN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#' Converts an OCN object to an igraph object.
#'
#' This helper function converts an optimal channel network (OCN) created with
#' the R package OCNet to an igraph object, ready to use with the Cantal model.
#'
#' @param OCN an optimal channel network (OCN) created with the R package OCNet.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph set_vertex_attr
#' 
#' #' @author Adapted from Claire Jacquet
#' 
#' @return a igraph object.
#'
#' @examples
#' library(OCNet)
#' set.seed(1)
#' dimX = 15
#' dimY = 15
#' cellsize = 500
#' thrA = 5*cellsize^2
#' OCN <- create_OCN(dimX, dimY, cellsize = cellsize)
#' OCN <- landscape_OCN(OCN)
#' OCN <- aggregate_OCN(OCN, thrA = thrA)
#' graph <- OCN2graph(OCN)
#'
#' @export
#' 
convert_OCN <- function(OCN, out_format, out_SSNdir, idcol = 'patch') {
  out_OCN_formatted <- list()
  
  #Create graph
  #(besides outlets which do not have downstream edges)
  if ('igraph' %in% out_format) {
    igraph <- data.table(from = 1:OCN$RN$nNodes,
                         to = OCN$RN$downNode,
                         weight = OCN$RN$leng,
                         DA = OCN$RN$A) %>%
      .[to != 0,] %>% #remove outlet
      .[order(to),] %>%
      graph_from_data_frame %>% 
      set.vertex.attribute(name ='x', 
                           value = OCN$RN$X) %>%
      set.vertex.attribute(name ='y', 
                           value = OCN$RN$Y) %>%
      set.vertex.attribute(name = idcol,
                           value = 1:gorder(.))
    
    out_OCN_formatted <- c(out_OCN_formatted, igraph=list(igraph))
    
  } 
  
  if ("SSN" %in% out_format) {
    SSN <- OCN_to_SSN(OCN = OCN,
                      level = "RN",
                      obsSites = 1:OCN$RN$nNodes,
                      #predDesign, 
                      #predSites,
                      path = out_SSNdir,
                      randomAllocation = FALSE,
                      importToR = TRUE)
    
    ssn_df <- SSN::getSSNdata.frame(SSN)
    ssn_df[, idcol] <- ssn_df[, "locID"]
    SSN <- SSN::putSSNdata.frame(ssn_df, SSN)
    
    out_OCN_formatted <- c(out_OCN_formatted, SSN=SSN)
  }
  return(out_OCN_formatted)
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ generate_OCN_formatted ~~~~~~~~~~~~~~~~~~~~~~~~####

#' Generate OCN
#' @author inspiration and parameters from Claire Jacquet (OID: 10.1111/OIK.09372)
#' 
#' 
generate_OCN_formatted <- function(patches, out_format, out_SSNdir,
                                   cellsize = 0.5, dimX = 25, dimY = 25,
                                   outletPos = 3, expEnergy = 0.1, 
                                   coolingRate = 0.3, slope0 = 0.05,
                                   plot=TRUE) {
  set.seed(1)
  OCN <- create_OCN(dimX, dimY, 
                    cellsize = cellsize, 
                    outletPos = outletPos, 
                    expEnergy = expEnergy, 
                    coolingRate = coolingRate)
  OCN <- landscape_OCN(OCN, slope0 = slope0)
  
  #Determine the threshold area that generates the number of patches closest to "patches" argument
  thrA <- as.data.table(
    OCNet::find_area_threshold_OCN(OCN=OCN, thrValues=seq(0.5,50,0.05))) %>%
    .[which.min(abs(nNodesRN-patches)), thrValues]
  
  #thrA = 5*cellsize^2
  OCN <- aggregate_OCN(OCN, thrA = thrA)
  
  if (plot) {
    draw_thematic_OCN(rep(1, OCN$RN$nNodes),
                      OCN,
                      drawNodes=T,
                      addLegend=F,
                      cex=1,
                      backgroundColor = NULL)
  }
  
  graph <- convert_OCN(OCN = OCN, 
                       out_format = out_format,
                       out_SSNdir = out_SSNdir)
  
  return(graph)
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ landscape_generate ~~~~~~~~~~~~~~~~~~~~~~~~~####
#' Generate landscape
#'
#' Generates a landscape for metacommunity simulations
#'
#' @param patches number of patches to include in landscape
#' @param xy optional dataframe with x and y columns for patch coordinates
#' @param plot option to show plot of landscape
#'
#' @return landscape with x and y coordinates
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' landscape_generate()
#'
#' @export
#'
landscape_generate <- function(patches = 100, xy, plot = TRUE) {
  if (missing(xy)){
    
    if(patches > 10000) stop("Maximum number of patches is 10000.")
    #Generate random coordinates
    positions_linear = sample(0:9999, patches) #Sample without replacement
    landscape = data.frame(x = floor(positions_linear / 100)+1, 
                           y = (positions_linear %% 100) + 1)
    
    #Cluster landscape based on euclidean distance among sites
    clusters <- hclust(dist(landscape),method = "ward.D2")
    
    #Re-order patches based on nearest neighbors
    landscape <- landscape[clusters$order, ]
    #re-assign rowname ID
    rownames(landscape) <- 1:patches
    
  } else {
    landscape <- xy
  }
  if (plot == TRUE){
    plot(landscape, pch = 19)
  }
  return (landscape)
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ TODO: check dispersal matrix ~~~~~~~~~~~~~~~~####
#Check dispersal matrix
# disp_mat <- disp_mat
# rownames(disp_mat) <- 1:nrow(disp_mat)
# colnames(disp_mat) <- 1:ncol(disp_mat)
# if (is.matrix(disp_mat) == FALSE) stop ("disp_mat is not a matrix")
# if (nrow(disp_mat) != nrow(landscape) | ncol(disp_mat) != nrow(landscape)) 
#   stop ("disp_mat does not have a row and column for each patch in landscape")


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ TOCOMMENT: plot dispersal matrix ~~~~~~~~~~~~####
plot_dispersal_matrix <- function(disp_mat, print=TRUE) {
  g <- as.data.frame(disp_mat) %>%
    dplyr::mutate(to.patch = rownames(disp_mat)) %>%
    tidyr::gather(key = from.patch, value = dispersal, -to.patch) %>%
    dplyr::mutate(from.patch = as.numeric(as.character(from.patch)),
                  to.patch = as.numeric(as.character(to.patch))) %>%
    ggplot2::ggplot(ggplot2::aes(x = from.patch, y = to.patch, fill = dispersal))+
    ggplot2::geom_tile()+
    scale_fill_viridis_c()
  
  if (print) {print(g)}
  
  return(g)
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ compute_distmat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#### 
compute_distmat <- function(landscape, torus, idcol='patch') {
  if (inherits(landscape, 'igraph')) {
    dist_mat <- igraph::distances(landscape)
    colnames(dist_mat) <- get.vertex.attribute(landscape, 'patch')
    rownames(dist_mat) <- get.vertex.attribute(landscape, 'patch')
  } else {
    #Compute distance among sites
    if(torus == TRUE){
      dist_mat <- as.matrix(som.nn::dist.torus(coors = landscape))
    } else {
      dist_mat <- as.matrix(dist(landscape))
    }
  }
  
  return(dist_mat)
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ dispersal_matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#' Generate Dispersal Matrix
#'
#' Generates dispersal matrix for metacommunity simulations
#'
#' @param landscape landscape generated by landscape_generate()
#' @param torus whether to model the landscape as a torus
#' @param disp_mat optional matrix with each column specifying the probability that an individual disperses to each other patch (row)
#' @param kernel_exp the exponential rate at which dispersal decreases as a function of the distance between patches
#' @param plot option to show plot of environmental variation
#'
#' @return matrix with dispersal probabilities
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' dispersal_matrix(landscape_generate())
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom som.nn dist.torus
#'
#' @export
#'
#'
dispersal_matrix <- function(landscape, torus = TRUE, 
                             kernel_exp = 0.1, plot = TRUE){
  
  dist_mat <- compute_distmat(landscape = landscape,
                              torus = torus)
  
  disp_mat <- exp(-kernel_exp * dist_mat) #Exponential decrease in dispersal with distance (see equation 4; 0.1 is Li)
  diag(disp_mat) <- 0
  disp_mat <- apply(disp_mat, 1, function(x) x / sum(x)) #Standard into relative probability
  
  
  if (sum(colSums(disp_mat) > 1.001) > 0)
    warning ("dispersal from a patch to all others exceeds 100%. 
  Make sure the rowSums(disp_mat) <= 1")
  if (sum(colSums(disp_mat) < 0.999) > 0) 
    warning ("dispersal from a patch to all others is less than 100%. 
           Some dispersing individuals will be lost from the metacommunity")
  
  if (plot == TRUE){
    plot_dispersal_matrix(disp_mat) 
  }
  
  return (disp_mat)
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ plot_SSNenv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
plot_SSNenv <- function(snn, df_sim, 
                        timesteps_to_plot = 1:10, tscol = 'timestep') {
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
  p <- ggplot(nets) +
    geom_path(aes(X1, X2, group = slot, size = afv_cat), lineend = 'round',
              linejoin = 'round', col = 'lightblue')+
    geom_point(data = dplyr::filter(df_sim, 
                                    get(tscol) %in% timesteps_to_plot),
               aes(x = coords.x1, y = coords.x2, col = y, shape = point),
               size = 1)+
    scale_size_manual(values = seq(0.2,2,length.out = 5))+
    facet_wrap(~get(tscol), nrow = 2)+
    scale_color_viridis(option = 'C')+
    scale_shape_manual(values = c(16,15))+
    xlab("x-coordinate") +
    ylab("y-coordinate")+
    theme_bw()
  
  print(p)
  return(p)
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ simulate_envSSN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
simulate_envSSN <- function(ssn, timesteps, plot = FALSE, seed) {
  #Code from https://www.kaggle.com/code/edsans/ssnbayes-simulated/notebook
  #Tutorial with SSNBayes package
  
  ## Set some useful options for modelling
  RNGkind(sample.kind = "Rounding")
  
  ## Set the seed for reproducibility
  if (!missing(seed)) {
    set.seed(seed)
  }
  
  ## Create stream distance matrices
  SSN::createDistMat(ssn, o.write=TRUE)
  
  ## Extract the data.frames for the observed and prediction location data
  rawDFobs <- SSN::getSSNdata.frame(ssn, Name = "Obs")
  
  ## Extract the geographic coordinates from the SpatialStreamNetwork
  ## object and add to data.frames
  obs_data_coord <- data.frame(ssn@obspoints@SSNPoints[[1]]@point.coords)
  obs_data_coord$pid <- as.numeric(rownames(obs_data_coord))
  rawDFobs <- rawDFobs %>% left_join(obs_data_coord, by = c("pid"),
                                     keep = FALSE)
  rawDFobs$point <- "Obs" ## Create label for observed points
  
  ## Generate continous covariate at observed and prediction locations
  if (!missing(seed)) {
    set.seed(seed)
  }
  rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
  #rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
  #rawDFobs[,"X3"] <- rnorm(length(rawDFobs[,1]))
  
  ## Put the new covariates back in the SpatialStreamNetwork object
  ssn <- SSN::putSSNdata.frame(rawDFobs, ssn, Name = 'Obs')
  
  ## Simulate the response variable at observed and prediction locations
  if (!missing(seed)) {
    set.seed(seed)
  }
  sim.out <- SSN::SimulateOnSSN(ssn.object = ssn,
                                ObsSimDF = rawDFobs, ## observed data.frame
                                #PredSimDF = rawDFpred, ## prediction data.frame
                                #PredID = "preds", ## name of prediction dataset
                                formula = ~ X1,
                                coefficients = c(10, 2), ## regression coefficients
                                CorModels = c("Exponential.taildown"), ## covariance model
                                use.nugget = TRUE, ## include nugget effect
                                CorParms = c(3, 10, .1)) ## covariance parameters
  
  
  ## Extract the SpatialStreamNetwork object from the list returned by
  ## SimulateOnSSN and extract the observed and prediction site
  ## data.frames. Notice the new column Sim_Values in the data.frames
  sim.ssn <- sim.out$ssn.object
  df_sim <- SSN::getSSNdata.frame(sim.ssn,"Obs") %>%
    setDT %>%
    replicate(timesteps, ., simplify = FALSE) %>%
    rbindlist
  
  ## Expand data.frames to include t=timesteps days per location
  df_sim[, timestep := rep(1:timesteps,
                           each = (.N/timesteps))] # Set timestep variable
  
  #Create autoregressive model
  if (!missing(seed)) {
    set.seed(seed)
  }
  phi <- 0.8 ## lag 1 autocorrelation value
  ar1_sim <- nlme::corAR1(form = ~ timestep, value = phi) # can also use corExp function
  AR1 <- nlme::Initialize(ar1_sim, 
                          data = data.frame(timestep = unique(df_sim$timestep)))
  
  ## Create a vector of AR1 errors for each date and expand to all all locations
  #NB AR1 error
  epsilon <- t(chol(corMatrix(AR1))) %*% rnorm(length(unique(df_sim$timestep)), 
                                               0, 2) 
  epsilon <- rep(epsilon, each = length(unique(df_sim$locID))) +
    rnorm(length(epsilon)*length(unique(df_sim$locID)), 0, 0.25) # for all the locations
  
  epsilon_df <- data.table(timestep = rep(unique(df_sim$timestep), 
                                          each = length(unique(df_sim$locID))),
                           locID = rep(unique(df_sim$locID), 
                                       times = length(unique(df_sim$timestep))),
                           epsilon = epsilon)
  
  df_sim <- df_sim %>% left_join(epsilon_df, 
                                 by = c('timestep' = 'timestep', 
                                        'locID' = 'locID'))
  
  ## Create a new simulated response variable, y, with errors added
  df_sim$env1 <- df_sim$Sim_Values + df_sim$epsilon 
  
  #Rename columns to match the rest of the workflow
  setnames(df_sim, 
           c("locID", "coords.x1", "coords.x2", "timestep"),
           c("patch", "x", "y", "time")
  )
  
  #Standardize all values to fall between 0 and 1
  df_sim[, env1 := vegan::decostand(env1, method = "range")]
  
  if (plot) {
    ## Create line plots of the response over time for training and test datasets
    ts_plot <- ggplot(df_sim) +
      geom_line(aes(x = time, y = env1, group = patch), alpha = 0.4) +
      ylab("Simulated environment")+
      #geom_line(aes(x = timestep, y=Sim_Values, group = locID), alpha=0.4, color='red') +
      theme_bw()
    
    print(ts_plot)
  }
  
  return(df_sim[, c("patch", "time", "env1"), with=F])
}



### ~~~~~~~~~~~~~~~~~~~~~~~~~~ env_generate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
# ###############################################################################Separate checking of env_df
# } else {
#   if(all.equal(names(env_df), c("env1", "patch", "time")) != TRUE) 
#     stop("env_df must be a dataframe with columns: env1, patch, time")
# }
# ###############################################################################Separate checking of env_df


#' Generate Environment
#'
#' Generates density independent environmental conditions for metacommunity simulations
#'
#' @param landscape landscape generated by landscape_generate()
#' @param env_df optional dataframe with environmental conditions with columns: env1, patch, time
#' @param env1Scale scale of temporal environmental autocorrelation between -2 (anticorrelated) and 2 (correlated), default is 2
#' @param timesteps number of timesteps to simulate
#' @param plot option to show plot of environmental variation
#'
#' @return dataframe with x and y coordinates, time, and environmental conditions
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' env_generate(landscape_generate())
#'
#' @importFrom synchrony phase.partnered
#' @import ggplot2
#'
#' @export
#'
env_generate <- function(landscape, env1Scale = 500,
                         timesteps = 1000, spatial_autocor = TRUE, 
                         torus=TRUE, plot = TRUE) {
  
  if (inherits(landscape, 'igraph')) {
    patches <- gorder(landscape)
    landscape <- as.data.table(get.vertex.attribute(landscape))
  } else if (inherits(landscape, 'data.frame')) {
    patches <- nrow(landscape)
  }
  
  repeat {
    if (spatial_autocor) {
      env_df <- simulate_envSSN(ssn = landscape,
                                timesteps = timesteps,
                                plot = FALSE)
      
      #(to better fill 0-1 space?)
      ecum <- ecdf(env_df$env1)
      env_cum <- ecum(env_df$env1)
      env_df$env1 <- env_cum
    } else { #i.e. if spatial_autocor == FALSE
      print("This version differs from Thompson et al. 2020 in that it does not produce spatially autocorrelated environmental variables.")
      
      env_df <- lapply(1:patches, function(i) {
        #Create two time series w the phase partnered algorithm (Vasseur (2007):
        # specific autocorrelation gamma,
        # cross-correlation rho #default = 1,
        # mean mu
        # standard deviation sigma 
        env1 = phase.partnered(n = timesteps, 
                               gamma = env1Scale, 
                               mu = 0.5, 
                               sigma = 0.25)$timeseries[,1]
        
        #Standardize all values to fall between 0 and 1
        env_df <- data.frame(env1 = vegan::decostand(env1,method = "range"),
                             patch = i, 
                             time = 1:timesteps)
        return(env_df)
      }) %>% rbindlist
    }
    
    #Get values at first step for burn-in
    env.initial <- env_df[env_df$time == 1,]
    #Make sure there's a sufficient range of env values for burn-in
    range_ini <- (max(env.initial$env1)-min(env.initial$env1))
    print(range_ini)
    if(range_ini > 0.5) {break}
  }
  
  if(plot == TRUE){
    g<-ggplot2::ggplot(env_df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
      ggplot2::geom_line()+
      scale_color_viridis_d(guide = "none")
    print(g)
  }
  return(env_df)
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ env_traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#' Generate Species Env. Traits
#'
#' Generates species specific traits for density independent environmental responses
#'
#' @param species number of species to simulate
#' @param max_r intrinsic growth rate in optimal environment, can be single value or vector of length species
#' @param min_env minium environmental optima
#' @param max_env minium environmental optima
#' @param env_niche_breadth standard deviation of environmental niche breadth, can be single value or vector of length species
#' @param optima optional values of environmental optima, should be a vector of length species
#' @param optima_spacing "even" or "random" to specify how optima should be distributed
#' @param plot option to show plot of environmental variation
#'
#' @return dataframe an optima and niche breadth for each species
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' env_traits(species = 10)
#'
#' @export
#'
env_traits <- function(species, max_r = 5, min_env = 0, max_env = 1, 
                       env_niche_breadth = 0.5, optima, plot = TRUE, 
                       optima_spacing = "random"){
  
  #Determine environmental optima for species
  if (missing(optima)){
    if(optima_spacing == "even"){ #Linear increase in environmental optima
      optima <- seq(from = 0,to = 1,length = species)
    }
    if(optima_spacing == "random"){ #Random assignment of environmental optima
      optima <- runif(n = species, min = min_env, max = max_env)
    }
  } else {
    if(length(optima)!=species) stop("optima is not a vector of length species")
    if(class(optima)!="numeric") stop("optima is not a numeric vector")
  }
  
  env_traits_df <- data.frame(species = 1:species, 
                              optima = optima, 
                              env_niche_breadth = env_niche_breadth,
                              max_r = max_r)
  
  if(plot == TRUE){
    matplot(sapply(X = 1:species, FUN = function(x) {
      exp(-((env_traits_df$optima[x]-seq(min_env, max_env, length = 30))/(2*env_traits_df$env_niche_breadth[x]))^2)
    })*rep(max_r,each = 30), type = "l", lty = 1, ylab = "r", xlab = "environment", ylim = c(0,max(max_r)))
    
  }
  return(env_traits_df)
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ species_int_mat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#' Generate Species Interaction Matrix
#'
#' Generates density dependent matrix of per capita competition
#'
#' @param species number of species to simulate
#' @param intra intraspecific competition coefficient, single value or vector of length species
#' @param min_inter min interspecific comp. coefficient
#' @param max_inter max interspecific comp. coefficient
#' @param int_mat option to supply externally generated competition matrix
#' @param comp_scaler value to multiply all competition coefficients by
#' @param plot option to show plot of competition coefficients
#'
#' @return species interaction matrix
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' env_traits(species = 10)
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#'

#Comment Mathis: the default parameters generate a mixed competitive structure,
#whereby alpha_ijvalues of αij are drawn from a uniform 
#distribution in the range [0, 1.5], resulting in a combination of 
#species pairs for which competition is stabilising, destabilising, 
#or where one of the two is competitively dominant
species_int_mat <- function(species, intra = 1, min_inter = 0, max_inter = 1.5, 
                            int_mat, comp_scaler = 0.05, plot = TRUE){
  if (missing(int_mat)){
    int_mat <- matrix(runif(n = species*species, 
                            min = min_inter,
                            max = max_inter), 
                      nrow = species, ncol = species)
    diag(int_mat) <- intra
    
    #scale all competition coefficients by multiplying them by 0.05, 
    #to allow for higher equilibrium abundances.
    int_mat <- int_mat * comp_scaler 
  } else {
    if (is.matrix(int_mat) == FALSE)
      stop("int_mat must be a matrix")
    if (sum(dim(int_mat) != c(species,species))>0)
      stop("int_mat must be a matrix with a row and column for each species")
    if (is.numeric(int_mat) == FALSE)
      stop("int_mat must be numeric")
  }
  
  if (plot == TRUE){
    colnames(int_mat)<- 1:species
    g <- as.data.frame(int_mat) %>%
      dplyr::mutate(i = 1:species) %>%
      tidyr::gather(key = j, value = competition, -i) %>%
      dplyr::mutate(i = as.numeric(as.character(i)),
                    j = as.numeric(as.character(j))) %>%
      ggplot2::ggplot(ggplot2::aes(x = i, y = j, fill = competition))+
      ggplot2::geom_tile()+
      scale_fill_viridis_c(option = "E")
    
    print(g)
  }
  return(int_mat)
}

### 

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ compare_predobs_env_traits ~~~~~~~~~~~~~~~~~~####

compare_predobs_env_traits <- function(MCsim, subn) {
  
  dynamics_df <- MCsim$dynamics_df
  env_traits_df <- MCsim$env_traits_df
  
  dynamics_df_sim <- setDT(dynamics_df)[time>=0,]
  
  dynamics_df_sim <- dynamics_df_sim[order(time),
                                     Nlag1 := shift(N, n=1, type='lag'),
                                     by=.(species, patch)] %>%
    .[Nlag1 > 0, r_obs := N/Nlag1] %>%
    .[, N_stand := (N-min(N))/(max(N)-min(N))]
  
  min_env <- min(dynamics_df_sim$env)
  max_env <- max(dynamics_df_sim$env)
  
  env_range <- data.table(env = seq(min_env, max_env, length = 30)) 
  species <- length(unique(dynamics_df_sim$species))
  
  env_traits_curves <- 
    cbind(
      env_range,
      lapply(X = 1:species, FUN = function(x) {
        r <- env_traits_df$max_r*exp(-((env_traits_df$optima[x]-env_range$env)/
                                         (2*env_traits_df$env_niche_breadth[x]))^2)
        data.table(species=x,
                   r)
      }) %>%
        rbindlist
    ) 
  env_traits_curves[, r_stand := (r-min(r))/(max(r)-min(r))]
  
  if (!missing(subn)) {
    dynamics_df_sim <- dynamics_df_sim[sample(dynamics_df_sim[, .N],
                                              subn,
                                              replace=F),]
  }
  
  
  # rcompare_plot<- ggplot(data=dynamics_df_sim,
  #                        aes(x=100*env, y=r_obs, color=factor(species))) + 
  #   geom_point(alpha=1/10) +
  #   geom_line(data=env_traits_curves, aes(y=r), size=1.5, linetype='dashed') + 
  #   facet_wrap(~species) +
  #   geom_smooth(size=1.5) + 
  #   scale_y_sqrt(breaks=c(0, 1, 2, 3, 4, 5, 10, 20), limits=c(0,20)) +
  #   theme_bw()
  # 
  
  rcompare_plot<- ggplot(data=dynamics_df_sim,
                         aes(x=100*env, y=N_stand, color=factor(species))) + 
    geom_point(alpha=1/10) +
    geom_line(data=env_traits_curves, aes(y=r_stand),
              color='black', size=1.5, linetype='dashed') + 
    facet_wrap(~species) +
    geom_smooth(size=1.5, color='black') + 
    scale_y_sqrt(breaks=seq(0,1, 0.1)) +
    theme_bw()
  
  return(rcompare_plot)
}
