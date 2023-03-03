### ~~~~~~~~~~~~~~~~~~~~~~~~~~ convert_OCN_to_igraph ~~~~~~~~~~~~~~~~~~~~~~####
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
convert_OCN_to_igraph <- function(OCN) {
  # Get geographic position
  nodes_coords <- expand.grid(c(1:OCN$dimX), 
                              c(1:OCN$dimY))[which(OCN$FD$toRN!=0),] %>%
    data.table %>%
    .[, lapply(.SD, function(x) x*OCN$cellsize)] %>%
    setnames(c("X", "Y"))

  #Create graph
  #(besides outlets which do not have downstream edges)
  graph <- data.table(from = 1:OCN$RN$nNodes,
                      to = OCN$RN$downNode,
                      weight = OCN$RN$leng,
                      DA = OCN$RN$A) %>%
    cbind(nodes_coords) %>%
    .[to != 0,] %>% #remove outlet
    graph_from_data_frame
  
  return(graph)
}

### ~~~~~~~~~~~~~~~~~~~~~~~~~~ generate_OCNigraph ~~~~~~~~~~~~~~~~~~~~~~~~####

#' Generate OCN
#' @author adapted from Claire Jacquet (OID: 10.1111/OIK.09372)
#' 
#' 
generate_OCNigraph <- function(patches, cellsize = 0.5, dimX = 25, dimY = 25,
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
  
  graph <- convert_OCN_to_igraph(OCN)
  
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
dispersal_matrix <- function(landscape, torus = TRUE, disp_mat, 
                             kernel_exp = 0.1, plot = TRUE){
  
  if (inherits(landscape, 'igraph')) {
    dist_mat <- igraph::distances(landscape)
  } else {
    #Compute distance among sites
    if(torus == TRUE){
      dist_mat <- as.matrix(som.nn::dist.torus(coors = landscape))
    } else {
      dist_mat <- as.matrix(dist(landscape))
    }
  }
  
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


### ~~~~~~~~~~~~~~~~~~~~~~~~~~ env_generate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
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
env_generate <- function(landscape, env_df, env1Scale = 2, 
                         timesteps = 1000, plot = TRUE){
  if (missing(env_df)){
    repeat {
      env_df <- data.frame()
      
      if (inherits(landscape, 'igraph')) {
        patches <- gorder(landscape)
      } else if (inherits(landscape, 'data.frame')) {
        patches <- nrow(landscape)
      }
      
      for(i in 1:patches){
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
        env_df <- rbind(env_df, 
                        data.frame(env1 = vegan::decostand(env1,method = "range"),
                                   patch = i, 
                                   time = 1:timesteps))
      }
      #Get values at first step for burn-in
      env.initial <- env_df[env_df$time == 1,]
      
      #Make sure there's a sufficient range of env values for burn-in
      if((max(env.initial$env1)-min(env.initial$env1)) > 0.6) {break}
    }
  } else {
    if(all.equal(names(env_df), c("env1", "patch", "time")) != TRUE) 
      stop("env_df must be a dataframe with columns: env1, patch, time")
  }
  
  if(plot == TRUE){
    g<-ggplot2::ggplot(env_df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
      ggplot2::geom_line()+
      scale_color_viridis_d(guide = "none")
    print(g)
  }
  
  return(env_df)
  print("This version differs from Thompson et al. 2020 in that it does not produce spatially autocorrelated environmental variables.")
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
    .[!is.na(r_obs)] %>%
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
    dynamics_df_sim <- dynamics_df_sim[sample(dynamics_df_sim[!is.na(r_obs), .N],
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
