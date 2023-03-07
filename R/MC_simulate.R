#' Simulate Metacommunity Dynamics
#'
#' @param patches number of patches to simulate
#' @param species number of species to simulate
#' @param dispersal dispersal probability between 0 and 1
#' @param plot option to show plot of landscape
#' @param torus whether to model the landscape as a torus
#' @param kernel_exp the exponential rate at which dispersal decreases as a function of the distance between patches
#' @param env1Scale scale of temporal environmental autocorrelation between -2 (anticorrelated) and 2 (correlated), default is 2
#' @param timesteps number of timesteps to simulate
#' @param burn_in length of burn in period
#' @param initialization length of initial period before environmental change begins
#' @param max_r intrinsic growth rate in optimal environment, can be single value or vector of length species
#' @param min_env minium environmental optima
#' @param max_env minium environmental optima
#' @param env_niche_breadth standard deviation of environmental niche breadth, can be single value or vector of length species
#' @param optima_spacing "even" or "random" to specify how optima should be distributed
#' @param intra intraspecific competition coefficient, single value or vector of length species
#' @param min_inter min interspecific comp. coefficient
#' @param max_inter max interspecific comp. coefficient
#' @param comp_scaler value to multiply all competition coefficients by
#' @param extirp_prob probability of local extirpation for each population in each time step (should be a very small value, e.g. 0 or 0.002)
#'
#' @param landscape optional dataframe with x and y columns for patch coordinates
#' @param disp_mat optional matrix with each column specifying the probability that an individual disperses to each other patch (row)
#' @param env_df optional dataframe with environmental conditions with columns: env1, patch, time
#' @param env_optima optional values of environmental optima, should be a vector of length species
#' @param int_mat optional externally generated competition matrix

#' @return list that includes metacommunity dynamics, landscape coordinates, environmental conditions, species environmental traits, dispersal matrix, and the competition matrix
#'
#' @author Patrick L. Thompson, \email{patrick.thompson@@zoology.ubc.ca}
#'
#' @examples
#' simulate_MC(patches = 6, species = 10, dispersal = 0.001, min_inter = 1, max_inter = 1, env_niche_breadth = 10)

#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
simulate_MC <- function(
  species, patches, dispersal = 0.01,
  plot = TRUE,
  torus = FALSE, kernel_exp = 0.1,
  env1Scale = 500, timesteps = 1200, burn_in = 800, 
  initialization = 200, max_r = 5,
  min_env = 0, max_env = 1, env_niche_breadth = 0.5, optima_spacing = "random",
  intra = 1, min_inter = 0, max_inter = 1, comp_scaler = 0.05,
  extirp_prob = 0,
  landscape, disp_mat, env_df, env_traits_df, int_mat){
  

  #Get landscape structure (coordinate of patches)
  if (missing(landscape)){
    landscape <- landscape_generate(patches = patches, 
                                    plot = plot)
  } else {
    if (inherits(landscape, 'igraph')) {
      patches <- gorder(landscape)
      
    } else if (inherits(landscape, c('data.frame',
                                     "SpatialStreamNetwork"))) {
      patches <- nrow(landscape)
    }
  }
  
  #Get dispersal matrix (defining the mean relative proportion of individuals 
  #that disperse from each site to every other site)
  if (missing(disp_mat)){
    disp_mat <- dispersal_matrix(landscape = landscape, 
                                 torus = torus, 
                                 kernel_exp = kernel_exp, 
                                 plot = plot)
  } 
  
  #Get environmental conditions
  if (missing(env_df)){
    env_df <- env_generate(landscape = landscape, 
                           env1Scale = env1Scale, 
                           timesteps = timesteps+burn_in, 
                           plot = plot)
  } else {
    timesteps <- length(unique(env_df$time))- burn_in
  }
  
  #Get species traits matrix 
  #(max growth rate + abiotic environment optima and niche breadth
  if (missing(env_traits_df)){
    env_traits_df <- env_traits(species = species, 
                                max_r = max_r, 
                                min_env = min_env, 
                                max_env = max_env, 
                                env_niche_breadth = env_niche_breadth, 
                                optima_spacing = optima_spacing, 
                                plot = plot)
  } else {
    species <- nrow(env_traits_df)
  }
  
  #Get species interaction matrix
  if (missing(int_mat)){
    int_mat <- species_int_mat(species = species, 
                               intra = intra, 
                               min_inter = min_inter, 
                               max_inter = max_inter, 
                               comp_scaler = comp_scaler, 
                               plot = TRUE)
  } 
  
  dynamics_df <- data.frame()
  
  #Seed all species in all patches based on Poisson distribution with lambda=0.5  
  N <- matrix(rpois(n = species*patches, 
                    lambda = 0.5), 
              nrow = patches,
              ncol = species)
  
  pb <- txtProgressBar(min = 0, 
                       max = initialization + burn_in + timesteps, 
                       style = 3)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in 1:(initialization + burn_in + timesteps)){ 

    if(i <= initialization){ #for initialization steps
      if(i %in% seq(10,100, by = 10)){ #Every 10 time steps
        N <- N + matrix(rpois(n = species*patches, lambda = 0.5), 
                        nrow = patches, ncol = species) #Recruitment event with lambda = 0.5
      }
      env <- env_df$env1[env_df$time == 1] #Keep the same environmental conditions
    } else {
      env <- env_df$env1[env_df$time == (i-initialization)] #Set environmental conditions to the right time step
    }
    
    #~~~~ Compute density-independent growth rate (eq. 3; p1319) ~~~~~~~~~~~~~~~
    #Compute for each species the difference between its env optima and the env
    #Each column is a species, each row is a patch
    z_minus_envxt <- t((env_traits_df$optima - 
                          matrix(rep(env, each = species),
                                 nrow = species, ncol = patches)))
    #Fill in eq. 3: density-independent growth for each species in all patches
    r <- max_r * exp(-(z_minus_envxt/(2*env_traits_df$env_niche_breadth))^2)
    
    #~~~~ Compute the population size before accounting for dispersal ~~~~~~~~~~
    #(eq. 2; p1319)
    N_hat <- N*r/(1+N%*%int_mat)
    N_hat[N_hat < 0] <- 0
    N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), 
                    ncol = species, nrow = patches)
    
    if (dispersal > 0) {
      #~~~~~~~~~ Compute emigration matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #(equal probability of dispersal for each species)
      E <- matrix(rbinom(n = patches * species, 
                         size = N_hat, 
                         prob = rep(dispersal, each = patches)), 
                  nrow = patches,
                  ncol = species)
      
      #~~~~~~~~~ Compute immigration matrix (with equation 4) ~~~~~~~~~~~~~~~~~~~~
      dispSP <- colSums(E) #Total number of emigrating individuals for each species
      
      I_hat_raw <- disp_mat%*%E 
      
      #Probability that a dispersing individual of each species immigrates to each site
      I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw)) 
      I_hat[is.nan(I_hat)] <- 1
      
      #For each species
      I <- sapply(1:species, function(x) {
        if(dispSP[x]>0){
          table(factor(sample(x = patches, 
                              size = dispSP[x], 
                              replace = TRUE, 
                              prob = I_hat[,x]), 
                       levels = 1:patches))
        } else {
          rep(0, patches)}
      })
      
      #Compute final number of individuals per species per patch (eq. 1)
      N <- N_hat - E + I
    } else {
      N <- N_hat
    }

    #Implement stochastic extirpation 
    #(in the case of a competition-colonisation trade-off set-up 
    #see Simulation runs (4), p1321)
    N[rbinom(n = species * patches, 
             size = 1, 
             prob = extirp_prob)>0] <- 0
    
    
    dynamics_df <- rbind(dynamics_df, 
                         data.frame(N = c(N), 
                                    patch = 1:patches, 
                                    species = rep(1:species, each = patches),
                                    env = env, 
                                    time = i-initialization-burn_in))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  dynamics_df <- left_join(dynamics_df, env_traits_df)
  env_df$time_run <- env_df$time - burn_in
  
  env_df_init <- data.frame(env1 = env_df$env1[env_df$time == 1], 
                            patch = 1:patches, 
                            time = NA, 
                            time_run = rep(seq(-(burn_in + initialization), 
                                               -burn_in), 
                                           each = patches))
  env_df <- rbind(env_df_init, env_df)
  
  if(plot == TRUE){
    sample_patches <- sample(1:patches, 
                             size = min(c(patches,6)), 
                             replace = FALSE)
    g <- dynamics_df %>%
      filter(time %in% seq(min(dynamics_df$time),
                           max(dynamics_df$time), 
                           by =10)) %>%
      filter(patch %in% sample_patches) %>%
      ggplot(aes(x = time, y = N, group = species, color = optima))+
      geom_line()+
      facet_wrap(~patch)+
      scale_color_viridis_c()+
      geom_path(data = filter(env_df, patch %in% sample_patches), 
                aes(y = -5, x = time_run, color = env1, group = NULL), 
                size = 3)
    print(g)
  }
  
  return(list(dynamics_df = dynamics_df, landscape = landscape, 
              env_df = env_df,env_traits_df = env_traits_df, 
              disp_mat = disp_mat, int_mat = int_mat))
}
