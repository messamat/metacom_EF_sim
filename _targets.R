source("R/packages.R")
source("R/sim_setup.R")
source("R/MC_simulate.R")

rootdir = rprojroot::find_root(rprojroot::has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')


#--- Define targets plan ------------------------------------------------------
list(
  tar_target(
    OCNigraph,
    generate_OCNigraph(patches=100, 
                       cellsize = 0.5, 
                       dimX = 25, 
                       dimY = 25, 
                       plot = T)
  ),
  
  tar_target(
    env_traits_df,
    env_traits(
      species = 1,
      max_r = 5,
      min_env = 0,
      max_env = 1,
      env_niche_breadth = 0.5,
      plot = F,
      optima_spacing = 'random'
    )
  ),
  
  tarchetypes::tar_map(
    values = tibble(
      in_dispersal = c(0, 0.01, 0.05)
    ),
    
    tar_target(
      sim_sp1,
      simulate_MC(timesteps = 1200, 
                  burn_in = 800,
                  initialization = 200,
                  intra = 1,
                  min_inter = 0, 
                  max_inter = 0.5,
                  dispersal = in_dispersal,
                  landscape = OCNigraph,
                  env_traits_df = env_traits_df,
                  )
    ),
    
    tar_target(
      FER_sp1,
      compare_predobs_env_traits(sim_sp1, subn=10000),
    )
  )
)
