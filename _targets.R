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
    sim_sp1_0dis_randenv,
    simulate_MC(species = 1,
                timesteps = 1200, 
                burn_in = 800,
                initialization = 200,
                intra = 1,
                max_r = 5,
                env_niche_breadth = 0.5,
                min_inter = 0, 
                max_inter = 0.5,
                dispersal = 0,
                landscape = OCNigraph)
  ),
  
  tar_target(
    FER_sp1_0dis_randenv,
    compare_predobs_env_traits(sim_sp1_0dis_randenv),
  )
)
