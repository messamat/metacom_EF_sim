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
    basic_sim,
    simulate_MC(species = 2,
                landscape = OCNigraph)
  )
)