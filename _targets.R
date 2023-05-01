source("R/packages.R")
source("R/sim_setup.R")
source("R/MC_simulate.R")

#--- Define targets plan ------------------------------------------------------
list(
  tar_target(
    OCN_formatted_list,
    generate_OCN_formatted(
      patches = 100, 
      cellsize = 0.5, 
      dimX = 25, 
      dimY = 25, 
      plot = T,
      out_format = list('SSN', 'igraph'),
      out_SSNdir = file.path(paste0('simSSN_', 
                                    format(Sys.Date(), "%Y%m%d%H%M")))
    )
  ),
  
  tar_target(
    env_traits_df,
    env_traits(
      species = 2,
      max_r = 5,
      min_env = 0,
      max_env = 1,
      env_niche_breadth = 0.25,
      plot = T,
      optima_spacing = c(0.40, 0.60)
    )
  ),
  
  tar_target(
    env_df,
    env_generate(
      landscape = OCN_formatted_list$SSN,
      env1Scale = 500, 
      timesteps = 2000, #Includes burn-in
      spatial_autocor = TRUE, 
      plot = TRUE
    )
  ),
  
  tar_target(
    env_map,
    plot_SSNenv(
      ssn = OCN_formatted_list$SSN,
      df_sim = env_df,
      timesteps_to_plot = 1:10,
      tscol = 'time'
    )
  ),
  
  tar_target(
    env_plot,
    ggplot(env_df[time %in% 1:100 & patch %in% 1:20,],
           aes(x=time, y=env1, group=patch, color=factor(patch))) +
      geom_line(alpha=0.7) +
      theme_bw() +
      theme(legend.position = 'none')
  ),
  
  tar_target(
    dist_mat,
    compute_distmat(landscape = OCN_formatted_list$igraph)
  ),
  
  tar_target(
    int_mat,
    species_int_mat(species = nrow(env_traits_df),
                    intra = 1,
                    min_inter = 0,
                    max_inter = 0.5, 
                    comp_scaler = 0.05,
                    plot = TRUE)
  ),
  
  tarchetypes::tar_map(
    values = expand.grid(
      data.table(
        in_kernel = 0.1,
        in_dispersal = c(0.01, 0.1, 0.5)
      )
    ) %>%
      setDT %>%
      rbind(data.table(in_kernel = 0, 
                       in_dispersal = 0)
            ) %>%
      .[, scenario_name := paste0("kernel_", in_kernel, 
                                  "_dispersal_", in_dispersal)] %>%
      unique
    ,
    
    names = "scenario_name",
      
    tar_target(
      disp_mat,
      dispersal_matrix(
        landscape = OCN_formatted_list$igraph,
        kernel_exp = in_kernel, 
        plot = TRUE
      )
    ),
    
    tar_target(
      sim_sp2,
      simulate_MC(#timesteps = 1200, 
        burn_in = 800,
        initialization = 200,
        intra = 1,
        min_inter = 0, 
        max_inter = 1,
        dispersal = in_dispersal,
        landscape = OCN_formatted_list$igraph,
        disp_mat = disp_mat,
        env_df = env_df,
        env_traits_df = env_traits_df,
        int_mat = int_mat,
      )
    ),
    
    tar_target(
      sim_sp1,
      simulate_MC(#timesteps = 1200, 
        burn_in = 800,
        initialization = 200,
        intra = 1,
        min_inter = 0, 
        max_inter = 1,
        dispersal = in_dispersal,
        landscape = OCN_formatted_list$igraph,
        disp_mat = disp_mat,
        env_df = env_df,
        env_traits_df = env_traits_df,
        species = 1
      )
    ),
    
    #Compare with one vs. two species and with three different levels of dispersal
    tar_target(
      FER_sp2,
      compare_predobs_env_traits(sim_sp2, subn=10000),
    )
    
    #Add barriers
    
    
    #Check at the scale of a site over time vs. a snapshot of all sites
  
    
  )
)
