subpatch <- 20

#--------------- Get scenarios output data ------------------------------------
scenario_list <- c(
  # 'sim_sp1_kernel_0_dispersal_0',
  # 'sim_sp1_kernel_0.1_dispersal_0.01',
  # 'sim_sp1_kernel_0.1_dispersal_0.1',
  # 'sim_sp1_kernel_0.1_dispersal_0.5',
  #'sim_sp2_kernel_0_dispersal_0',
  'sim_sp2_kernel_0.1_dispersal_0.01_nbarriers_0',
  'sim_sp2_kernel_0.1_dispersal_0.1_nbarriers_0',
  'sim_sp2_kernel_0.1_dispersal_0.5_nbarriers_0',
  'sim_sp2_kernel_0.1_dispersal_0.01_nbarriers_10',
  'sim_sp2_kernel_0.1_dispersal_0.1_nbarriers_10',
  'sim_sp2_kernel_0.1_dispersal_0.5_nbarriers_10'
)


simdt <- lapply(scenario_list, function(scenario) {
  print(scenario)
  MCsim <- tar_read_raw(scenario)
  scenario_split <- strsplit(scenario, '_')[[1]]
  
  return(
    data.table(
      cbind(MCsim$dynamics_df),
      scenario = scenario,
      sp_num = as.numeric(substr(scenario_split[2], 3,4)),
      kernel = as.numeric(scenario_split[4]),
      dispersal_p = as.numeric(scenario_split[6]),
      nbarriers = as.numeric(scenario_split[8])
    )
  )
}) %>% 
  rbindlist
dynamics_df_sim <- setDT(simdt)[time>=0,]

#----------------------- Sample data -------------------------------------------
compute_gampreds <- function(in_dynamics_df_sim, k) {
  env_preds <- data.frame(env=seq(0, 1, 0.01))
  scenario_species_combinations <- unique(in_dynamics_df_sim[
    , c('scenario', 'species', 'dispersal_p', 'nbarriers'), with=F])
  
  gam_preds <- mapply(function(in_scenario, in_species) {
    in_dt <- in_dynamics_df_sim[scenario==in_scenario & 
                                  species==in_species,] 
    
    gam_form <- (N ~ s(env, bs = "cs", k=k))
    gam_trained_50 <- qgam(gam_form, #cubic splines are good when lots of data points
                           data=in_dt, 
                           qu =0.5,
                           multicore=T)
    
    gam_trained_10 <- qgam(gam_form, #cubic splines are good when lots of data points
                           data=in_dt,
                           qu =0.1,
                           multicore=T)

    gam_trained_90 <- qgam(gam_form, #cubic splines are good when lots of data points
                           data=in_dt,
                           qu =0.9,
                           multicore=T)
    
    preds <- data.table(
      fit_50 = predict(gam_trained_50, env_preds, se.fit=FALSE),
      fit_10 = predict(gam_trained_10, env_preds, se.fit=FALSE),
      fit_90 = predict(gam_trained_90, env_preds, se.fit=FALSE)
    )
    
    # setnames(preds, c('predN', 'predN_se'))
    # preds$predN_95CIlow <- (preds$predN - 1.96 * preds$predN_se)
    # preds$predN_95CIhigh <- (preds$predN + 1.96 * preds$predN_se)
    
    return(
      cbind(
        env=env_preds,
        preds,
        data.table(
          scenario = in_scenario,
          species = in_species)
      )
    )
  },
  in_scenario = scenario_species_combinations$scenario, 
  in_species = scenario_species_combinations$species,
  SIMPLIFY=FALSE
  ) %>%
    rbindlist
  
  gam_preds[, `:=`(fit_50_stand = 100*fit_50/max(fit_50),
                   fit_10_stand = 100*fit_10/max(fit_50),
                   fit_90_stand = 100*fit_90/max(fit_50)
  )]
  predcols = c('fit_50_stand', 'fit_10_stand', 'fit_90_stand')
  
  gam_preds <- gam_preds[, (predcols) := lapply(.SD, 
                                                function(x) fifelse(x<0, 0, x)), 
                         .SDcols = predcols] %>%
    merge(scenario_species_combinations, by=c('scenario', 'species'))
  
  return(gam_preds)
}


plot_scenarios <- function(in_dynamics_df_sim, in_gam_preds,
                           display_points=FALSE, gam=TRUE) {
  in_dynamics_df_sim <- in_dynamics_df_sim[order(time),
                                           `:=`(
                                             Nlag1 = shift(N, n=1, type='lag'),
                                             envlag1 = shift(env, n=1, type='lag')
                                           ),
                                           by=.(species, patch, scenario,
                                                dispersal_p, nbarriers)] %>%
    .[Nlag1 > 0, r_obs := N/Nlag1] %>%
    .[, N_stand := N/max(in_gam_preds$fit_50)]
  #.[, N_stand := (N-min(N))/(max(N)-min(N))]
  
  min_env <- min(in_dynamics_df_sim$env)
  max_env <- max(in_dynamics_df_sim$env)
  
  env_range <- data.table(env = seq(min_env, max_env, length = 30)) 
  species <- length(unique(in_dynamics_df_sim$species))
  
  
  env_traits_df <- lapply(scenario_list, function(scenario) {
    MCsim <- tar_read_raw(scenario)
    return(MCsim$env_traits_df)
  }) %>% 
    rbindlist %>%
    unique
  
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
  
  
  rcompare_plot<- ggplot(data=in_dynamics_df_sim,
                         aes(x=100*env, group=factor(species), 
                             fill=factor(species),
                             color=factor(species))) + 
    geom_line(data=env_traits_curves, 
              aes(x=100*env, y=100*r_stand),
              #color='black', 
              size=1, linetype='dashed', alpha=0.75) +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,120),
                       expand=c(0,0), 
                       breaks=seq(0,100, 25)) +
    scale_color_manual(values=c('#fbbd08ff', '#094561bf')) +
    scale_fill_manual(values=c('#fbbd08ff', '#094561bf')) +
    labs(x='Environment', 
         y=stringr::str_wrap('Relative population size (solid) and habitat suitability (dashed)',
                   width=40)
         )+
    facet_grid(dispersal_p~nbarriers, labeller=label_both) +
    theme_bw() +
    theme(legend.position = 'none')
  
  if (display_points) {
    rcompare_plot <- rcompare_plot + geom_point(alpha=0.1, aes(y=100*N_stand))
  }
  
  if (gam) {
    rcompare_plot <- rcompare_plot + geom_ribbon(data=in_gam_preds, 
                                                 aes(ymin=fit_10_stand, ymax=fit_90_stand),
                                                 alpha=0.3, color=NA, size=1) +
      geom_line(data=in_gam_preds, 
                aes(y=fit_50_stand),
                size=1, alpha=0.75) 
  }
  
  return(rcompare_plot)
}

set.seed(0)
patch_sample <- sample(unique(dynamics_df_sim$patch), subpatch)
dynamics_df_sim_10patches <- dynamics_df_sim[patch %in% patch_sample,]
# if (!missing(subpatch)) {
#   patch_sample <- sample(unique(dynamics_df_sim$patch), subpatch)
#   dynamics_df_sim <- dynamics_df_sim[patch %in% patch_sample,]
# }
# if (!missing(subn)) {
#   dynamics_df_sim <- dynamics_df_sim[sample(dynamics_df_sim[, .N],
#                                             subn,
#                                             replace=F),]
# }

gampreds_10patches <- compute_gampreds(
  in_dynamics_df_sim = dynamics_df_sim_10patches,
  k=3)

plot_10patches <- plot_scenarios(
  in_gam_preds = gampreds_10patches[dispersal_p %in% c(0.1, 0.5),],
  in_dynamics_df_sim = dynamics_df_sim_10patches[dispersal_p %in% c(0.1, 0.5),],
  display_points = FALSE,
  gam=TRUE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust=-0.5)
        )

dynamics_df_sim_1patch <- dynamics_df_sim[patch == 5,]
# gampreds_1patch <- compute_gampreds(
#   in_dynamics_df_sim = dynamics_df_sim_1patch,
#   k=1)
plot_1patch <- plot_scenarios(
  in_gam_preds = gampreds_1patch[dispersal_p == 0.1,],
  in_dynamics_df_sim = dynamics_df_sim_1patch[dispersal_p == 0.1,],
  display_points = TRUE,
  gam=FALSE) +
  theme(axis.title.y = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.background = element_blank()
  ) + 
  geom_smooth(aes(y=100*N_stand), se=FALSE, method='loess', 
              span=1, fullrange = T)

plot_10patches / plot_1patch + plot_layout(height=c(3,1))
