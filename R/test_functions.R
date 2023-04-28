subpatch <- 10

scenario_list <- c(
  'sim_sp1_kernel_0_dispersal_0',
  'sim_sp1_kernel_0.1_dispersal_0.01',
  'sim_sp1_kernel_0.1_dispersal_0.1',
  'sim_sp1_kernel_0.1_dispersal_0.5',
  'sim_sp2_kernel_0_dispersal_0',
  'sim_sp2_kernel_0.1_dispersal_0.01',
  'sim_sp2_kernel_0.1_dispersal_0.1',
  'sim_sp2_kernel_0.1_dispersal_0.5'
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
      dispersal_p = as.numeric(scenario_split[6])
    )
  )
}) %>% 
  rbindlist

dynamics_df_sim <- setDT(simdt)[time>=0,]

patch_sample <- sample(unique(dynamics_df_sim$patch), subpatch)
dynamics_df_sim <- dynamics_df_sim[patch %in% patch_sample,]
# if (!missing(subpatch)) {
#   patch_sample <- sample(unique(dynamics_df_sim$patch), subpatch)
#   dynamics_df_sim <- dynamics_df_sim[patch %in% patch_sample,]
# }
# if (!missing(subn)) {
#   dynamics_df_sim <- dynamics_df_sim[sample(dynamics_df_sim[, .N],
#                                             subn,
#                                             replace=F),]
# }

#Compute GAM by scenario and species
env_preds <- data.frame(env=seq(0, 1, 0.01))
scenario_species_combinations <- unique(dynamics_df_sim[
  , c('scenario', 'species'), with=F])

gam_preds <- mapply(function(in_scenario, in_species) {
  in_dt <- dynamics_df_sim[scenario==in_scenario & 
                             species==in_species,] 
  
  gam_form <- (N ~ s(env, bs = "cs", k=4))
  gam_trained_50 <- qgam(gam_form, #cubic splines are good when lots of data points
                      data=in_dt, 
                      qu =0.5)
  
  gam_trained_10 <- qgam(gam_form, #cubic splines are good when lots of data points
                         data=in_dt, 
                         qu =0.1)
  
  gam_trained_90 <- qgam(gam_form, #cubic splines are good when lots of data points
                         data=in_dt, 
                         qu =0.9)
  
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

gam_preds[, `:=`(predN_stand = 100*predN/max(predN),
                 predN_95CIlow_stand = 100*predN_95CIlow/max(predN),
                 predN_95CIhigh_stand = 100*predN_95CIhigh/max(predN)
                 )]

dynamics_df_sim <- dynamics_df_sim[order(time),
                                   `:=`(
                                     Nlag1 = shift(N, n=1, type='lag'),
                                     envlag1 = shift(env, n=1, type='lag')
                                   ),
                                   by=.(species, patch, scenario)] %>%
  .[Nlag1 > 0, r_obs := N/Nlag1] %>%
  .[, N_stand := N/max(gam_preds$predN)]
  #.[, N_stand := (N-min(N))/(max(N)-min(N))]

min_env <- min(dynamics_df_sim$env)
max_env <- max(dynamics_df_sim$env)

env_range <- data.table(env = seq(min_env, max_env, length = 30)) 
species <- length(unique(dynamics_df_sim$species))


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


rcompare_plot<- ggplot(data=dynamics_df_sim,
                       aes(x=env, group=factor(species), color=factor(species))) + 
  #geom_point(alpha=1/10) +
  geom_ribbon(data=gam_preds, aes(ymin=predN_95CIlow_stand, 
                                  ymax=predN_95CIhigh_stand)) +
  geom_line(data=gam_preds, aes(y=predN_stand)) +
  facet_wrap(~scenario) +
  geom_line(data=env_traits_curves, aes(x=env, y=100*r_stand),
            #color='black', 
            size=1.5, linetype='dashed') +
  theme_bw()
 

  #scale_y_sqrt(breaks=seq(0,1, 0.1)) +

rcompare_plot


# rcompare_plot<- ggplot(data=dynamics_df_sim,
#                        aes(x=100*env, y=r_obs, color=factor(species))) + 
#   geom_point(alpha=1/10) +
#   geom_line(data=env_traits_curves, aes(y=r), size=1.5, linetype='dashed') + 
#   facet_wrap(~species) +
#   geom_smooth(size=1.5) + 
#   scale_y_sqrt(breaks=c(0, 1, 2, 3, 4, 5, 10, 20), limits=c(0,20)) +
#   theme_bw()
# 







