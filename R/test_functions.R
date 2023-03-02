tar_load(basic_sim)
dynamics.df <- basic_sim$dynamics.df
env_traits.df <- basic_sim$env_traits.df

compare_predobs_env_traits <- function(dynamics.df, env_traits.df) {
  setDT(dynamics.df)[, Nlag1 := shift(N, n=1, type='lag')] %>%
    .[Nlag1 > 0, r_obs := N/Nlag1]
  
  min_env <- min(dynamics.df$env)
  max_env <- max(dynamics.df$env)
  
  
  env_range <- data.table(env = seq(min_env, max_env, length = 30)) 
  env_traits_curves <- 
    cbind(
      env_range,
      lapply(X = 1:species, FUN = function(x) {
        r <- env_traits.df$max_r*exp(-((env_traits.df$optima[x]-env_range$env)/
                                         (2*env_traits.df$env_niche_breadth[x]))^2)
        data.table(species=x,
                   r)
      }) %>%
        rbindlist
    )

  dynamics.df_sub <- dynamics.df[sample(dynamics.df[!is.na(r_obs), .N], 
                                        10000,
                                        replace=F),]
  
  ggplot(data=env_traits_curves,
         aes(x=30*env, y=r, color=factor(species))) + 
    geom_line()
  
  ggplot(data=dynamics.df,
         aes(x=30*env, y=r_obs, color=factor(species))) + 
    geom_point(alpha=1/5) +
    geom_line(data=env_traits_curves, aes(y=r)) + 
    facet_wrap(~species) +
    geom_smooth(color='black') + 
    scale_y_sqrt(breaks=c(0, 1, 2, 3, 4, 5, 10, 20), limits=c(0,20))
}