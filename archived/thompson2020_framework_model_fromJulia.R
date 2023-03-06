source("model_packages.R")

reps = 10

for (rep in 1:reps) {
  env_df <- fread(file.path(datdir_landscape, 
                               paste0("env_", rep, ".csv")))
  disp_mat <- fread(file.path(datdir_landscape, 
                                 paste0("disp_mat_", rep, ".csv")),
                    header=T)
  disp_mat <- as.array(disp_mat[, 2:ncol(disp_mat)]) ########TO figure out################## ERROR
  landscape <- fread(file.path(datdir_landscape, 
                                  paste0("landscape_", rep, ".csv")))
  time <- fread(file.path(datdir_landscape, 
                             paste0("time_", rep, ".csv")))
  
  S <- 50 #Number of species
  dominants <- trunc(S * 0.3) #Number of dominant species (see p1321 in Thompson et al. 2020)
  M <- dim(disp_mat)[1] #Number of patches (e.g., eq4)
  
  #################################### JUST CASTING ENV_DV: TO REWRITE###########
  #Create zero matrix where:
  #   the number of rows is equal to the number of time steps (in the env data provided)
  #   the number of columns is equal to the number of patches
  env_mat2 <- matrix(0, nrow = max(env_df$time), ncol = M) 
  for (i in 1:max(env_df$time)) {
    env_mat2[i,] <- env_df$env1[env_df$time == i]
  }
  ###############################################################################
  
  burn_in <- time$burn_in[1]
  generations <- max(env_df$time) - burn_in
  Tmax <- burn_in + generations
  
  #set up landscape
  #plot(env_mat2, legend = FALSE)
  
  r <- 5 #rmax (eq3)
  
  #Draw a random value from a standard normal distribution for each species 
  #across all patches
  #Format it into a matrix so that each row is a patch and each column is a species.
  z <- t(replicate(M, rnorm(S))) 
  
  #Create vector of range of dispersal rates with exponential intervals
  a_Vect <- exp(seq(log(1e-5), log(1), length = 16))
  a_Vect <- a_Vect[a_Vect < 1]
  
  #Create vector of range of abiotic niche breadth sigma with exponential intervals
  sigma_niche_V <- exp(seq(log(0.001), log(10), length = 13))
  
  ##############################################################################Not sure yet what this does
  alpha_V <- c(0.5, 1.5) # c(0, 0.5, 1, 1.5) surely 
  
  #Generate matrix of density-dependent competitive effects (alpha_ij) depending on
  #various assumptions (see simulation runs, p1321)
  for (k in 1:(length(alpha_V) + 2)) {
    #alpha <- matrix(rnorm(S^2, mean = 0.07, sd = 0.01), ncol = S, nrow = S) # matrix(rlnorm(S^2, meanlog = -6.9078, sdlog = 2), ncol = S, nrow = S)
    if (k <= length(alpha_V)) {
      if (alpha_V[k] == 0) { 
        #No interaction: alpha_ij == 0
        alpha <- matrix(0, ncol = S, nrow = S)
        alpha_write <- as.character(alpha_V[k])
      } else {
        #Mixed structure: alpha_ijvalues of αij are drawn from a uniform 
        #distribution in the range [0, 1.5], resulting in a combination of 
        #species pairs for which competition is stabilising, destabilising, 
        #or where one of the two is competitively dominant
        alpha <- matrix(runif(S^2, min = 0, max = alpha_V[k]), ncol = S, nrow = S) 
        alpha_write <- as.character(alpha_V[k])
      }
    } else if (k == (length(alpha_V) + 1)) { 
      #Equal structure: alpha_ij == alpha_ii
      alpha <- matrix(1, ncol = S, nrow = S)
      alpha_write <- "equal"
    } else {
      #Competition-colonisation trade-off: values of αij are drawn from a 
      #uniform distribution in the range [0, 1], except for 30% of species, 
      #which are considered dominant species.
      alpha <- matrix(runif(S^2, min = 0, max = 1), 
                      ncol = S, nrow = S)
      alpha_hold <- matrix(runif(S^2, min = 0, max = 1), 
                           ncol = S, nrow = S)
      alpha[1:dominants, ] <- matrix(runif(dominants*S, min = 1, max = 1.5),
                                     ncol = S, nrow = dominants)
      alpha[lower.tri(alpha)] <- alpha_hold[lower.tri(alpha)] #identical(lower.tri(alpha) > 0, lower.tri(alpha)) == 0
      alpha_write <- "patch_dynamics"
    }

    diag(alpha) <- rep(1.0, S) #Always keep alpha_ii at 1 for all species
    alpha <- alpha * 0.05 #scale all competition coefficients by multiplying them by 0.05, to allow for higher equilibrium abundances.
    

    for (i in 1:length(a_Vect)) { #For each level of dispersal
      a <- a_Vect[i] #Get level of dispersal
      
      for (j in 1:length(sigma_niche_V)) { #For each abiotic niche breadth
        sigma_niche <- sigma_niche_V[j] #Get abiotic niche breadth
        
        #M <- 100 #Number of patches already defined earlier on
        
        #initialise the simulation by seeding each habitat patch with populations 
        #of each species drawn from a Poisson distribution where λ = 0.5. 
        #Thus, species start in a random subset of patches and at different abundances.
        N <- matrix(rpois(M*S, 0.5), M, S)
        
        #burn_in <- 2000 #burn_in already defined earlier on
        #Tmax <- 5000 #Tmax already defined earlier on.
        
        seedV <- as.integer(seq(burn_in/10, burn_in/2, length.out = 20))
        sampV <- as.integer(seq(burn_in+800, Tmax, by = 20))
        N_save <- N
        lambda_save <- matrix(0, M, S)
        env_save <- z
        den_save <- matrix(0, M, S)
        env_match_save <- matrix(0, M, S)
        
        for (gen in 1:Tmax) {
          if (any(seedV == gen)) {
            N <- N + rpois(M*S, 0.5) * 1.0
          }
          
          x <- matrix(rep(env_mat2[gen,], S), ncol=S, byrow=T)
          env <- exp(-((x - z) / (2.0 * sigma_niche)) ^ 2.0)
          density <- N * alpha
          lambda_v <- r * N * (1.0 / (1.0 + density)) * env
          lambda_v[lambda_v < 0.0] <- 0.0
          N <- rpois(length(lambda_v), lambda_v)
          
          if (k < (length(alpha_V) + 2)) {
            emigrants <- rbinom(length(N), N, a)
            immigrants_exp <- disp_mat %*% emigrants
            immigrants_S <- apply(emigrants, 2, sum)
            immigrants <- matrix(0, nrow=M, ncol=S)
            for (l in 1:S) {
              immigrants[, l] <- as.numeric(names(sort(table(wsample(1:M, immigrants_exp[, l] / sum(immigrants_exp[, l]), immigrants_S[l]))))) - 1
            }
          } else {
            emigrants <- rbinom(length(N), N, a)
            emigrants[, 1:dominants] <- rbinom(dominants*M, N[, 1:dominants], a*0.1)
            immigrants_exp <- disp_mat %*% emigrants
            immigrants_S <- apply(emigrants, 2, sum)
            immigrants <- matrix(0, nrow=M, ncol=S)
            for (l in 1:S) {
              immigrants[, l] <- as.numeric(names(sort(table(wsample(1:M, immigrants_exp[, l] / sum(immigrants_exp[, l]), immigrants_S[l]))))) - 1
            }
            N[rbinom(M*S, 1, 0.002) > 0] <- 0
          }
          
          emigrants_sum <- apply(emigrants, 2, sum)
          immigrants_sum <- apply(immigrants, 2, sum)
          N <- N - emigrants + immigrants
          N[N < 0.0] <- 0.0
          
          if (any(sampV == gen)) {
            N_save <- array(c(N_save, N), dim=c(dim(N_save)[1], dim(N_save)[2], dim(N_save)[3]+1))
            lambda_save <- array(c(lambda_save, lambda_v), dim=c(dim(lambda_save)[1], dim(lambda_save)[2], dim(lambda_save)[3]+1))
            env_save <- array(c(env_save, x), dim=c(dim(env_save)[1], dim(env_save)[2], dim(env_save)[3]+1))
            env_match_save <- array(c(env_match_save, env), dim=c(dim(env_match_save)[1], dim(env_match_save)[2], dim(env_match_save)[3]+1))
            den_save <- array(c(den_save, density), dim=c(dim(den_save)[1], dim(den_save)[2], dim(den_save)[3]+1))
          }
          N = N*1.0
        }
        
        #To correct:        
        # N_save <- N_save[, , 2:ncol(N_save)]
        # λ_save <- λ_save[, , 2:ncol(λ_save)]
        # env_save <- env_save[, , 2:ncol(env_save)]
        # den_save <- den_save[, , 2:ncol(den_save)]
        # env_match_save <- env_match_save[, , 2:ncol(env_match_save)]
        
        # Model_df_1 <- data.frame(
        #   N = as.vector(N_save),
        #   lambda = as.vector(λ_save),
        #   density = as.vector(den_save),
        #   env_match = as.vector(env_match_save),
        #   env = as.vector(env_save),
        #   Species = rep(1:S, each = M, length.out = length(sampV) * S),
        #   Time = rep(1:length(sampV), each = S * M),
        #   Patch = rep(1:M, length.out = length(sampV) * S * M),
        #   z = matrix(rep(z[1, ], each = M, length.out = length(sampV) * M), nrow = length(sampV) * M),
        #   dispersal = a,
        #   sig_niche = sigma_niche,
        #   alpha = toString(alpha_write)
        # )
        
        # Model_df_1 <- Model_df_1[Model_df_1$N > 0, ]
        # 
        # if (a == a_Vect[1] && sigma_niche == sigma_niche_V[1] && k == 1) {
        #   Model_df <- Model_df_1
        # } else {
        #   Model_df <- rbind(Model_df, Model_df_1)
        # }
      }
    }
  }
}




