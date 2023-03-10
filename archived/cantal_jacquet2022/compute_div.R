library(betapart)
#' Compute diversity indices
#'
#'Compute the average and standard deviation of alpha, beta and gamma diversity over the last steps of the simulation generated by the function drying_model of the package cantal.
#' Alpha diversity corresponds to local species richness in all patches of the network, gamma diversity corresponds to regional species richness and beta diversity corresponds to overall dissimilarity in species composition of local communities.
#' @param res list with  the simulation results generated by the function "drying_model" of the package cantal
#' @param nstep number of simulation steps over which the diversity indices are calculated
#' @param index.family family of dissimilarity indices ("sorensen", "jaccard","bray" or "ruzicka") to compute overall beta diversity, calculated with the function beta.multi (Sorensen and Jaccard dissimilarity) or the function beta.multi.abund (Bray_Curtis or Ruzicka multiple-site dissimilarity) of the package betapart
#' @importFrom betapart beta.multi.abund
#' @importFrom betapart beta.multi
#' @importFrom stats sd
#' @return The function returns a list with the mean and standard deviation of alpha, beta and gamma diversity over the last "nstep" of the simulation:
#' \item{alpha_mean}{alpha diversity, averaged over the last "nstep" steps of the simulation}
#' \item{alpha_sd}{standard deviation of alpha diversity over the last "nstep" steps of the simulation}
#' \item{beat_mean}{beta diversity, averaged over the last "nstep" steps of the simulation}
#' \item{beta_sd}{standard deviation of beta diversity over the last "nstep" steps of the simulation}
#' \item{gamma_mean}{gamma diversity, averaged over the last "nstep" steps of the simulations}
#' \item{gamma_sd}{standard deviation of gamma diversity over the last "nstep" steps of the simulation}
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
#' res = cantal_model(graph=OCN2graph(OCN), S=10, iK=500, d=0.1, m=0.25, dm=0, d_max = 2000,
#'    LDD=0.0001, int=0.8, dryingmonths=6, dryingscenario=1, burnsteps=10,
#'    nbyears=1)
#' compute_div(res)
#'
#' @export
compute_div <- function(res,nstep=20,index.family="sorensen"){
  gamma_t=NULL
  beta_t=NULL
  for(i in c((length(res$metacom)-(nstep-1)):length(res$metacom))){
    gamma_t = c(gamma_t,length(which(apply(res$metacom[[i]],2,sum)!=0))) # Regional species richness (gamma diversity) at time k
    beta_mat = res$metacom[[i]]
    beta_mat[which(beta_mat!=0)]=1
    if(index.family=="jaccard"){
      beta_t = c(beta_t,betapart::beta.multi(beta_mat,index.family = "jaccard")$beta.JAC) # value of the overall beta diversity, measured as Jaccard dissimilarity
    }
    if(index.family=="sorensen"){
      beta_t = c(beta_t,betapart::beta.multi(beta_mat,index.family ="sorensen")$beta.SOR) # value of the overall beta diversity, measured as Sorensen dissimilarity
    }
    if(index.family=="bray"){
      beta_t = c(beta_t,betapart::beta.multi.abund(res$metacom[[i]],index.family = "bray")$beta.BRAY) # value of the overall dissimilarity, measured as Bray-Curtis multiple-site dissimilarity
    }
    if(index.family=="ruzicka"){
      beta_t = c(beta_t,betapart::beta.multi.abund(res$metacom[[i]],index.family = "ruzicka")$beta.RUZ) # value of the overall dissimilarity, measured as Ruzicka multiple-site dissimilarity
    }
  }
  list(
  alpha_mean = apply(res$simulation[(nrow(res$simulation)-(nstep-1)):nrow(res$simulation),],2,mean),
  alpha_sd = apply(res$simulation[(nrow(res$simulation)-(nstep-1)):nrow(res$simulation),],2,sd),
  beta_mean = mean(beta_t),
  beta_sd = sd(beta_t),
  gamma_mean = mean(gamma_t),
  gamma_sd = sd(gamma_t)
   )
}
