library(igraph)
library(raster)

#' Run the Cantal model.
#'
#' TODO longer description.
#' The model is stochastic and this simulator uses pseudo-random number generation
#' for initializing the drying location if the random scenario is choosen and for
#' the dynamics of mortality and reproduction. Please consider correctly initializing
#' the RNG seed for garantying reproductibility and for exploring the variability of
#' the model.
#'
#' The given river network is represented by an igraph object with expected attributes
#' X, Y and width for the geographic position of each node and the width of the area.
#'
#' @param graph a river network represented by an igraph object, see below
#' @param S number of species
#' @param iK Carrying capacity in local patches
#' @param d Mortality and reproduction rate
#' @param m Proportion of propagules that disperse
#' @param dm dispersal mode: "drifting", "swimming", "flying"
#' @param d_max Maximum dispersal distance for flying organisms (in meters)
#' @param LDD Long distance dispersal
#' @param int drying intensity (fraction of patches that dry-up)
#' @param dryingmonths drying duration in months
#' @param dryingscenario scenario of drying location: "random" (random location), "upstream" (upstream patches dry first), "downstream" (downstream patches dry first), or a vector of patches indexes (int parameter will be ignored)
#' @param burnsteps duration of the burning simulation phase (weeks)
#' @param nbyears number of simulation steps (years)
#' @importFrom stats rbinom rmultinom rpois
#' @importFrom raster pointDistance
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph neighbors
#' @importFrom igraph vertex_attr
#' @importFrom igraph V
#' @return a list with:
#'  \item{burnin}{gives for each timestep the species richness for each patch during the burn-in phase}
#'  \item{metacom}{the individuals amount for each patch (row), for each species (column) at each timestep of the simulation phase}
#'  \item{simulation}{gives for each timestep the species richness for each patch during the simulation}
#'  \item{drying}{the list of patches concerned by the drying scenario}
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
#' graph <- OCN2graph(OCN)
#' cantal_model(graph=graph, S=10, iK=500, d=0.1, m=0.25, dm="drifting", d_max = 2000,
#'    LDD=0.0001, int=0.8, dryingmonths=6, dryingscenario="upstream", burnsteps=10,
#'    nbyears=1)
#'
#' @export
cantal_model <- function(graph, S, iK, d=0.1, m=0.25, dm, d_max = 2000, LDD=0.001, int=0, dryingmonths, dryingscenario, burnsteps, nbyears) {
  if (is.null(vertex_attr(graph,"X"))) stop("Missing X vector in graph object.")
  if (is.null(vertex_attr(graph,"Y"))) stop("Missing Y vector in graph object.")
  if (is.null(vertex_attr(graph,"width"))) stop("Missing width vector in graph object.")
  dryingweeks <- round((dryingmonths*365/12)/7) # Conversion from months to weeks
  regional_pool <- rep(1/S,S) # Species frequency in the regional pool
  #########################################
  # Compute flow-directed, watercourse and overland distances between patches of the river network
  AD <- as.matrix(as_adjacency_matrix(graph)) # Adjacency matrix of the river network (0/1)
  n_patches <- length(V(graph)) # Number of patches in the river network
  ED <- pointDistance(cbind(V(graph)$X,V(graph)$Y), lonlat=F) # Euclidian overland distances between patches of a dimX*dimY habitat matrix
  diag(ED)=NA
  graph2 = graph
  if(dm=="swimming") {
    D1 = AD+t(AD)
    graph2 = graph_from_adjacency_matrix(D1,mode="directed") # Watercourse distances between patches
  } else if(dm=="drifting") {
    D0 = AD
    graph2 = graph_from_adjacency_matrix(D0,mode="directed") # Flow-directed distances between patches
  }
  #########################################
  # Scenarios of drying location
  if (is.numeric(dryingscenario)) {
    d_patches = dryingscenario
  } else {
    width = sqrt(V(graph)$width)
    if (dryingscenario == "random") {
      d_patches = sample(c(1:n_patches),size=(n_patches*int)) # Random location of drying events
    } else if (dryingscenario == "upstream") {
      d_patches = order(width,decreasing=F)[1:(n_patches*int)] # Upstream location of drying events
    } else if (dryingscenario == "downstream") {
      d_patches = order(width,decreasing=T)[1:(n_patches*int)] # Downstream location of drying events
    }
  }
  #########################################
  # Simulation
  metacom = t(rmultinom(n_patches,size=iK,prob=regional_pool))
  deadcom = matrix(0,nrow=n_patches,ncol=S)
  deadcoms=NULL
  outS = NULL
  outS2 = NULL
  out_metacom = list()
  out_metacom = NULL
  steps = c(-1*(burnsteps-seq(burnsteps)),seq(nbyears*52))
  for (istep in steps) {
    # Mortality (a fraction d of the individuals)
    deadcom[] <- 0
    for (i in 1:n_patches){
      for (j in 1:S){
        if (metacom[i,j]>0){
          deadcom[i,j] = min(metacom[i,j], rbinom(1,as.numeric(metacom[i,j]),prob=d))
        }
      }
    }
    if (istep<=0) {deadcoms=c(deadcoms,deadcom)}
    # Update mortality
    metacom = metacom-deadcom
    ####################################################
    # Drying event if not in burning phase and in drying period
    if ((istep>0) && ((istep-1)%%52 < dryingweeks)) {
      metacom[d_patches,] = 0
    }
    ####################################################
    # Production and dispersal of propagules (d -> rate of propagule production, m -> fraction of propagule that disperse)
    propagule_cloud = matrix(0,nrow=n_patches,ncol=S)
    for (i in 1:n_patches){
      if(dm!="flying"){ # Drifting or Swimming organisms
        nodes = which(graph2[i]!=0)
        n_nodes = sum(graph2[i])
        if (length(neighbors(graph2,i))==0) { n_nodes = n_nodes+1 }  # Addition of a neighbour node  at the outlet to avoid n_nodes = 0.
      }
      if(dm=="flying"){ # Flying organisms
        nodes =  which(ED[i,]<=d_max)
        n_nodes = length(nodes)
      }
      for (j in 1:S){
        # Computation of propagule cloud: propagules are produced at a rate d are distributed evenly to the neighbouring nodes (1/n_nodes)
        if (length(nodes)!=0){
          propagule_cloud[nodes,j] = propagule_cloud[nodes,j]+ rpois(1,d*m*(1/n_nodes)*metacom[i,j])
        }
        # Computation of local recruitment: propagules are produced at a rate d and a fraction (1-m) stays in the node
        propagule_cloud[i,j] = propagule_cloud[i,j] + rpois(1,d*(1-m)*metacom[i,j]) + rpois(1,LDD*regional_pool[j])
      }
    }
    # Recruitment
    n_recruits = 0
    for (i in 1:n_patches){
      if (sum(metacom[i,])<iK){
        n_recruits = rpois(1,iK-sum(metacom[i,]))
        prop=propagule_cloud[i,]
        if (sum(prop)!=0){
          propagule_cloud[i,] = c(t(rmultinom(1,size=n_recruits,prob=prop)))
        }
      }
    }
    metacom = metacom + propagule_cloud
    # Compute stats
    pa = matrix(0,nrow=n_patches,ncol=S)
    pa[which(metacom!=0)]=1
    stats=apply(pa,1,sum)
    if (istep<=0) { # burning phase
      outS = rbind(outS,stats)
    } else {
      outS2 = rbind(outS2,stats)
      out_metacom[[istep]] = metacom
    }
  }
  #########################################
  # Return outputs
  # Burn in statistics: species richness only
  rownames(outS) <- paste("t",c(1:burnsteps),sep="")
  # Metacommunity abundance and species richness at each time step
  rownames(outS2) <- paste("t",c(1:(nbyears*52)),sep="")
  list(burnin=outS, metacom=out_metacom, simulation=outS2, drying=d_patches)
}
