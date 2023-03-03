library(igraph)
library(raster)
#' Plot perceived graph
#'
#' Plotting of a spatial network as perceived by an organism with a given dispersal mode and maximum dispersal distance.
#'  Three modes of dispersal commonly observed in riverine communities are available: (i) drifting, with downstream flow-directed dispersal along the waterway (e.g., plant seeds and planktonic larvae), (ii) swimming, with bidirectional dispersal along the waterway (e.g., fishes), and (iii) flying, with overland dispersal (insects with flying adult stages).
#'
#' For drifting dispersal, individuals dispersed to the patch situated downstream only.
#'
#' For swimming dispersal, individuals dispersed to patches that are connected by water flow, which could be located upstream or downstream.
#'
#' For flying dispersal, individuals dispersed to the patches situated closer than a threshold overland distance set by the user (dm = 2km by default, corresponding to a maximum flying distance that has been reported in the literature for mayflies (Ephemeroptera) (Kovats et al. 1996)).
#'
#' @param OCN an optimal channel network (OCN) created with the R package OCNet
#' @param dm dispersal mode: "drifting", "swimming" or "flying"
#' @param d_max maximum dispersal distance for flying organisms (in meters)
#' @param vertex.size graphical parameter of the plot.igraph function
#' @param vertex.color graphical parameter of the plot.igraph function
#' @param vertex.label graphical parameter of the plot.igraph function
#' @importFrom raster pointDistance
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @importFrom igraph simplify
#' @importFrom igraph igraph_opt
#' @return Returns NULL, invisibly.
#'
#' @seealso
#' igraph.plotting for the detailed description of the plotting parameters.
#' @examples
#' library(OCNet)
#' set.seed(1)
#' dimX = 25
#' dimY = 25
#' cellsize = 500
#' thrA = 5*cellsize^2
#' OCN <- create_OCN(dimX, dimY, cellsize = cellsize)
#' OCN <- landscape_OCN(OCN)
#' OCN <- aggregate_OCN(OCN, thrA = thrA)
#' par(mar=rep(0,4))
#' draw_perceived_graph(OCN,dm="drifting",d_max=1000)
#' @export
draw_perceived_graph <- function(OCN,dm,d_max,vertex.size=5,vertex.color="#309BCD",vertex.label=NA){
  if(dm=="drifting") {
    graph = igraph::graph_from_adjacency_matrix(
      as.matrix(OCN$RN$W),
      mode="undirected") # Flow-directed distances between patches
  }
  if(dm=="swimming") {
    graph = igraph::graph_from_adjacency_matrix(
      (as.matrix(OCN$RN$W) + t(as.matrix(OCN$RN$W))),
      mode="undirected") # Watercourse distances between patches
  }
  if(dm=="flying") {
    ED <- raster::pointDistance(expand.grid(c(1:OCN$dimX),
                                            c(1:OCN$dimY)),
                                lonlat=F) # Euclidian overland distances between patches of a dimX*dimY habitat matrix
    ED_OCN <- ED[which(OCN$FD$toRN!=0),
                 which(OCN$FD$toRN!=0)]*OCN$cellsize # Euclidian overland distances of the OCN, multiplied by OCN cellsize
    diag(ED_OCN)=0
    ED_OCN[which(ED_OCN>(d_max))]=NaN # Dispersal distance determines which nodes are connected
    
    graph=igraph::graph_from_adjacency_matrix(ED_OCN,mode="undirected")
    graph=igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, 
                           edge.attr.comb = igraph_opt("edge.attr.comb"))
  }
  igraph::plot.igraph(graph, 
                      layout = matrix(c(OCN$RN$X,OCN$RN$Y), 
                                     ncol = 2, 
                                     nrow = OCN$RN$nNodes),
                      vertex.size = vertex.size,
                      vertex.color = vertex.color,
                      vertex.label=NA)
}
