library(igraph)
library(sf)

#' Converts an OCN object to an igraph object.
#'
#' This helper function converts an optimal channel network (OCN) created with
#' the R package OCNet to an igraph object, ready to use with the Cantal model.
#'
#' @param OCN an optimal channel network (OCN) created with the R package OCNet.
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph set_vertex_attr
#' @return a igraph object.
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
#'
#' @export
OCN2graph <- function(OCN) {
  graph = graph_from_adjacency_matrix(as.matrix(OCN$RN$W))
  # Put geographic position
  nodes_coords = expand.grid(c(1:OCN$dimX),c(1:OCN$dimY))[which(OCN$FD$toRN!=0),] # selon code Claire
  graph = set_vertex_attr(graph, "X", value = nodes_coords[,1]*OCN$cellsize)
  graph = set_vertex_attr(graph, "Y", value = nodes_coords[,2]*OCN$cellsize)
  # Put width
  graph = set_vertex_attr(graph, "width", value = OCN$RN$A)
  graph
}

#' Converts a Hydrorivers shapefile to an igraph object.
#'
#' This helper function converts a river network stored in a shapefile from
#' [Hydrorivers](https://www.hydrosheds.org/page/hydrorivers)
#' to an igraph object, ready to use with the Cantal model.
#'
#' @param filename The filename of the shapefile to load.
#' @param MAIN_RIV_ID The ID of the main river of the watershed to load.
#'
#' @importFrom sf st_read
#' @return a igraph object.
#'
#' @examples
#' \dontrun{
#' graph <- hydroriversSHP_to_graph("HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu.shp", "20487373")
#' }
#' @export
hydroriversSHP_to_graph  <- function(filename, MAIN_RIV_ID) {
  shp = st_read(filename,  query = paste0("SELECT HYRIV_ID,UPLAND_SKM,NEXT_DOWN FROM \"HydroRIVERS_v10_eu\" WHERE MAIN_RIV=\'",MAIN_RIV_ID,"\'"))
  hydroriversSF_to_graph(shp)
}

#' Converts a Hydrorivers data loaded in SimpleFeature dataset to an igraph object.
#'
#' This helper function converts a river network loaded from
#' [Hydrorivers](https://www.hydrosheds.org/page/hydrorivers)
#' to an igraph object, ready to use with the Cantal model.
#'
#' @param sfdata The SimpleFeature dataset (as return by sf::st_read)
#'
#' @importFrom sf st_cast
#' @importFrom sf st_sf
#' @importFrom sf st_sfc
#' @importFrom sf st_coordinates
#' @importFrom sf st_centroid
#' @importFrom igraph graph_from_data_frame
#' @return a igraph object.
#'
#' @examples
#' \dontrun{
#'  sfdata = st_read(filename,  query = paste0("SELECT HYRIV_ID,UPLAND_SKM,NEXT_DOWN FROM \"HydroRIVERS_v10_eu\" WHERE MAIN_RIV=\'",MAIN_RIV_ID,"\'"))
#'  graph <- hydroriversSF_to_graph(sfdata)
#' }
#' @export
hydroriversSF_to_graph  <- function(sfdata) {
  # "Simple Feature" Dataframe to store coordinates of two bounds of each network reach
  bounds = data.frame(matrix(ncol = 2, nrow = 0))
  colnames(bounds) = c('name','geometry')
  # "Simple Feature" Dataframe to store centroids of the reaches of the network
  vertices <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(vertices) = c('name', 'X', 'Y', 'width', 'outlet')
  # Loop on reaches
  for(i in seq_len(nrow(sfdata))) {
    # reduce the polylines to simple segments
    points_i = st_cast(sfdata[i,]$geometry, "POINT")
    start = points_i[1]
    end = points_i[length(points_i)]
    bounds = rbind(bounds, st_sf(name=sfdata[i,]$HYRIV_ID, geometry = st_sfc(start)))
    bounds = rbind(bounds, st_sf(name=sfdata[i,]$HYRIV_ID, geometry = st_sfc(end)))
    centroid = st_coordinates(st_centroid(sfdata[i,]$geometry))
    vertices[i,]=c(sfdata[i,]$HYRIV_ID, centroid[1], centroid[2], sfdata[i,]$UPLAND_SKM, sfdata[i,]$NEXT_DOWN==0)
  }
  # Now generating the edges
  outlet = sfdata[sfdata$NEXT_DOWN ==0,]$HYRIV_ID
  remainingReaches = bounds
  currentReaches = c(outlet)
  edges <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(edges) = c('from','to')
  while (length(currentReaches)>0){
    nextReaches = c()
    for (reach in currentReaches) {
      reachcoo = remainingReaches[remainingReaches$name==reach,]
      # this reach will not more be considered for next steps
      remainingReaches = remainingReaches[remainingReaches$name!=reach,]
      # processing reaches that share one point with current
      neighbours = remainingReaches[remainingReaches$geometry %in% reachcoo$geometry,]
      for (i in seq_len(nrow(neighbours))) {
        nextReaches = c(nextReaches, neighbours[i,]$name)
        edges = rbind(edges, c(neighbours[i,]$name, reach))
      }
      # removing processed confluence points
      remainingReaches = remainingReaches[!(remainingReaches$geometry %in% reachcoo$geometry),]
    }
    currentReaches = nextReaches
  }
  # Returns the resulting graph
  graph_from_data_frame(edges, directed=TRUE, vertices=vertices)
}
