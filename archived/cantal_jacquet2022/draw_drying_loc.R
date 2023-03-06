library(OCNet)
#' Plotting location of drying events in a river network.
#'
#' @param OCN an optimal channel network (OCN) created with the R package OCNet.
#' @param drying_location location of drying events: a numeric vector with the status of each patch of the network (1 if dry, 0 otherwise) or a character string ("random","upstream" or"downstream") specifying the location of drying events. In the latter case, the parameter drying_intensity must be filled.
#' @param drying_intensity fraction of patches affected by drying events. Required if drying location is not a numeric vector.
#' @param colPalette colors displayed for affected and unaffected patches.
#' @param backgroundColor color used in the background of the figure.
#' @param cex a numerical value giving patch size (equivalent to parameter cex of function points).
#' @importFrom OCNet draw_thematic_OCN
#' @importFrom graphics title
#' @return Returns NULL, invisibly.
#'
#' @seealso
#' Function draw_thematic_OCN of the package OCNet for the detailed description of the plotting parameters.
#' @examples
#' library(OCNet)
#' set.seed(1)
#' dimX = 20
#' dimY = 20
#' cellsize = 500
#' thrA = 5*cellsize^2
#' OCN <- create_OCN(dimX, dimY, cellsize = cellsize)
#' OCN <- landscape_OCN(OCN)
#' OCN <- aggregate_OCN(OCN, thrA = thrA)
#' draw_drying_loc(OCN,drying_location="random",drying_intensity=0.1)
#'
#' @export
draw_drying_loc <- function(OCN,drying_location,drying_intensity=0,colPalette=c("#309BCD","#EFA528"),backgroundColor = NULL,cex = 1.4)
{
  if(is.numeric(drying_location)){
    patch_status=drying_location
    if(sum(patch_status)==OCN$RN$nNodes){colPalette=colPalette[2]}
    }
else{
  if(drying_location=='random') {
    patch_status = rbinom(OCN$RN$nNodes,1,prob=drying_intensity)# Random location of dry patches: 0/1 is wet/dry
    }
  if(drying_location=='upstream') {
    patch_status = numeric(OCN$RN$nNodes)
    if(drying_intensity!=0){
    patch_status[order(OCN$RN$A,decreasing=F)[1:(OCN$RN$nNodes*drying_intensity)]]=1  # Upstream location of dry patches:
    }
  }
  if(drying_location=='downstream') {
    patch_status = numeric(OCN$RN$nNodes)
    if(drying_intensity!=0){
    patch_status[order(OCN$RN$A,decreasing=T)[1:(OCN$RN$nNodes*drying_intensity)]]=1  # Downstream location of dry patches:
    }
  }
  if(drying_intensity==1){
    colPalette=colPalette[2]
    }
}
  OCNet::draw_thematic_OCN(patch_status,OCN,drawNodes=T,backgroundColor = backgroundColor,cex = 1.4, discreteLevels = T, colPalette=colPalette,addLegend = T,colLevels = NULL)
  drying_int=round(length(which(patch_status!=0))/length(patch_status),2)
  title(main=list(paste("Drying intensity =",drying_int),cex=1.5))
  }
