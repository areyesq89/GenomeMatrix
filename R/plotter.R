#' @title Plot upper half of a matrix rotated -45 degrees.
#' @aliases GenomeMatrix
#' @description
#' This function inputs a matrix along with an optional \code{\link{GenomicRanges}}
#' object and plots the upper half of the matrix rotated -45 degrees. 
#'
#' @param mat A n x n symmetric matrix to plot.
#' @param granges An optional genomic ranges object of length n representing the genomic coordinates the 
#' rows and columns of the matrix.
#' @param plotGR A \code{\link{GenomicRanges}} object indicating the genomic region to be plotted. 
#' This option is ignored if granges is not provided.
#' @param extend Proportion of the plotting region to extend both upstream and downstream. These 
#' extra regions are normally be plotted, but if white spaces appear in the upper corners of the 
#' plot, consider increasing the value of this parameter. 
#' @param heightProp Proportion of the plotting region (x-axis) to plot in the y-axis.
#' @param zlim Limit for the maximum value of z (color) to plot. Larger values will be truncated
#' to this parameter.
#' @param colorBias Bias parameter for color palette (see ?colorRampPalette).
#' @param highlight Optional. A \code{\link{GenomicRanges}} object indicating a region to highlight. 
#' If specified, A box would be drawn around this region.
#' @param colRamp A name of a color palette (see ?brewer.pal).
#' @param smoothFilt An optional smoothing function for the matrix. 
#'
#' @examples
#' library(HiTC)
#' data(Dixon2012_IMR90, package="HiCDataHumanIMR90")
#' mat <- as.matrix( hic_imr90_40@.Data[[1]]@intdata )
#' granges <- hic_imr90_40@.Data[[1]]@xgi
#' pl <- matrixPlotter( log2(mat+1), granges, plotGR=GRanges("chr1", IRanges( 50000000, 60000000 ) ), zlim=5 ) 
#' \dontrun{
#' print(pl)
#' }
#'
#' @return A ggplot2 object.
#' @author Alejandro Reyes
#' @export
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
matrixPlotter <- function( mat, granges=NULL, plotGR=NULL, extend=0.5, heightProp=1/3,
                            zlim=NULL, colorBias=1, highlight=NULL, colRamp="YlOrRd",
                            smoothFilt=NULL ){
  stopifnot( is.matrix( mat ) )
  stopifnot( dim( mat )[1] == dim( mat )[2] )
  if( is.null( granges ) ){
    granges <- GRanges("Z",
                       IRanges(
                         start=seq_len( nrow(mat) ),
                         end=seq_len( nrow(mat) ) + 1 ) )
    res <- 1
    plotGR <- range( granges )
    extend <- 0
  }
  stopifnot( is( granges, "GenomicRanges" ) )
  stopifnot( dim( mat )[1] == length( granges ) )
  stopifnot( length( unique( width( granges )[dim(mat)-1] ) ) == 1 )
  res <- unique( width( granges )[dim(mat)-1] ) - 1
  if( !is.null( smoothFilt ) ){
    mat <- smoothFilt( mat )
  }
  originalWidth <- width( plotGR )
  extend <- width( plotGR ) * extend
  xstart <- start( plotGR )
  xend <- end( plotGR )
  matLong <- melt( mat )
  matLong <- matLong[matLong$Var2 >= matLong$Var1,]
  binSize <- res
  binNum <- (xend - xstart)/binSize
  start(plotGR) <- start(plotGR) - extend
  end(plotGR) <- end(plotGR) + extend
  bins <- subjectHits( findOverlaps( plotGR, granges ) )
  if( !length( bins ) > 0 ){
    stop("No overlaps found\n")
  }
  #matLongSub <- dplyr::filter( matLong, Var1 %in% bins, Var2 %in% bins )
  matLongSub <- matLong[matLong$Var1 %in% bins & matLong$Var2 %in% bins,]
  matLongSub$binID <- sprintf( "id%d", seq_len( nrow( matLongSub ) ) )
  binDiffs <- matLongSub$Var2 - matLongSub$Var1
  centers <- rowMeans(
    data.frame(
      start(granges)[matLongSub$Var1],
      end(granges)[matLongSub$Var2] ) )
  maxys <- 1 + ( binDiffs - 1 )*0.5
  minys <- maxys - 1
  coordsDF <- rbind(
    data.frame(
      x=centers - (binSize/2),
      y=rowMeans( data.frame(maxys, minys) ),
      val=matLongSub$value,
      bin=matLongSub$binID ),
    data.frame(
      x=centers,
      y=maxys,
      val=matLongSub$value,
      bin=matLongSub$binID ),
    data.frame(
      x=centers + (binSize/2),
      y=rowMeans( data.frame(maxys, minys) ),
      val=matLongSub$value,
      bin=matLongSub$binID ),
    data.frame(
      x=centers,
      y=minys,
      val=matLongSub$value,
      bin=matLongSub$binID ) )
  coordsDF$x <- coordsDF$x - binSize/2
  xstart <- xstart - binSize/2
  xend <- xend - binSize/2
  coordsDF$y <- pmax( 0, coordsDF$y )
  coordsDF$y <- coordsDF$y * binSize
  ylim2 <- originalWidth*heightProp
  cols <- colorRampPalette( brewer.pal( 9, colRamp ), bias=colorBias )( 50 )
  if( !is.null( zlim ) ){
    coordsDF$val <- pmax( -zlim, pmin( coordsDF$val, zlim ) )
  }
  p <- ggplot(coordsDF, aes(x = x, y = y)) +
    geom_polygon(aes(fill = val, group = bin)) +
    ylab("") +
    scale_fill_gradientn( colours=cols ) +
    coord_fixed( ratio=1, xlim=c(xstart, xend), ylim=c(0, ylim2) ) +
    theme( axis.text.y=element_blank(), legend.position="none" )
  if( !is.null( highlight ) ){
    highlight <- subsetByOverlaps( highlight, plotGR )
    start(highlight) <- start( highlight ) - binSize/2
    end(highlight) <- end (highlight) - binSize/2
    p <- p + geom_rect(
      data = as.data.frame( highlight ), inherit.aes=FALSE,
      aes( xmin=start, xmax=end, ymin=-Inf, ymax=Inf ),
      colour="black", fill=NA, alpha=0.5 )
  }
  p
}
