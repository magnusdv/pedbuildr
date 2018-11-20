#' @importFrom graphics par plot text title
#' @export
plotPeds = function(pedlist,  titles = NULL,
                    nrow = floor(sqrt(length(pedlist))),
                    ncol = ceiling(length(pedlist)/nrow), ...) {

  op = par(mfrow = c(nrow, ncol)); on.exit(par(op))
  for(i in seq_along(pedlist)) {
    ped = pedlist[[i]]

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    if(is.pedList(ped)) {
      plot(-1:1,-1:1, type="n", axes = F, xlab="", ylab="")
      mess ="Disconnected pedigree.\n\nPlot separately with\n`pedtools::plotPedList()`."
      text(0, 0, mess, cex = 1.3)
      title(tit)
      next
    }

    labs = labels(ped)
    origs = labs[!startsWith(labs, "p")]

    labs[!labs %in% origs] = ""
    mar = if(is.null(tit)) c(1.5,1.5,1.5,1.5) else c(1.5,1.5,3,1.5)

    plot(ped, shaded = origs, col = list(red=origs),
         margin = mar, id.lab = labs, title = tit, ...)
  }
}


