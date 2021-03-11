#' @importFrom graphics par plot text title
plotPeds = function(pedlist, titles = NULL, nrow = NA, ...) {

  L = length(pedlist)
  if(is.na(nrow))
    nrow = if(L<6) 1 else floor(sqrt(L))
  ncol = ceiling(length(pedlist)/nrow)

  op = par(mfrow = c(nrow, ncol))
  on.exit(par(op))

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
    origs = labs[!grepl("^e[1-9]", labs)]

    mar = if(is.null(tit)) c(1.5,1.5,1.5,1.5) else c(1.5,1.5,3,1.5)

    plot(ped, hatched = origs, col = list(red = origs),
         margin = mar, labs = origs, title = tit, ...)
  }
}
