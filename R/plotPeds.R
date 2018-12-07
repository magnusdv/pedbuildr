#' @importFrom graphics par plot text title
#' @export
plotPeds = function(pedlist,  titles = NULL, nrow = NA, ...) {

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
    origs = labs[!startsWith(labs, "p")]

    labs[!labs %in% origs] = ""
    mar = if(is.null(tit)) c(1.5,1.5,1.5,1.5) else c(1.5,1.5,3,1.5)

    plot(ped, shaded = origs, col = list(red=origs),
         margin = mar, id.lab = labs, title = tit, ...)
  }
}

#' @export
plotBestPeds = function(x, top = 6, ...) {
  stopifnot2(is.numeric(top), length(top) == 1, top > 0)
  if(!isTRUE(all(c("pedlist", "logliks", "alleleMatrix") %in% names(x)))) {
    stop2("`x` is not a proper output of `reconstruct()`")
  }

  # Sort
  ord = order(x$logliks, decreasing = T)
  logliks = x$logliks[ord]
  pedlist = x$pedlist[ord]

  if(top > length(pedlist))
    top = length(pedlist)

  titles = paste("Loglik =", round(logliks, 2))

  plotPeds(pedlist[1:top], titles = titles[1:top], ...)
}
