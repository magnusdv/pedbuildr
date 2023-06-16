#' @importFrom graphics par plot text title
plotPeds = function(pedlist, labs = NA, highlight = NA, titles = NULL, nrow = NA, ...) {

  labsAll = unlist(labels(pedlist[[1]]))

  if(identical(labs, NA)) {
    labs = labsAll[!grepl("^e[1-9]", labsAll)]
  }

  if(identical(highlight, NA)) {
    highlight = labsAll[!grepl("^e[1-9]", labsAll)]
  }

  hatched = highlight
  col = list(red = highlight)

  L = length(pedlist)
  if(is.na(nrow))
    nrow = if(L<6) 1 else floor(sqrt(L))
  ncol = ceiling(length(pedlist)/nrow)

  op = par(mfrow = c(nrow, ncol))
  on.exit(par(op))

  for(i in seq_along(pedlist)) {
    ped = pedlist[[i]]
    if(is.pedList(ped))
      ped = unclass(ped)

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    mar = 1.5

    plot(ped, labs = labs, hatched = hatched, col = col,
         margin = mar, title = tit, ...)
  }
}
