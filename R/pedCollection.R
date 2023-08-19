#' @export
`[.pedCollection` = function(x, i) {
  y = unclass(x)[i]
  class(y) = "pedCollection"
  y
}


#' @export
print.pedCollection = function(x, ...) {
  n = length(x)
  isCon = vapply(x, is.ped, FUN.VALUE = TRUE)

  cat(sprintf("List of %d pedigree%s, ", n, if(n==1) "" else "s"))
  if(n == 1)
    cat(if(isCon[1]) "connected" else "disconnected")
  else if(all(isCon))
    cat("all connected")
  else
    cat(sprintf("%d connected, %d disconnected", sum(isCon), n - sum(isCon)))

  cat("\n")
}




#' @importFrom graphics box par plot text title
#' @export
plot.pedCollection = function(x, labs = NA, highlight = NA, titles = NULL,
                              nrow = NA, frames = FALSE, ...) {
  pedlist = x
  labsAll = unname(unlist(labels(pedlist[[1]])))

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
  on.exit(par(op), add = TRUE)

  for(i in seq_along(pedlist)) {
    ped = pedlist[[i]]
    if(is.pedList(ped))
      ped = unclass(ped)

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    mar = 1.5

    plot(ped, labs = labs, hatched = hatched, col = col,
         margin = mar, title = tit, ...)
    if(frames)
      box("figure")
  }
}
