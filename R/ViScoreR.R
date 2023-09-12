
#' Plot 2-d embedding of single-cell data
#'
#' Plots the embedding result of a dimension-reduction algorithm, with cell labels.
#' If any 'unassigned' labels are specified, the unassigned points are plotted in light grey in the background.
#'
#' @param X numeric matrix: row-wise coordinate matrix of 2-dimensional embedding
#' @param annot string/factor vector: vector of labels per row of \code{X}
#' @param unassigned string (vector): labels in `annot` that belong to cells of unknown identity/debris/doublets (default: none)
#' @param palette string vector: palette of colours per population to use for plotting (default palette is provided)
#' @param pointsize numeric: size of each point in the plot (default: 0.5)
#' @param no_legend logical: disable displaying the legend mapping colours to populations? (default: FALSE)
#' @param font_family optional string: non-default font family to use for plotting
#'
#' @export
PlotEmbedding <- function(
    X,
    annot,
    unassigned = c(),
    palette = .Palette(),
    pointsize = 0.5,
    no_legend = FALSE,
    font_family = NULL
) {
  
  ## Check arguments ----
  if (!is.matrix(X) || !is.numeric(X)) {
    stop('Invalid `X` value')
  }
  if ((!is.character(annot) && !is.factor(annot))) {
    stop('Invalid `annot` value')
  }
  if (length(annot) != nrow(X)) {
    stop('Length of `annot` must match number of rows of `X`')
  }
  if (length(unassigned) > 0 && !is.character(unassigned)) {
    stop('Invalid `unassigned` value')
  }
  au <- unassigned[!unassigned %in% levels(annot)]
  if (length(au) > 0) {
    warning('The following `unassigned` populations are not present in `annot`:\n\t', paste0(au, collapse = ',\n\t'))
  }
  if (!is.character(palette) || !is.vector(palette)) {
    stop('Invalid `palette` value')
  }
  if (!is.numeric(pointsize) || pointsize <= 0) {
    stop('Invalid `pointsize`: ', pointsize)
  }
  if (!is.logical(no_legend)) {
    stop('Invalid `no_legend` value')
  }
  if (!is.null(font_family) && (!is.character(font_family) || length(font_family) != 1)) {
    stop('Invalid `font_family` value') 
  }
  annot <- as.factor(annot)
  pops <- levels(annot)
  pops <- pops[!pops %in% exclude]
  lpa <- length(palette)
  lpo <- length(pops)
  if (lpa < lpo) {
    stop('Provided `palette` only contains ', lpa, ' colours, but there are ', lpo, ' labelled populations')
  }
  
  ## Identify and separate unassinged cells ----
  u <- annot %in% unassigned
  d0 <- data.frame()
  if (any(u)) {
    d0 <- data.frame(X[u, , drop = FALSE], 'Population' = as.factor('unassigned'))
  }
  d1 <- data.frame(X[!u, , drop = FALSE], 'Population' = annot[!u])
  
  ## Prep plotting data ----
  colnames(d0)[1:2] <- c('X', 'Y')
  colnames(d1)[1:2] <- c('X', 'Y')
  
  ## Create plot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = d0, colour = 'grey', mapping = ggplot2::aes(x = X, y = Y), size = pointsize) +
    ggplot2::geom_point(data = d1, mapping = ggplot2::aes(x = X, y = Y, col = Population), size = pointsize) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
  
  ## Resolve legend and font ----
  if (no_legend) {
    p <- p + ggplot2::theme(legend.position = 'none')
  }
  if (!is.null(font_family)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }
  
  p
}

#' Plot likeness distributions
#'
#' Plots the likeness distributions (from xNPE) for a specific population in original data and an embedding.
#'
#' @param xnpe xNPE object, result of the \code{ViScoreR::xNPE} function call
#' @param pop string: name of population
#' @param bin_width integer: width of bins (like-neighbour count values from 0 to k+1) for aggregating the data (default: 10)
#' @param font_family optional string: font family to use for plotting
#'
#' @export
PlotLD <- function(
  xnpe,
  pop,
  bin_width = 10,
  font_family = NULL
) {
  
  ## Check arguments ----
  if (class(xnpe) != 'xNPE') {
    stop('Invalid `xnpe` value')
  }
  if (!is.character(pop) || length(pop) != 1) {
    stop('Invalid `pop` value')
  }
  if (!pop %in% names(attributes(xnpe)$distr_hd)) {
    stop('Population "', pop, '" not found in `xnpe`')
  }
  if (!.IsInt(bin_width) || bin_width < 1) {
    stop('Invalid `bin_width` value')
  }
  if (!is.null(font_family) && (!is.character(font_family) || length(font_family) != 1)) {
    stop('Invalid `font_family` value') 
  }
  
  ## Extract likeness distributions ----
  dhd <- attributes(xnpe)$distr_hd[[pop]]
  dld <- attributes(xnpe)$distr_ld[[pop]]
  
  ## Resolve binning ----
  bin <- function(x, bin_width = 1) {
    if (bin_width == 1) {
      return(x)
    }
    n <- length(x)
    
    i1 <- seq(from = 1, to = n, by = bin_width)
    i2 <- c(seq(from = bin_width, to = n, by = bin_width), n)
    if (length(i2) > length(i1)) {
      i2 <- i2[1:(length(i2)-1)]
    }
    mapply(function(i,j) sum(x[i:j]), i1, i2)
  }
  dhd <- bin(dhd, bin_width)
  dld <- bin(dld, bin_width)
  
  ## Prep plotting data ----
  d <- rbind(
    data.frame('Value' = bin_width*seq_along(dhd), 'Count' = dhd, 'Data' = as.factor('Original')),
    data.frame('Value' = bin_width*seq_along(dld), 'Count' = dld, 'Data' = as.factor('Embedding'))
  )
  
  ## Create plot ----
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Value, y = Count, fill = Data)) +
    ggplot2::geom_density(stat = 'identity', alpha = 0.5, linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c('lightblue', 'pink')) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste0(pop, ' likeness distributions')) +
    ggplot2::xlab('Like-neighbour count') +
    ggplot2::ylab('Frequency') +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  ## Resolve font ----
  if (!is.null(font_family)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }
  
  p
}

#' Plot xNPE per population
#'
#' Plots one or more xNPE results for a single dataset.
#' xNPE values for each population are shown, and can be compared across multiple embedding methods.
#'
#' @param ... one or more xNPE objects, results of \code{ViScoreR::xNPE} function calls
#' @param methods string vector: name of each xNPE result provided
#' @param exclude string vector: labelled populations to exclude (default: empty)
#' @param polar logical: use polar coordinates for plotting? (default: FALSE)
#' @param palette string vector: palette of colours per population to use for plotting (default palette is provided)
#' @param font_family optional string: non-default font family to use for plotting
#'
#' @export
PlotxNPE <- function(
  ...,
  methods = c('default'),
  exclude = c(),
  polar = FALSE,
  palette = .Palette(),
  font_family = NULL
) {
  
  args <- list(...)
  
  ## Check arguments ----
  if (any(sapply(args, class) != 'xNPE')) {
    stop('At least one of the provided xNPE objects is not valid')
  }
  la <- length(args)
  if (la > 1) {
    pops <- lapply(seq_len(la), function(idx) names(attributes(args[[idx]])$distr_hd))
    for (idx in 1:(la-1)) {
      if (any(pops[[idx]] != pops[[idx+1]])) {
        stop('Population names found in xNPE objects are inconsistent')
      }
    }
    pops <- pops[[1]]
  } else {
    pops <- names(attributes(args[[1]])$distr_hd)
  }
  if (!is.character(methods)) {
    stop('Invalid `methods` value')
  }
  if (length(methods) != length(args)) {
    stop('`methods` vector must be as long as the number of xNPE vectors')
  }
  if (length(exclude) > 0 && !is.character(exclude)) {
    stop('Invalid `exclude` value')
  }
  au <- exclude[!exclude %in% pops]
  if (length(au) > 0) {
    warning('The following `exclude` populations were not found in an xNPE object:\n\t', paste0(au, collapse = ',\n\t'))
  }
  if (!is.logical(polar) || length(polar) != 1) {
    stop('Invalid `polar` value')
  }
  if (!is.null(font_family) && (!is.character(font_family) || length(font_family) != 1)) {
    stop('Invalid `font_family` value') 
  }
  
  ## Prep plotting data ----
  d <- do.call(
    rbind,
    lapply(
      seq_along(args),
      function(idx) data.frame(
        'Population' = as.factor(names(args[[idx]])),
        'Method' = as.factor(methods[idx]),
        'Value' = as.vector(args[[idx]])
      )
    )
  )
  d <- d[!d$Population %in% exclude, ]
  
  if (length(methods) > 1) {
    title <- 'xNPE comparison'
  } else {
    title <- 'xNPE'
  }
  
  ## Generate plot ----
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Population, y = Value, fill = Method)) +
    ggplot2::geom_bar(position = 'dodge', stat = 'identity', colour = 'black', linewidth = .1) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::ggtitle(title)
  
  ## Resolve Cartesian vs polar coords ----
  if (polar) {
    p <- p + ggplot2::coord_polar(theta = 'x') +
      ggplot2::theme(
        text = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(size = 12)
      )
  } else {
    p <- p + ggplot2::theme(
      text = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust= 0.5),
      axis.text = ggplot2::element_text(size = 12)
    )
  }
  
  ## Resolve font family ----
  if (!is.null(font_family)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }
  
  p
}

#' Compute neighbourhood composition
#'
#' Computes composition of near neighbourhood of a labelled population in terms of other populations in the data.
#' To that end, a separate neighbour search that excludes unassigned cells and cells from the reference population is conducted.
#' This means, for instance, that a neighbourhood of 'Monocytes' excluding 'Monocytes' themselves, as well as 'Debris/Doublets' is searched.
#'
#' The computed neighbourhood composition is then relevant for a specific dataset, a population within it and either the original high-dimensional data or a specific embedding of it.
#'
#' @param X numeric matrix: row-wise coordinate matrix of high-dimensional data or lower-dimensional embedding of it
#' @param annot string/factor vector: vector of labels per row of \code{X}
#' @param pop string: name of reference cell population
#' @param k integer: size of neighbourhood (default: 1000)
#' @param stepsize integer: increment in neighbourhood size per which to aggregate neighbourhood composition
#' @param exclude string vector: labelled populations to exclude (default: empty)
#'
#' @export
NeighbourhoodComposition <- function(
  X,
  annot,
  pop,
  k = 1000,
  stepsize = 10,
  exclude = c()
) {
  
  if (!is.matrix(X) || !is.numeric(X)) {
    stop('Invalid `X` value')
  }
  if ((!is.character(annot) && !is.factor(annot))) {
    stop('Invalid `annot` value')
  }
  annot <- as.factor(annot)
  if (length(annot) != nrow(X)) {
    stop('Length of `annot` must match number of rows of `X`')
  }
  if (!is.character(pop) || length(pop) != 1) {
    stop('Invalid `pop` value')
  }
  if (!pop %in% annot) {
    stop('Population "', pop, '" not found in `annot`')
  }
  if (!.IsInt(k)) {
    stop('Invalid `k` value')
  }
  if (!.IsInt(stepsize) || stepsize < 1) {
    stop('Invalid `stepsize` value')
  }
  if (length(exclude) > 0 && !is.character(exclude)) {
    stop('Invalid `exclude` value')
  }
  au <- exclude[!exclude %in% levels(annot)]
  if (length(au) > 0) {
    warning('The following `exclude` populations were not found in `annot`:\n\t', paste0(au, collapse = ',\n\t'))
  }
  
  ## Exclude self & other specified pops ----
  idcs         <- which(annot == pop)
  idcs_exclude <- c(idcs, which(annot %in% exclude))
  
  ## Resolve included pops ----
  pops <- as.vector(sort(unique(annot)))
  pops <- as.factor(pops)
  idx_self <- which(pops == pop)
  pops <- pops[!pops %in% c(exclude, pop)]
  
  ## Separate reference & query for k-NN ----
  annot_ref    <- annot[-idcs_exclude]
  X_ref   <- X[-idcs_exclude, ]
  X_query <- X[idcs, ]
  
  ## Get nearest relevant neighbours ----
  nn <-
    BiocNeighbors::queryKNN(
      X = X_ref,
      query = X_query,
      k = k,
      get.distance = FALSE,
      BNPARAM = BiocNeighbors::AnnoyParam()
    )$index
  
  ## Divide by sections ----
  step1 <- seq(from = 1, to = k, by = stepsize)
  step2 <- c(step1[-1] - 1, k)
  nn_section <- function(from, to) {
    x <- annot_ref[nn[, from:to]]
    sapply(pops, function(pop) sum(x == pop))
  }
  
  ## Get labels ----
  res <- do.call(rbind, lapply(
    seq_along(step1),
    function(idx) nn_section(from = step1[idx], to = step2[idx]))
  )
  
  ## Save metadata ----
  colnames(res) <- pops
  attributes(res)$steps <- step2
  attributes(res)$populations <- pops
  attributes(res)$self <- pop
  attributes(res)$idx_self <- idx_self
  
  class(res) <- 'NeighbourhoodComposition'
  res
}

#' Plot neighbourhood composition 
#'
#' Plots neighbourhood composition for a population of interest, as computed by function \code{NeighbourhoodComposition}.
#' A stacked area plot that is smoothed along the x-axis (neighbourhood size) is produced.
#' (The plot is not cumulative along the x-axis.)
#'
#' @param nc neighbourhood composition object, result of \code{ViScoreR::NeighbourhoodComposition} function call
#' @param stepsize integer: increment in neighbourhood size per which to aggregate neighbourhood composition
#' @param exclude string vector: labelled populations to exclude (default: empty)
#' @param k optional integer: size of neighbourhood cut-off
#' @param palette string vector: palette of colours per population to use for plotting (default palette is provided)
#' @param method optional string: name of the dimensionality reduction method or original dataset to include in title (default: none)
#' @param no_legend logical: disable displaying the legend mapping colours to populations? (default: FALSE)
#' @param font_family optional string: non-default font family to use for plotting
#'
#' @export
PlotNC <- function(
  nc,
  stepsize = 10,
  exclude = c(),
  k = NULL,
  palette = .Palette(),
  method = NULL,
  no_legend = FALSE,
  font_family = NULL
) {
  
  ## Check arguments ----
  if (class(nc) != 'NeighbourhoodComposition') {
    stop('Invalid `nc` value')
  }
  if (!is.null(k) && !.IsInt(k)) {
    stop('Invalid `k` value')
  }
  if (!.IsInt(stepsize) || stepsize < 1) {
    stop('Invalid `stepsize` value')
  }
  if (length(exclude) > 0 && !is.character(exclude)) {
    stop('Invalid `exclude` value')
  }
  pops <- levels(attributes(nc)$populations)
  au <- exclude[!exclude %in% pops]
  if (length(au) > 0) {
    warning('The following `exclude` populations are not present in `nc`:\n\t', paste0(au, collapse = ',\n\t'))
  }
  if (!is.character(palette) || !is.vector(palette)) {
    stop('Invalid `palette` value')
  }
  if (!is.null(method) && (!is.character(method) || length(method) != 1)) {
    stop('Invalid `method` value')
  }
  if (!is.logical(no_legend) || length(no_legend) != 1) {
    stop('Invalid `no_legend` value')
  }
  if (!is.null(font_family) && (!is.character(font_family) || length(font_family) != 1)) {
    stop('Invalid `font_family` value') 
  }
  
  ## Gather metadata ----
  steps <- attributes(nc)$steps
  pops  <- attributes(nc)$populations
  self  <- attributes(nc)$self
  idx_self <- attributes(nc)$idx_self
  nc    <- nc / rowSums(nc)
  
  ## Prep plotting data ----
  class(nc) <- 'matrix'
  nc    <- tidyr::pivot_longer(
    as.data.frame(nc),
    cols = tidyr::everything(),
    names_to = 'Population',
    values_to = 'Share'
  )
  npops <- length(pops)
  K <- rep(steps, each = npops)
  nc <- cbind('K' = K, nc)
  
  ## Resolve plot title ----
  title <- paste0(self, ' neighbourhood composition')
  if (!is.null(method)) {
    title <- paste0(title, ' in ', method) 
  }
  
  ## Create plot ----
  p <- ggplot2::ggplot(nc, ggplot2::aes(x = K, y = Share, fill = Population)) +
    ggplot2::geom_area() +
    ggplot2::scale_fill_manual(values = palette[-idx_self]) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 12)
    )
  attributes(p)$nc <- nc
  
  ## Resolve legend ----
  if (no_legend) {
    p <- p + ggplot2::theme(legend.position = 'none')
  }
  
  ## Resolve font family ----
  if (!is.null(font_family)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }
  
  p
}

#' Compute xNPE (Extended Neighbourhood-Proportion-Error)
#'
#' Computes xNPE (based in the Neighborhood Proportion Error) for dimensionality reduction (DR) result.
#' The xNPE gives a shape distortion score for each labelled population, computed in terms of how the proportion of alike versus not-alike cells changes in the reference population's neighbourhood as we go from the original data to the DR embedding.
#' 
#' It is advisable to pre-compute the k-nearest-neighbour graph for the original data as well as the embedding and provide them to this function.
#' (Use the \code{ViScoreR::ComputekNN} function for this.)
#'
#' @param hd numeric matrix: row-wise coordinate matrix of high-dimensional data
#' @param ld numeric matrix: row-wise coordinate matrix of low-dimensional embedding of \code{hd}
#' @param annot string/factor vector: vector of labels per row of \code{hd}, \code{ld}
#' @param knn_hd optional integer matrix: k-nearest-neighbour matrix of \code{hd}, result of \code{ViScoreR::ComputekNN} function call
#' @param knn_hd optional integer matrix: k-nearest-neighbour matrix of \code{hd}, result of \code{ViScoreR::ComputekNN} function call
#' @param k integer: size of neighbourhood if \code{knn_hd} or \code{knn_ld} are unspecified
#' @param exclude string vector: labelled populations to exclude (default: empty)
#'
#' @references Konstorum A, Jekel N, Vidal E, Laubenbacher R (2018) Comparative Analysis of Linear and Nonlinear Dimension Reduction Techniques on Mass Cytometry Data, bioRxiv 273862; doi: https://doi.org/10.1101/273862
#'
#' @export
xNPE <- function(
  hd,
  ld,
  annot,
  knn_hd = NULL,
  knn_ld = NULL,
  k = NULL,
  exclude = c()
) {
  
  ## Check arguments ----
  if (!is.matrix(hd) || !is.numeric(hd)) {
    stop('Invalid `hd` value')
  }
  if (!is.matrix(ld) || !is.numeric(ld)) {
    stop('Invalid `ld` value')
  }
  if (nrow(hd) != nrow(ld)) {
    stop('`hd` and `ld` must have same number of rows')
  }
  if ((!is.character(annot) && !is.factor(annot))) {
    stop('Invalid `annot` value')
  }
  if (length(annot) != nrow(hd)) {
    stop('Length of `annot` must match number of rows of `hd` and `ld`')
  }
  if (!is.null(knn_hd) && (!is.matrix(knn_hd) || !.IsInt(knn_hd))) {
    stop('Invalid `knn_hd` value')
  }
  if (nrow(knn_hd) != nrow(hd)) {
    stop('Row count of `knn_hd` must match row count of `hd`')
  }
  if (!is.null(knn_ld) && (!is.matrix(knn_ld) || !.IsInt(knn_ld))) {
    stop('Invalid `knn_ld` value')
  }
  if (nrow(knn_ld) != nrow(ld)) {
    stop('Row count of `knn_ld` must match row count of `ld`')
  }
  if (!is.null(k) && (!.IsInt(k) || k < 1)) {
    stop('Invalid `k` value')
  }
  if (length(exclude) > 0 && !is.character(exclude)) {
    stop('Invalid `exclude` value')
  }
  pops <- unique(annot)
  au <- exclude[!exclude %in% pops]
  if (length(au) > 0) {
    warning('The following `exclude` populations are not present in `annot`:\n\t', paste0(au, collapse = ',\n\t'))
  }
  
  ## Resolve k-NNGs ----
  if (is.null(k)) {
    k_hd <- NULL
    k_ld <- NULL
    if (!is.null(knn_hd)) {
      k_hd <- ncol(knn_hd)
      k <- k_hd
    }
    if (!is.null(knn_ld)) {
      k_ld <- ncol(knn_ld)
      if (!is.null(k)) {
        if (k_ld < k) {
          k <- k_ld
        }
      } else {
        k <- ncol(knn_ld)
      }
    }
    if (is.null(k_hd) && is.null(k_ld)) {
      stop('At least one of `knn_hd`, `knn_ld` or `k` must be specified')
    }
  }
  if (!is.null(knn_hd)) {
    knn_hd <- knn_hd[, 1:k]
  } else {
    message('Computing k-NNG for HD data')
    knn_hd <- ComputeKNN(hd, k = k)
  }
  if (!is.null(knn_ld)) {
    knn_ld <- knn_ld[, 1:k]
  } else {
    message('Computing k-NNG for LD data')
    knn_ld <- ComputeKNN(ld, k = k)
  }
  
  ## Get likeness distributions ----
  annot <- as.factor(annot)
  distr_hd <- .LikenessDistributions(annot, knn_hd)
  distr_ld <- .LikenessDistributions(annot, knn_ld)
  dhd <- distr_hd
  dld <- distr_ld
    
  distr_hd <- lapply(distr_hd, function(x) x / sum(x))
  distr_ld <- lapply(distr_ld, function(x) x / sum(x))

  distr_hd <- lapply(distr_hd, function(d) matrix(c(d, 1:(k+1)), ncol = 2))
  distr_ld <- lapply(distr_ld, function(d) matrix(c(d, 1:(k+1)), ncol = 2))
  
  ## Compute distribution dissimilarities ----
  res <- mapply(emdist::emd, distr_hd, distr_ld)
  
  ## Exclude unwanted populations ----
  res <- res[!names(res) %in% exclude]
  
  attributes(res)$distr_hd <- dhd
  attributes(res)$distr_ld <- dld
  
  class(res) <- 'xNPE'
  res
}

#' Compute k-nearest-neighbour graph
#'
#' Finds k nearest neighbours to each point in a set of points using an approximate algorithm.
#'
#' @param X numeric matrix: row-wise coordinate matrix
#' @param k integer: size of neighbourhood (default: 1000)
#'
#' @references Lun A (2023). BiocNeighbors: Nearest Neighbor Detection for Bioconductor Packages. R package version 1.18.0.
#'
#' @export
ComputeKNN <- function(
    X,
    k = 1000
) {
  
  ## Check arguments ----
  if (!is.matrix(X) || !is.numeric(X)) {
    stop('Invalid `X` value')
  }
  if (!.IsInt(k) || k < 1) {
    stop('Invalid `k` value')
  }
  
  ## Find nearest neighbours ----
  BiocNeighbors::findKNN(X, k = k, BNPARAM = BiocNeighbors::AnnoyParam())$index
}


.LikenessDistributions <- function(
  annot, knn
) {
  k <- ncol(knn)
  
  ## Get counts of like neighbours per point ----
  neighbours <- matrix(annot[knn], ncol = k)
  neighbours <- split(neighbours, row(neighbours))
  
  like <- mapply(function(x, y) sum(x == y), neighbours, annot)
  like <- as.vector(like)
  
  ## Split counts by population ----
  like_by_pop <- split(like, annot)
  
  ## Compute counts distribution for each population ----
  distr <- lapply(like_by_pop, function(x) as.vector(table(factor(x, levels = 0:k))))
  
  distr
}

.Palette <- function() {
  c(
    '#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D',
    '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#FDBF6F', '#FF7F00',
    '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', '#7FC97F', '#BEAED4', '#FDC086',
    '#FFFF99', '#386CB0', '#F0027F', '#BF5B17', '#666666', '#FBB4AE', '#B3CDE3',
    '#CCEBC5', '#DECBE4', '#FED9A6', '#FFFFCC', '#E5D8BD', '#FDDAEC', '#F2F2F2'
  )
}

.IsInt <- function(x) x == round(x)
