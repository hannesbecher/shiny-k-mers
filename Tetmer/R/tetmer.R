


# classes, methods, and related functionality -----------------------------------------------------


#' A stored Tetmer fit result
#'
#' Stores the result of a single Tetmer model fit, including the fit type,
#' model type, fitted parameters, parameter ranges used in optimisation, the
#' x range of the spectrum used for fitting, convergence information, and
#' the version of Tetmer used.
#'
#' @slot fitType A string indicating whether the fit was produced manually
#'   (\code{"manual"}) or by optimisation (\code{"auto"}).
#' @slot model A string indicating the model type (e.g. "d", "tal", "tau")
#' @slot par A named numeric vector of fitted parameters
#' @slot ranges A named list of parameter ranges used in optimisation
#' @slot xrange A numeric vector of length 2 giving the lower and upper
#'   multiplicity bounds of the spectrum used for fitting
#' @slot convergence A numeric indicating the convergence code from \code{optim}
#' @slot value A numeric indicating the minimised objective function value
#' @slot k A numeric indicating the k-mer length used
#' @slot version A string indicating the Tetmer version used
#'
#' @export
#' @importFrom methods setClass
setClass("tetmerFit", slots = list(
  fitType     = "character",
  model       = "character",
  par         = "numeric",
  ranges      = "list",
  xrange      = "numeric",
  convergence = "numeric",
  value       = "numeric",
  k           = "numeric",
  version     = "character"
), prototype = list(fitType = "auto"))

#' A named k-mer spectrum class
#'
#' @slot name A string indicating the name of the spectrum (or sample)
#' @slot data A \code{data.frame} with numeric columns \code{mult} and \code{count}.
#' @slot k A \code{numeric} indicating the k-mer length.
#' @slot fits A \code{list} of \code{tetmerFit} objects storing fit results.
#'
#' @export
#' @importFrom methods setClass
setClass("spectrum", slots = list(
  name = "character",
  data = "data.frame",
  k    = "numeric",
  fits = "list"
))

#' Show method for spectrum objects
#'
#' @param object A \code{spectrum} object
#' @export
#' @importFrom methods show
setMethod("show", "spectrum", function(object){
  cat("Tetmer spectrum:\n")
  cat("  Name:", object@name, "\n")
  if(object@k > 0){
    cat("  k-mer length:", object@k, "\n")
  } else {
    cat("  k-mer length: not specified\n")
  }
  cat("  Multiplicity range:", min(object@data$mult),
      "--", max(object@data$mult), "\n")
  nFits <- length(object@fits)
  if(nFits == 0){
    cat("  Fits: none\n")
  } else {
    cat("  Fits:", nFits, "stored\n")
    for(i in seq_len(nFits)){
      fit <- object@fits[[i]]
      cat("    Fit", i, "-- model:", fit@model,
          "| type:", fit@fitType,
          "| convergence:", fit@convergence,
          "| Tetmer version:", fit@version, "\n")
    }
  }
})

#' Show method for tetmerFit objects
#'
#' @param object A \code{tetmerFit} object
#' @export
setMethod("show", "tetmerFit", function(object){
  cat("Tetmer fit result:\n")
  cat("  Fit type:", object@fitType, "\n")
  cat("  Model:", object@model, "\n")
  cat("  k-mer length:", object@k, "\n")
  cat("  x range used:", object@xrange[1], "--", object@xrange[2], "\n")
  cat("  Tetmer version:", object@version, "\n")
  cat("  Convergence:", object@convergence, "\n")
  cat("  Objective value:", object@value, "\n")
  cat("  Parameters:\n")
  for(nm in names(object@par)){
    cat("   ", nm, ":", object@par[nm], "\n")
  }
  if(length(object@ranges) > 0){
    cat("  Ranges:\n")
    for(nm in names(object@ranges)){
      label <- switch(nm,
                      log10theta = "log10 theta",
                      log10Mbp = "log10 GS (Mbp)",
                      nm)
      cat("   ", label, ":", paste(object@ranges[[nm]], collapse = " "), "\n")
    }
  }
})


#' Add a fit to a spectrum object
#'
#' Appends a \code{tetmerFit} object to the \code{fits} slot of a
#' \code{spectrum} object.
#'
#' @param spec A \code{spectrum} object
#' @param fit A \code{tetmerFit} object
#'
#' @return An updated \code{spectrum} object with the fit appended
#' @export
addFit <- function(spec, fit){
  spec@fits <- c(spec@fits, list(fit))
  return(spec)
}



#' Retrieve a stored fit from a spectrum, with consistent validation
#'
#' Shared lookup logic used by \code{plot.spectrum} and
#' \code{residuals.spectrum} so that \code{fitIndex} handling stays
#' consistent between them.
#'
#' @param x A \code{spectrum} object
#' @param fitIndex Requested fit index. \code{NULL} is treated as \code{1}.
#' @param caller Name used in warning/error messages
#' @param requireFit If \code{TRUE}, raise an error rather than a warning
#'   when no valid fit can be found (used where a fit is mandatory, e.g.
#'   \code{residuals.spectrum}, since there is nothing to compute without
#'   one). If \code{FALSE} (used by \code{plot.spectrum}), \code{NA}
#'   suppresses fit lookup and returns \code{NULL} silently, and an
#'   unsatisfiable explicit request warns rather than errors.
#'
#' @return Either a \code{tetmerFit} object, or \code{NULL} if
#'   \code{requireFit} is \code{FALSE} and no fit could be found or was
#'   requested.
#' @keywords internal
getStoredFit <- function(x, fitIndex, caller, requireFit = FALSE){
  if(is.null(fitIndex)) fitIndex <- 1

  if(requireFit){
    if(length(x@fits) == 0){
      stop(caller, ": spectrum '", x@name, "' has no stored fits.", call. = FALSE)
    }
    if(length(fitIndex) != 1 || is.na(fitIndex) ||
       fitIndex < 1 || fitIndex > length(x@fits)){
      stop(caller, ": fitIndex ", fitIndex,
           " does not exist in this spectrum's fits.", call. = FALSE)
    }
    return(x@fits[[fitIndex]])
  }

  if(length(fitIndex) == 1 && is.na(fitIndex)) return(NULL)

  if(fitIndex >= 1 && fitIndex <= length(x@fits)){
    return(x@fits[[fitIndex]])
  }
  if(!(isTRUE(fitIndex == 1) && length(x@fits) == 0)){
    # Only warn for an explicit request that can't be satisfied --
    # stay silent for the implicit default of 1 on a spectrum with
    # no stored fits, so existing plot(spec) calls keep working.
    warning(caller, ": fitIndex ", fitIndex,
            " does not exist in this spectrum's fits -- plotting data only.")
  }
  NULL
}

#' Prepare spectrum for plotting and fitting
#'
#' @param spe A \code{spectrum} object
#'
#' @keywords internal
#'
#' @return A \code{spectrum} object.
prepareSpectrum <- function(spe){
  sp <- spe
  maxMult <- sp@data[nrow(sp@data), 1]
  if(maxMult < 500) maxMult <- 500
  allMults <- data.frame(mult=1:maxMult)
  sp@data <- merge(allMults, sp@data, by = "mult", all.x = TRUE)
  sp@data[is.na(sp@data[, 2]), 2] <- 0
  # Initialise fits slot if not present (e.g. old spectrum objects)
  if(!.hasSlot(sp, "fits")) sp@fits <- list()
  return(sp)
}


#' Write a k-mer spectrum to a text file
#'
#' Writes a \code{spectrum} object to a plain text file in a human-readable
#' format. The file uses \code{#TETMER}-prefixed comment lines to store
#' metadata (name, k, and any stored fits) followed by two columns:
#' multiplicity and count.
#' The format is compatible with \code{read.spectrum} and can be inspected
#' with any text editor or scrolled through using \code{less}. See the
#' spectrum and fit classes vignette for the \code{#TETMER format: 1}
#' specification.
#'
#' @param x A \code{spectrum} object
#' @param file A string giving the path to the output file, or a connection
#' @param ... Additional arguments (currently unused)
#'
#' @return \code{NULL}, invisibly
#' @export
#' @importFrom utils write.table packageVersion
#' @examples
#' \dontrun{
#' result <- tetmer(E030)
#' write.spectrum(result, "E030_fitted.txt")
#' }
write.spectrum <- function(x, file, ...){
  con <- file(file, open = "wt")
  on.exit(close(con))

  # -- Metadata header --
  writeLines("#TETMER format: 1", con)
  writeLines(paste0("#TETMER tetmer.version: ",
                    packageVersion("Tetmer")), con)
  writeLines(paste0("#TETMER name: ", x@name), con)
  writeLines(paste0("#TETMER k: ",    x@k),    con)

  # -- Stored fits --
  nFits <- length(x@fits)
  writeLines(paste0("#TETMER fits: ", nFits), con)
  if(nFits > 0){
    for(i in seq_len(nFits)){
      fit <- x@fits[[i]]
      writeLines(paste0("#TETMER fit.", i, ".fitType: ",     fit@fitType),     con)
      writeLines(paste0("#TETMER fit.", i, ".model: ",       fit@model),       con)
      writeLines(paste0("#TETMER fit.", i, ".k: ",           fit@k),           con)
      writeLines(paste0("#TETMER fit.", i, ".version: ",     fit@version),     con)
      writeLines(paste0("#TETMER fit.", i, ".convergence: ", fit@convergence), con)
      writeLines(paste0("#TETMER fit.", i, ".value: ",       fit@value),       con)
      writeLines(paste0("#TETMER fit.", i, ".xrange: ",
                        paste(fit@xrange, collapse = " ")),               con)
      for(nm in names(fit@par)){
        writeLines(paste0("#TETMER fit.", i, ".par.", nm, ": ", fit@par[nm]), con)
      }
      for(nm in names(fit@ranges)){
        vals <- fit@ranges[[nm]]
        writeLines(paste0("#TETMER fit.", i, ".range.", nm, ": ",
                          paste(vals, collapse = " ")), con)
      }
    }
  }

  # -- Spectrum data (two columns: mult count) --
  writeLines("# mult count", con)
  write.table(x@data, file = con,
              col.names = FALSE, row.names = FALSE,
              sep = " ", quote = FALSE)

  invisible(NULL)
}

#' Read in a k-mer spectrum
#'
#' Uses \code{read.table()} to read in a k-mer spectrum. For files written by
#' \code{write.spectrum()}, \code{#TETMER}-prefixed metadata lines are used to
#' restore the spectrum name, k-mer length, and any stored \code{tetmerFit}
#' objects.
#'
#' @param f A string indicating the path to the spectrum file
#' @param nam (optional) A string indicating the name of the spectrum
#' @param k (optional) A numeric indicating the k-mer length
#' @param cropAt (optional) An integer specifying the maximum multiplicity
#'   value to retain. K-mers with multiplicity above this value are
#'   discarded. Defaults to 1000.
#' @param no0 (optional) A logical value. If \code{FALSE} (default),
#'   \code{prepareSpectrum} is run on the data to insert missing
#'   multiplicity values. If \code{TRUE}, this step is skipped.
#' @param ... keyword arguments to be passed to \code{read.table}
#' @return A \code{spectrum} object
#' @export
#' @importFrom methods new .hasSlot
#' @importFrom utils read.table
#' @examples
#' \dontrun{testdf <- data.frame(mult=1:100, count = 100:1)}
#' \dontrun{tf <- tempfile()}
#' \dontrun{write.table(testdf, tf)}
#' \dontrun{sp <- read.spectrum(tf, nam = "Test spectrum")}
#' \dontrun{unlink(tf)}
#' \dontrun{plot(sp)}
read.spectrum <- function(f,
                          nam="MySpectrum",
                          k=0,
                          cropAt=1000,
                          no0=FALSE,
                          ...){
  namProvided <- !missing(nam)
  kProvided <- !missing(k)
  metadata <- parseTetmerSpectrumMetadata(readTetmerSpectrumHeader(f))
  if(!namProvided && !is.null(metadata$name)){
    nam <- metadata$name
  }
  if(!kProvided && !is.null(metadata$k)){
    k <- metadata$k
  }

  sp <- read.table(f, ...)
  sp <- sp[sp[,1] <= cropAt,]
  names(sp) = c("mult", "count")
  spc <- new("spectrum",
             name=nam,
             data=sp,
             k=k,
             fits=metadata$fits
  )
  if(no0==FALSE) {
    return(prepareSpectrum(spc))
  } else {
    return(spc)
  }

}

#' Read Tetmer metadata header lines from a spectrum file
#'
#' @param f A file path
#'
#' @return A character vector of leading comment lines
#' @noRd
readTetmerSpectrumHeader <- function(f){
  if(!is.character(f) || length(f) != 1){
    return(character())
  }

  con <- file(f, open = "rt")
  on.exit(close(con))
  header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if(length(line) == 0 || !grepl("^\\s*#", line)){
      break
    }
    header <- c(header, line)
  }
  header
}

#' Parse Tetmer metadata from spectrum comment lines
#'
#' @param header A character vector of leading comment lines
#'
#' @return A list with entries \code{name}, \code{k}, and \code{fits}
#' @noRd
parseTetmerSpectrumMetadata <- function(header){
  metadata <- list(name = NULL, k = NULL, fits = list())
  if(length(header) == 0){
    return(metadata)
  }
  if(!any(grepl("^\\s*#TETMER\\s+", header))){
    return(metadata)
  }

  fields <- parseTetmerHeaderFields(header)
  if(length(fields) == 0){
    return(metadata)
  }

  if("name" %in% names(fields)){
    metadata$name <- fields[["name"]]
  }
  if("k" %in% names(fields)){
    metadata$k <- parseTetmerNumeric(fields[["k"]])
  }
  if("fits" %in% names(fields)){
    nFits <- parseTetmerNumeric(fields[["fits"]])
    if(length(nFits) == 1 && !is.na(nFits) && nFits > 0){
      metadata$fits <- lapply(seq_len(nFits), parseTetmerFit, fields = fields)
    }
  }

  metadata
}

#' Parse colon-separated Tetmer header fields
#'
#' @param header A character vector of leading comment lines
#'
#' @return A named character vector of header values
#' @noRd
parseTetmerHeaderFields <- function(header){
  lines <- grep("^\\s*#TETMER\\s+", header, value = TRUE)
  lines <- sub("^\\s*#TETMER\\s+", "", lines)
  hasField <- grepl("^[^:]+:", lines)
  lines <- lines[hasField]
  if(length(lines) == 0){
    return(character())
  }

  keys <- trimws(sub(":.*$", "", lines))
  values <- trimws(sub("^[^:]+:\\s*", "", lines))
  stats::setNames(values, keys)
}

#' Parse one stored fit from Tetmer header fields
#'
#' @param index Integer fit index
#' @param fields A named character vector of header values
#'
#' @return A \code{tetmerFit} object
#' @noRd
parseTetmerFit <- function(index, fields){
  prefix <- paste0("fit.", index, ".")
  parPrefix <- paste0(prefix, "par.")
  rangePrefix <- paste0(prefix, "range.")
  fitType <- fields[[paste0(prefix, "fitType")]]
  if(is.null(fitType)){
    fitType <- "auto"
  }

  parKeys <- names(fields)[startsWith(names(fields), parPrefix)]
  fitPar <- vapply(fields[parKeys], parseTetmerNumeric, numeric(1))
  names(fitPar) <- sub(parPrefix, "", names(fitPar), fixed = TRUE)

  rangeKeys <- names(fields)[startsWith(names(fields), rangePrefix)]
  if(length(rangeKeys) == 0){
    ranges <- list()
  } else {
    ranges <- lapply(fields[rangeKeys], parseTetmerNumeric)
    names(ranges) <- sub(rangePrefix, "", names(ranges), fixed = TRUE)
  }

  new("tetmerFit",
      fitType     = fitType,
      model       = fields[[paste0(prefix, "model")]],
      par         = fitPar,
      ranges      = ranges,
      xrange      = parseTetmerNumeric(fields[[paste0(prefix, "xrange")]]),
      convergence = parseTetmerNumeric(fields[[paste0(prefix, "convergence")]]),
      value       = parseTetmerNumeric(fields[[paste0(prefix, "value")]]),
      k           = parseTetmerNumeric(fields[[paste0(prefix, "k")]]),
      version     = fields[[paste0(prefix, "version")]])
}

#' Parse numeric values from a Tetmer metadata value
#'
#' @param value A character scalar
#'
#' @return A numeric vector
#' @noRd
parseTetmerNumeric <- function(value){
  if(is.null(value)){
    return(numeric())
  }
  values <- strsplit(trimws(value), "\\s+")[[1]]
  if(all(values %in% c("NA", "NaN"))){
    return(rep(NA_real_, length(values)))
  }
  as.numeric(values)
}


#' Plot a k-mer spectrum
#'
#' Uses the plot function to generate a scatter plot of a k-mer spectrum,
#' optionally overlaying one of its stored model fits and/or the residual
#' between the spectrum and that fit.
#'
#' It may be useful to set \code{log="xy"} or to limit the plotting range
#' using \code{xlim} or \code{ylim}.
#'
#' @param x An object of class \code{spectrum}
#' @param main (optional) A string passed to \code{plot()}
#' @param xlab (optional) A string passed to \code{plot()}
#' @param ylab (optional) A string passed to \code{plot()}
#' @param fitIndex (optional) Which stored fit, if any, to overlay on the
#'   plot. \code{NULL} (the default) or \code{1} overlays the first stored
#'   fit, if one exists. \code{NA} suppresses fit plotting entirely (only
#'   meaningful when \code{residuals = FALSE}, since residual plotting
#'   requires a fit, and raises an error otherwise). Any other integer
#'   overlays the fit at that position in \code{x@fits}, if it exists; a
#'   warning is issued and no fit is drawn if it does not (again, unless
#'   \code{residuals} is requested, in which case this is an error since
#'   there is nothing to show without a fit). When a fit is plotted,
#'   \code{xlim} and \code{ylim} default to a range based on that fit's
#'   \code{xrange}, unless overridden via \code{...}. Fit peak positions
#'   are annotated with dashed vertical lines and \code{1x}, \code{2x},
#'   etc. labels, matching the Shiny app display. The fit x-range
#'   boundaries are shown as dashed blue vertical lines and described in
#'   the legend.
#' @param residuals Controls whether residuals (observed minus expected
#'   count, computed via \code{residuals.spectrum()}) are shown.
#'   \code{FALSE} (the default) plots only the data (and fit, if any), as
#'   before. \code{TRUE} replaces the main panel with the residual curve
#'   for the fit given by \code{fitIndex} -- a diagnostic view for fits
#'   obtained non-interactively, where there is no app window to inspect
#'   them in. \code{"overlay"} plots the data and fit as usual, with the
#'   residual additionally drawn on a secondary y-axis on the right-hand
#'   side of the plot.
#' @param colFit Colour for the fitted curve. Defaults to \code{"red"} for
#'   automatic fits and \code{"orange"} for manual fits.
#' @param colResid Colour for residual curves and the residual axis when
#'   \code{residuals = TRUE} or \code{residuals = "overlay"}. Defaults to
#'   \code{"darkgreen"}.
#' @param manualXRange Optional numeric vector of length 2 giving the
#'   multiplicity range to use when plotting a manual fit. If \code{NULL},
#'   manual fits default to \code{c(5, fit@xrange[2])}, clamped to a valid
#'   range. Automatic fits ignore this argument and use their stored
#'   \code{xrange}. This affects only plotting and residual diagnostics;
#'   it does not modify the stored fit.
#' @param ... other keyword arguments to be passed to \code{plot()};
#'   any \code{xlim} or \code{ylim} supplied here takes precedence over
#'   the defaults derived from the plotted fit
#'
#' @return NULL, invisibly
#' @export
#' @importFrom graphics grconvertY axis mtext axTicks
#' @examples
#' plot(E030, log="xy")
#' plot(E030, xlim=c(0,200), ylim=c(0,10000000))
#' fit  <- fitSpectrum(E030, model = "d", kcov = c(40, 80), vf = c(1, 10),
#'                      log10theta = c(-4, -1), log10Mbp = c(0, 3),
#'                      xrange = c(45, 200))
#' spec <- addFit(E030, fit)
#' plot(spec)                        # overlays the first (only) fit
#' plot(spec, fitIndex = NA)         # spectrum only, no fit overlay
#' plot(spec, residuals = TRUE)      # residual-only diagnostic plot
#' plot(spec, residuals = "overlay") # data + fit, with residual on secondary axis
plot.spectrum <- function(x,
                          main,
                          xlab,
                          ylab,
                          fitIndex = NULL,
                          residuals = FALSE,
                          colFit = NULL,
                          colResid = "darkgreen",
                          manualXRange = NULL,
                          ...){
  residualsMode <- if(isTRUE(residuals)){
    "only"
  } else if(identical(residuals, "overlay")){
    "overlay"
  } else if(isFALSE(residuals)){
    "none"
  } else {
    stop("plot.spectrum: `residuals` must be FALSE, TRUE, or \"overlay\".",
         call. = FALSE)
  }

  if(is.null(fitIndex)) fitIndex <- 1

  fit <- getStoredFit(x, fitIndex, caller = "plot.spectrum",
                      requireFit = residualsMode != "none")

  if(missing(main)){
    main <- if(residualsMode == "only"){
      paste0(x@name, " -- fit ", fitIndex, " residuals")
    } else {
      x@name
    }
  }
  if(missing(xlab)) xlab <- "K-mer multiplicity (coverage)"
  if(missing(ylab)){
    ylab <- if(residualsMode == "only") "Observed - expected k-mer count" else "K-mer count"
  }

  extraArgs <- list(...)
  legendCol <- function(col){
    if(is.numeric(col) && length(col) > 0 && is.finite(col[1]) && col[1] >= 1){
      grDevices::palette()[col[1]]
    } else {
      col[1]
    }
  }
  normalisePlotXRange <- function(xrange){
    xrange <- normaliseExpectedSpectrumXrange(xrange)
    multRange <- range(x@data$mult, na.rm = TRUE)
    if(xrange[1] < multRange[1] || xrange[2] > multRange[2]){
      stop("plot.spectrum: `manualXRange` must fall within the spectrum's multiplicity range.",
           call. = FALSE)
    }
    xrange
  }

  if(residualsMode == "overlay"){
    # The secondary residual axis on the right needs more room than R's
    # default right margin (2 lines) -- 5-6 character labels like "-20000"
    # plus the axis title were getting clipped at the device edge.
    oldPar <- graphics::par(mar = graphics::par("mar") + c(0, 0, 0, 3))
    on.exit(graphics::par(oldPar), add = TRUE)
  }

  fitInput <- NULL
  params   <- NULL
  fitY     <- NULL
  fitCol   <- NULL
  resid    <- NULL
  plotXRange <- NULL
  xmin <- xmax <- NULL

  if(!is.null(fit)){
    message("plot.spectrum: showing fit index ", fitIndex,
            " (model ", fit@model, ")")
    fitCol <- if(is.null(colFit)){
      if(identical(fit@fitType, "manual")) "orange" else "red"
    } else {
      colFit
    }

    currentVersion <- as.character(utils::packageVersion("Tetmer"))
    if(!identical(fit@version, currentVersion)){
      warning(
        "plot.spectrum: plotting a fit generated by Tetmer v", fit@version,
        " using current Tetmer v", currentVersion,
        "; model behaviour may differ.",
        call. = FALSE
      )
    }

    plotXRange <- fit@xrange
    if(identical(fit@fitType, "manual")){
      if(is.null(manualXRange)){
        manualLower <- min(5, fit@xrange[2])
        plotXRange <- c(manualLower, fit@xrange[2])
      } else {
        plotXRange <- normalisePlotXRange(manualXRange)
      }
    }
    xmin <- plotXRange[1]
    xmax <- plotXRange[2]

    fitInput <- list(mod = fit@model)
    params   <- getModelParameters(fitInput, optimised = list(par = fit@par))
    probs    <- getProbs(fitInput)
    factors  <- getFactors(fitInput)
    fitY     <- evalModel(probs, factors,
                          xmin = xmin, xmax = xmax,
                          kcov = params$cov, vf = params$vf,
                          theta = params$theta, gs = params$gs,
                          diverg = params$diverg, pallo = params$pallo)
    resid    <- residuals.spectrum(x, fitIndex = fitIndex, xrange = plotXRange)

    if(residualsMode != "only"){
      if(is.null(extraArgs$xlim)){
        extraArgs$xlim <- if(identical(fit@fitType, "manual")){
          c(xmin, xmax * 1.1)
        } else {
          c(0, xmax * 1.1)
        }
      }

      dataMax <- max(x@data$count[xmin:xmax], na.rm = TRUE)
      dataMin <- min(x@data$count[xmin:xmax], na.rm = TRUE)
      finiteFit <- fitY[is.finite(fitY)]
      fitMax <- if(length(finiteFit) > 0) max(finiteFit) else NA_real_
      useManualFitMax <- identical(fit@fitType, "manual") &&
        is.finite(fitMax) && fitMax > 0
      manualFitYMax <- if(useManualFitMax) 1.1 * fitMax else NA_real_

      if(residualsMode == "overlay"){
        if(!is.null(extraArgs$log) && grepl("y", extraArgs$log, fixed = TRUE)){
          stop("plot.spectrum: `log` cannot include \"y\" when `residuals = \"overlay\"`.",
               call. = FALSE)
        }
        # Scale the residual so it shares the panel with the data/fit curve,
        # anchoring residual = 0 at the *same* panel y-coordinate as
        # count = 0, rather than the previous affine rescaling (which
        # mapped range(resid, 0) onto the full panel range and let the two
        # zero points drift apart whenever residuals weren't symmetric
        # around zero). The residual gets up to ~40% of the data's own
        # vertical span to plot in.
        residMin  <- min(resid, 0, na.rm = TRUE)  # always <= 0
        residMax  <- max(resid, 0, na.rm = TRUE)  # always >= 0
        residSpan <- max(residMax, -residMin)
        countMax <- dataMax
        residScale <- if(residSpan > 0) (0.4 * countMax) / residSpan else 1

        if(is.null(extraArgs$ylim)){
          # Tight fit: panel floor is exactly low enough for the scaled
          # negative residual tail (or 0, if residuals are all >= 0); panel
          # ceiling is the data/fit curve's own max.
          extraArgs$ylim <- c(min(0, residMin * residScale), countMax)
        } else {
          # User supplied their own ylim -- keep it, but still need a scale
          # factor that keeps the residual within whichever side of zero
          # has less headroom, so it doesn't run off the panel.
          headroom   <- min(extraArgs$ylim[2], -extraArgs$ylim[1])
          residScale <- if(residSpan > 0 && is.finite(headroom) && headroom > 0){
            headroom / residSpan
          } else 1
        }
      } else if(is.null(extraArgs$ylim)){
        if(!is.null(extraArgs$log) && grepl("y", extraArgs$log, fixed = TRUE)){
          # prepareSpectrum() pads missing multiplicities with count = 0,
          # which is incompatible with a log y-axis; choose the smallest
          # positive count in the plotted range instead.
          positiveCounts <- x@data$count[xmin:xmax]
          positiveFit <- fitY
          positiveCounts <- positiveCounts[is.finite(positiveCounts) &
                                             positiveCounts > 0]
          positiveValues <- positiveCounts
          if(identical(fit@fitType, "manual")){
            positiveFit <- positiveFit[is.finite(positiveFit) & positiveFit > 0]
            positiveValues <- c(positiveCounts, positiveFit)
          }
          dataMin <- if(length(positiveValues) > 0) min(positiveValues) else 1
        } else if(useManualFitMax){
          dataMin <- 0
        }
        if(useManualFitMax){
          extraArgs$ylim <- c(dataMin, manualFitYMax)
        } else {
          extraArgs$ylim <- c(dataMin, dataMax)
        }
      }
    }
  }

  if(residualsMode == "only"){
    if(is.null(extraArgs$xlim)){
      extraArgs$xlim <- if(!is.null(fit) && identical(fit@fitType, "manual")){
        c(xmin, xmax * 1.1)
      } else {
        c(0, xmax * 1.1)
      }
    }
    extraArgs$col <- NULL
    plotArgs <- c(list(x = xmin:xmax, y = resid, type = "l",
                       main = main, xlab = xlab, ylab = ylab,
                       col = colResid),
                  extraArgs)
    do.call(plot, plotArgs)
    abline(h = 0, lty = 3, col = "grey40")
    abline(v = plotXRange, col = 4, lty = 2, lwd = 2)

    ypos <- grconvertY(0.9, from = "nfc", to = "user")
    drawPeakMarkers(cov = params$cov,
                    nPeaks = getPeakCount(fit@model),
                    ypos = ypos,
                    col = "grey",
                    cex = 1.2)

    return(invisible(NULL))
  }

  if(is.null(fit) && is.null(extraArgs$ylim)){
    dataMax <- max(x@data$count, na.rm = TRUE)
    dataMin <- min(x@data$count, na.rm = TRUE)
    if(!is.null(extraArgs$log) && grepl("y", extraArgs$log, fixed = TRUE)){
      positiveCounts <- x@data$count
      positiveCounts <- positiveCounts[is.finite(positiveCounts) & positiveCounts > 0]
      dataMin <- if(length(positiveCounts) > 0) min(positiveCounts) else 1
    }
    extraArgs$ylim <- c(dataMin, dataMax)
  }

  plotArgs <- c(list(count ~ mult, data = x@data,
                     main = main, xlab = xlab, ylab = ylab),
                extraArgs)
  do.call(plot, plotArgs)
  abline(h = 0, lty = 3, col = "grey40")

  if(!is.null(fit)){
    # Match app behavior by showing the fitted x-range boundaries.
    abline(v = plotXRange, col = 4, lty = 2, lwd = 2)

    points(
      xmin:xmax,
      fitY,
      col = fitCol, type = 'l', lty = 1, lwd = 2
    )

    # Match the app by showing expected peak multiplicities for the fitted model.
    ypos <- grconvertY(0.7, from = "nfc", to = "user")
    drawPeakMarkers(cov = params$cov,
                    nPeaks = getPeakCount(fit@model),
                    ypos = ypos,
                    col = "grey",
                    cex = 1.2)

    legendLabels <- c("Data", "Fit")
    dataCol <- if(!is.null(extraArgs$col)) extraArgs$col else "black"
    legendCols   <- c(legendCol(dataCol), legendCol(fitCol))
    legendLwd    <- c(1, 2)
    legendLty    <- c(0, 1)
    legendPch    <- c(1, NA)

    if(residualsMode == "overlay"){
      # Residual is drawn on the same panel using the zero-anchored
      # residScale computed above (alongside this panel's ylim), so
      # residual = 0 sits at exactly the same panel y-coordinate as
      # count = 0 -- both are shown by the single abline(h = 0) already
      # drawn above, rather than each axis having its own independent
      # notion of where zero sits.
      points(xmin:xmax, resid * residScale, col = colResid, type = "l", lty = 1, lwd = 2)
      # Use the SAME panel y-positions the left axis already chose (so both
      # axes have matching tick counts and aligned rows), just relabelled
      # in residual units by inverting the scale.
      leftTicks <- axTicks(2)
      usr <- graphics::par("usr")
      leftTicks <- leftTicks[leftTicks >= usr[3] & leftTicks <= usr[4]]
      axis(4, at = leftTicks, labels = signif(leftTicks / residScale, 2),
           col.axis = "black", col = colResid)
      mtext("Residual (data - fit)", side = 4, line = 3, col = colResid)

      legendLabels <- c(legendLabels, "Residual")
      legendCols   <- c(legendCols, legendCol(colResid))
      legendLwd    <- c(legendLwd, 2)
      legendLty    <- c(legendLty, 1)
      legendPch    <- c(legendPch, NA)
    }

    legend("topright",
           col = legendCols, lwd = legendLwd, lty = legendLty, pch = legendPch,
           legend = legendLabels)
  }

  invisible(NULL)
}



#' Plot the deviation between a fit and its observed spectrum
#'
#' A thin convenience wrapper around
#' \code{plot(x, fitIndex = fitIndex, residuals = TRUE, ...)}, kept under
#' its own name since it was the original entry point for this diagnostic.
#' New code may prefer calling \code{plot()} directly with
#' \code{residuals = TRUE} (residual only) or \code{residuals = "overlay"}
#' (data and fit, with the residual on a secondary axis).
#'
#' @param x A \code{spectrum} object with at least one stored fit
#' @param fitIndex Which stored fit to diagnose. Defaults to \code{1}. An
#'   error is raised if \code{x} has no stored fits, or if \code{fitIndex}
#'   does not refer to one of them.
#' @param ... Additional arguments passed on to \code{plot.spectrum}, e.g.
#'   \code{main}, \code{xlab}, \code{ylab}, or graphical parameters such
#'   as \code{xlim}/\code{ylim}.
#'
#' @return NULL, invisibly
#' @export
#' @examples
#' \dontrun{
#' fit  <- fitSpectrum(E030, model = "d", kcov = c(40, 80), vf = c(1, 10),
#'                      log10theta = c(-4, -1), log10Mbp = c(0, 3),
#'                      xrange = c(45, 200))
#' spec <- addFit(E030, fit)
#' diagPlot(spec)
#' }
diagPlot <- function(x, fitIndex = 1, ...){
  plot(x, fitIndex = fitIndex, residuals = TRUE, ...)
}



#' Residuals between a spectrum and one of its stored fits
#'
#' Computes the residual (observed minus expected k-mer count) between a
#' spectrum's observed data and one of its stored model fits, over a given
#' multiplicity range. This is the S3 method for the \code{spectrum} class
#' of the generic \code{stats::residuals}, and is used internally by
#' \code{plot.spectrum} when \code{residuals} is requested.
#'
#' @param object A \code{spectrum} object with at least one stored fit
#' @param fitIndex Which stored fit to compute residuals for. Defaults to
#'   \code{1}, the first stored fit. An error is raised if \code{object}
#'   has no stored fits, or if \code{fitIndex} does not refer to one of
#'   them.
#' @param xrange (optional) Numeric vector of length 2 giving the
#'   multiplicity range to compute residuals over. Defaults to the
#'   selected fit's own \code{xrange}.
#' @param ... Currently unused; present for compatibility with the
#'   \code{residuals} generic.
#'
#' @return A named numeric vector of residuals (observed minus expected
#'   k-mer count), named by multiplicity.
#' @export
#' @examples
#' \dontrun{
#' fit  <- fitSpectrum(E030, model = "d", kcov = c(40, 80), vf = c(1, 10),
#'                      log10theta = c(-4, -1), log10Mbp = c(0, 3),
#'                      xrange = c(45, 200))
#' spec <- addFit(E030, fit)
#' residuals(spec)
#' residuals(spec, xrange = c(45, 300))
#' }
residuals.spectrum <- function(object, fitIndex = 1, xrange = NULL, ...){
  fit <- getStoredFit(object, fitIndex, caller = "residuals.spectrum",
                      requireFit = TRUE)

  if(is.null(xrange)) xrange <- fit@xrange

  expected <- expectedSpectrum(fit, xrange = xrange)
  observed <- object@data$count[match(expected@data$mult, object@data$mult)]
  residual <- observed - expected@data$count
  names(residual) <- expected@data$mult
  residual
}



#' Build an expected spectrum from a model
#'
#' Generates a synthetic \code{spectrum} object representing the expected
#' k-mer counts under one of Tetmer's models.
#'
#' @param object Either a \code{tetmerFit} object or a model code such as
#'   \code{"d"} or \code{"tal"}.
#' @param ... Additional arguments passed to methods.
#' @param par For the \code{character} method, a named numeric vector of
#'   model parameters. Required names are \code{cov}, \code{vf},
#'   \code{theta}, and \code{haplSize}, plus \code{diverg} for
#'   allopolyploid models and \code{pallo} for \code{"tse"}.
#' @param xrange Numeric vector of length 2 giving the multiplicity range
#'   for the generated spectrum. For the \code{tetmerFit} method, defaults
#'   to \code{object@xrange}. For the \code{character} method, this must
#'   be supplied explicitly.
#' @param name (optional) Name for the returned \code{spectrum}. Defaults to
#'   the fit model or model code if omitted.
#' @param k (optional) k-mer length for the returned \code{spectrum}. For
#'   the \code{tetmerFit} method, defaults to \code{object@k}. For the
#'   \code{character} method, defaults to \code{0}.
#'
#' @return A \code{spectrum} object containing expected counts for the
#'   specified model.
#' @export
#' @importFrom methods setGeneric
#'
#' @examples
#' \dontrun{
#' fit <- fitSpectrum(E030, model = "d", kcov = c(5, 100), vf = c(1, 30),
#'                    log10theta = c(-3, 0), log10Mbp = c(0, 3),
#'                    xrange = c(45, 200),
#'                    maxAttempts = 1)
#' expectedSpectrum(fit)
#' expectedSpectrum("d",
#'                  par = c(cov = 30, vf = 1.5, theta = 0.04, haplSize = 200),
#'                  xrange = c(10, 80),
#'                  k = 21)
#' }
setGeneric("expectedSpectrum", function(object, ...) standardGeneric("expectedSpectrum"))

#' @rdname expectedSpectrum
#' @export
setMethod("expectedSpectrum", "tetmerFit",
          function(object, xrange = object@xrange, name = object@model, k = object@k){
            buildExpectedSpectrum(
              model = object@model,
              par = object@par,
              xrange = xrange,
              name = name,
              k = k
            )
          })

#' @rdname expectedSpectrum
#' @export
setMethod("expectedSpectrum", "character",
          function(object, par, xrange, name = object, k = 0){
            if(missing(par)){
              stop("`par` must be supplied when `object` is a model code.", call. = FALSE)
            }
            if(missing(xrange)){
              stop("`xrange` must be supplied when `object` is a model code.", call. = FALSE)
            }

            buildExpectedSpectrum(
              model = object,
              par = par,
              xrange = xrange,
              name = name,
              k = k
            )
          })

#' Build expected spectrum data from model parameters
#'
#' @param model Model code
#' @param par Named numeric vector of model parameters
#' @param xrange Numeric vector of length 2 giving multiplicity limits
#' @param name Spectrum name
#' @param k Numeric k-mer size
#'
#' @return A \code{spectrum} object
#' @keywords internal
buildExpectedSpectrum <- function(model, par, xrange, name, k){
  par <- normaliseExpectedSpectrumPar(par, model)
  xrange <- normaliseExpectedSpectrumXrange(xrange)
  modelInput <- list(mod = model)
  probs <- getProbs(modelInput)
  factors <- getFactors(modelInput)

  counts <- evalModel(
    probs, factors,
    xmin = xrange[1], xmax = xrange[2],
    kcov = par["cov"], vf = par["vf"],
    theta = par["theta"], gs = par["haplSize"],
    diverg = if("diverg" %in% names(par)) par["diverg"] else NULL,
    pallo = if("pallo" %in% names(par)) par["pallo"] else NULL
  )

  new("spectrum",
      name = name,
      data = data.frame(
        mult = seq.int(xrange[1], xrange[2]),
        count = as.numeric(counts)
      ),
      k = as.numeric(k),
      fits = list())
}

#' Normalise expected-spectrum parameters
#'
#' @param par Named numeric vector or coercible list of model parameters
#' @param model Model code
#'
#' @return Named numeric vector ordered for the given model
#' @keywords internal
normaliseExpectedSpectrumPar <- function(par, model){
  if(is.list(par)){
    par <- unlist(par, use.names = TRUE)
  }

  parNames <- names(par)
  par <- suppressWarnings(as.numeric(par))
  names(par) <- parNames

  if(is.null(names(par)) || any(!nzchar(names(par)))){
    stop("`par` must be a named numeric vector or named list.", call. = FALSE)
  }

  required <- getModelParameterNames(model)
  missingPar <- setdiff(required, names(par))
  if(length(missingPar) > 0){
    stop(
      paste0("Missing parameter(s) for model `", model, "`: ",
             paste(missingPar, collapse = ", ")),
      call. = FALSE
    )
  }

  extraPar <- setdiff(names(par), required)
  if(length(extraPar) > 0){
    stop(
      paste0("Unknown parameter(s) for model `", model, "`: ",
             paste(extraPar, collapse = ", ")),
      call. = FALSE
    )
  }

  if(anyNA(par[required])){
    stop("`par` must contain finite numeric values.", call. = FALSE)
  }

  orderedPar <- par[required]
  names(orderedPar) <- required
  orderedPar
}

#' Normalise expected-spectrum x-range
#'
#' @param xrange Numeric vector of length 2
#'
#' @return Integer vector of length 2
#' @keywords internal
normaliseExpectedSpectrumXrange <- function(xrange){
  xrange <- suppressWarnings(as.numeric(xrange))
  if(length(xrange) != 2 || anyNA(xrange)){
    stop("`xrange` must be a numeric vector of length 2.", call. = FALSE)
  }
  xrange <- as.integer(round(xrange))
  if(xrange[1] < 1 || xrange[2] < xrange[1]){
    stop("`xrange` must satisfy 1 <= xrange[1] <= xrange[2].", call. = FALSE)
  }
  xrange
}


# Default settings --------------------------------------------------------


# Default parameter values -- consolidated into a single internal list
# to avoid polluting the user namespace
.tetmerDefaults <- list(
  # Manual fitting defaults
  txmax   = 200,
  txmin   = 5,
  tymax   = 10000,
  tkcov   = 15,
  tvf     = 2,
  tth     = 0.04,
  tyadj   = 200,
  tdiverg = 30,
  pallo   = 0,
  # Auto fitting range defaults
  agsl    = 0,
  agsh    = 3,
  akcovl  = 10,
  akcovh  = 100,
  avfl    = 1,
  avfh    = 3,
  athl    = -2,
  athh    = 0.6,
  adivl   = 0.1,
  adivh   = 100,
  axrangel = 45,
  axrangeh = 200,
  apallol  = 0.01,
  apalloh  = 0.99
)

# these are wider than the default ranges, and are used to set the slider limits in the Shiny app
.defaultSliderRanges <- list(gsMin=-3,
                             gsMax=4,
                             kcovMin=5,
                             kcovMax=300,
                             vfMin=1,
                             vfMax=10,
                             thMin=-4,
                             thMax=1,
                             divMin=0.001,
                             divMax=100,
                             xrangeMin=1,
                             xrangeMax=500,
                             palloMin=0,
                             palloMax=1,
                             ymax=2)


# the model classes (types of ploidy) allowed, segmental allotetra is swiched off
# as it is overprameterised. May be switched on using allowSegTet()
modelClasses <-
  list("Diploid" = "d",
       "Triploid (aaa)" = "traaa",
       "Triploid (aab)" = "traab",
      "Tetraploid (aaaa)" = "tau",
      "Tetraploid (aabb)" = "tal"#,
      #"Tetraploid (seg.)" = "tse"
)


# Housekeeping ------------------------------------------------------------


# Suppress R CMD check notes for intentional package-level variables
utils::globalVariables(c("E028", "initialSpec"))


# Functions ---------------------------------------------------------------



#' Run the Tetmer app server
#'
#' @param input Default argument for a \code{shiny} server function. Leave empty.
#' @param output Default argument for a \code{shiny} server function. Leave empty.
#' @param session Default argument for a \code{shiny} server function. Leave empty.
#' @param initialSpec An object of class \code{spectrum} used to initialize the reactive values of the Shiny application.
#'
#' @return NULL
#' @importFrom stats optim runif
#' @importFrom graphics abline grconvertY legend points text
tetServer <- function(input, output, session, initialSpec) {

  # Reactive values -- session-scoped state, replaces all <<- assignments
  rv <- shiny::reactiveValues(
    spec      = initialSpec,
    optimised = NULL
  )

  computeAutofit <- function() {
    probs   <- getProbs(input)
    factors <- getFactors(input)
    minFun  <- makeMinFun(input, probs = probs, factors = factors)
    doOptimisation(input, rv$spec,
                   minFun       = minFun,
                   startingVals = getStartingVals(input))
  }

  shiny::observe({
    if (!identical(input$fitmod, "auto")) {
      rv$optimised <- NULL
      return()
    }

    # Only fit-related inputs are read here, so purely visual controls
    # can redraw the plot without rerunning optimisation.
    input$mod
    input$akcov
    input$avf
    input$ath
    input$ayadj
    input$axrange
    rv$spec

    if (input$mod %in% c("tal", "traab", "tse")) {
      input$adiv
    }
    if (input$mod == "tse") {
      input$apallo
    }

    rv$optimised <- computeAutofit()
    if (is.null(rv$optimised)) {
      shiny::showNotification(
        "Optimisation failed after multiple attempts -- try adjusting the parameter ranges.",
        type = "error", duration = 10
      )
    }
  })

  output$outText <- shiny::renderUI({
    if (input$fitmod == "man") {
      return(shiny::tags$pre(textOut(input, 0, rv$spec)))
    }

    if (is.null(rv$optimised)) {
      return(shiny::tags$pre(""))
    }

    textOutHtml(input, rv$optimised, rv$spec)
  })

  # Helper: render the plot given current inputs and reactive state
  renderTetmerPlot <- function() {
    probs   <- getProbs(input)
    factors <- getFactors(input)

    if (input$fitmod == "man") {
      plotSpecApp(input, rv$spec)
      addvertlines(input)
      pointsFit(input, probs = probs, factors = factors)
      return()
    }

    if (input$fitmod == "auto") {
      plotSpecApp(input, rv$spec)
      if (is.null(rv$optimised)) {
        return()
      }
      addvertlines(input, rv$optimised)
      pointsFit(input, rv$optimised, probs = probs, factors = factors)
      pointsExtrap(input, rv$optimised, probs = probs, factors = factors)
      pointsContam(input, rv$optimised, rv$spec,
                   probs = probs, factors = factors)
    }
  }

  output$plot <- shiny::renderPlot({
    renderTetmerPlot()
  })

  shiny::observeEvent(input$done, {
    # Build updated spectrum with fit stored. Manual mode stores the current
    # point values; auto mode stores the optimiser result when available.
    if(input$fitmod == "man" || !is.null(rv$optimised)){
      fit     <- makeFitRecord(input, rv$optimised, rv$spec@k)
      updated <- addFit(rv$spec, fit)
    } else {
      updated <- rv$spec
    }
    shiny::stopApp(updated)
  })

}


# User interface

makeUI <- function(){
  currentSliderRanges <- sliderRanges()
  shiny::fluidPage(shiny::titlePanel("Tetmer v2.3.2"),
                    "Fitting population parameters to k-mer spectra",
                    shiny::fluidRow(
                      shiny::column(8, shiny::plotOutput('plot')),
                      shiny::column(4,
                             shiny::uiOutput("outText"),
                             shiny::div(
                               style = "margin-top: 10px;",
                               shiny::actionButton("done", "Done", icon=shiny::icon("check"),
                                                   class="btn-success")
                             )
                      )
                    ),
                    shiny::fluidRow(
                      shiny::column(3,

                                              shiny::wellPanel(shiny::h4("1st: Select fitting mode and model"),
                                                        shiny::checkboxInput("showData", "Show data", value = TRUE),
                                                        shiny::radioButtons("fitmod", "Fitting mode",
                                                                     c("Manual" = "man",
                                                                       "Autofit" = "auto")
                                                        ),
                                                        shiny::radioButtons("mod", "Model", modelClasses

                                                        )
                                              )

                      ),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "input.fitmod == 'man'",
                                              shiny::wellPanel(
                                                shiny::h4("2nd: Adjust plotting area, make all data peaks visible"),
                                                shiny::numericInput('txmax', 'Max multiplicity', .tetmerDefaults$txmax),
                                                shiny::numericInput('tymax', 'y axis max (x1000)', .tetmerDefaults$tymax)
                                              ))),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "input.fitmod == 'man'",
                                              shiny::wellPanel(
                                                shiny::h4("3rd: Param ranges"),
                                                shiny::numericInput('tkcov', 'Monopl k-mer multiplicity (x1)', .tetmerDefaults$tkcov),
                                                shiny::numericInput('tvf', 'Variance factor (vf)', .tetmerDefaults$tvf,
                                                             min = 1, step = 0.1),
                                                shiny::numericInput('tth', 'theta', .tetmerDefaults$tth),
                                                shiny::numericInput('tyadj', 'Monopl non-rep GS (Mbp)', .tetmerDefaults$tyadj)
                                              ))),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "(input.fitmod == 'man') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              shiny::wellPanel(shiny::h4("4th: Only allopolyploids, adjust sub-genome split time"),
                                                        shiny::numericInput('tdiverg', 'T (in units of 2Ne)', .tetmerDefaults$tdiverg)
                                              )),
                             shiny::conditionalPanel(condition = "(input.fitmod == 'man') && (['tse'].includes(input.mod))",
                                              shiny::wellPanel(shiny::h4("5th: Only seg. allopolyploids, adjust p-allo"),
                                                        shiny::numericInput('pallo', 'p-allo', .tetmerDefaults$pallo)
                                              ))
                      ),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "input.fitmod == 'auto'",
                                              shiny::wellPanel(shiny::h4("2nd: Adjust the fitting area, make all data peaks visible"),
                                                        shiny::sliderInput("axrange", "x limits for fitting",
                                                                    min=currentSliderRanges$xrangeMin, max = currentSliderRanges$xrangeMax,
                                                                    value=c(.tetmerDefaults$axrangel, .tetmerDefaults$axrangeh)),
                                                        shiny::sliderInput("ymax", "y axis max (does not affect fit)",
                                                                    min=-2, max = currentSliderRanges$ymax,
                                                                    value=(-2 + currentSliderRanges$ymax)/2 + 1, step = (currentSliderRanges$ymax +2)/100 )
                                              ))),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "input.fitmod == 'auto'",
                                              shiny::wellPanel(shiny::h4("3rd: Param ranges"),

                                                        shiny::sliderInput('akcov', 'k-mer multiplicity (x1)',
                                                                    min=currentSliderRanges$kcovMin, max = currentSliderRanges$kcovMax,
                                                                    value=c(.tetmerDefaults$akcovl, .tetmerDefaults$akcovh)),
                                                        shiny::sliderInput('avf', 'Variance factor (vf)',
                                                                    min=currentSliderRanges$vfMin, max = currentSliderRanges$vfMax,
                                                                    value=c(.tetmerDefaults$avfl, .tetmerDefaults$avfh), step = 0.1),
                                                        shiny::sliderInput('ath', "log10 of theta",
                                                                    min=currentSliderRanges$thMin, max = currentSliderRanges$thMax, step = 0.05,
                                                                    value=c(.tetmerDefaults$athl, .tetmerDefaults$athh)),
                                                        shiny::sliderInput('ayadj', 'log10 monopl non-rep GS (Mbp)',
                                                                    min=currentSliderRanges$gsMin, max = currentSliderRanges$gsMax,
                                                                    value=c(.tetmerDefaults$agsl, .tetmerDefaults$agsh))
                                              ))),
                      shiny::column(3,
                             shiny::conditionalPanel(condition = "(input.fitmod == 'auto') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              shiny::wellPanel(shiny::h4("4th: Allopolyploids only, adjust sub-genome split time"),
                                                        shiny::sliderInput('adiv', 'T (in units of 2Ne)',
                                                                    min=currentSliderRanges$divMin, max=currentSliderRanges$divMax,
                                                                    value=c(.tetmerDefaults$adivl, .tetmerDefaults$adivh))
                                              )),
                             shiny::conditionalPanel(condition = "(input.fitmod == 'auto') && (['tse'].includes(input.mod))",
                                              shiny::wellPanel(shiny::h4("5th: Proportion of genome that is allopolyploid"),
                                                        shiny::sliderInput('apallo', 'p-allo',
                                                                    min=currentSliderRanges$palloMin, max=currentSliderRanges$palloMax,
                                                                    value=c(.tetmerDefaults$apallol, .tetmerDefaults$apalloh))
                                              ))

                      )
                    )

)
}



#' Run the interactive Tetmer app
#'
#' @param sp A \code{spectrum} object as generated by \code{read.spectrum}.
#'
#' @return A \code{spectrum} object. Pressing Done appends the current
#'   manual fit or current autofit result to the spectrum's \code{fits}
#'   slot. If the app is closed without pressing Done, the original input
#'   spectrum is returned unchanged.
#' @export
#'
#' @examples \dontrun{result <- tetmer(E028)}
#' \dontrun{result <- tetmer(E030)}
tetmer <- function(sp=E028){
  if(!requireNamespace("shiny", quietly = TRUE)){
    stop("The 'shiny' package is required to run the Tetmer app. ",
         "Please install it with: install.packages('shiny')",
         call. = FALSE)
  }
  # trim to a max multiplicity of 500 and fill in empyty counts with 0:
  spec   <- prepareSpectrum(sp)

  server <- function(input, output, session){
    tetServer(input, output, session, initialSpec = spec)
  }
  # runApp returns the value passed to shiny::stopApp(); if the app is closed
  # without pressing Done, runApp returns NULL.
  result <- shiny::runApp(shiny::shinyApp(ui = makeUI(), server = server))
  result <- resolveTetmerAppResult(result, spec)
  invisible(result)
}

#' Resolve tetmer app return value
#'
#' @param result Value returned by \code{shiny::runApp}
#' @param originalSpec Original \code{spectrum} used to launch the app
#'
#' @return A \code{spectrum} object
#' @keywords internal
resolveTetmerAppResult <- function(result, originalSpec){
  if(is.null(result)) return(originalSpec)
  result
}


#' Create a fit record from the current app state
#'
#' Constructs a \code{tetmerFit} object from the current Shiny input and,
#' for auto fits, an optimisation result.
#'
#' @param input Input from Shiny GUI
#' @param optimised Result from \code{doOptimisation}; leave \code{NULL}
#'   for manual fits
#' @param k Numeric k-mer length
#'
#' @return A \code{tetmerFit} object
#' @export
#' @importFrom methods new
#' @importFrom utils packageVersion
makeFitRecord <- function(input, optimised = NULL, k){
  if(is.null(optimised)){
    params <- getModelParameters(input)
    fitPar <- c(cov = params$cov,
                vf = params$vf,
                theta = params$theta,
                haplSize = params$gs)
    if(modelUsesDivergence(input$mod)) fitPar <- c(fitPar, diverg = params$diverg)
    if(modelUsesPallo(input$mod)) fitPar <- c(fitPar, pallo = params$pallo)

    return(new("tetmerFit",
      fitType     = "manual",
      model       = input$mod,
      par         = unlist(fitPar),
      ranges      = list(),
      xrange      = c(1, input$txmax),
      convergence = NA_real_,
      value       = NA_real_,
      k           = k,
      version     = as.character(packageVersion("Tetmer"))))
  }

  fitConfig <- asAutoFitConfig(input)
  new("tetmerFit",
      fitType     = "auto",
      model       = fitConfig$model,
      par         = optimised$par,
      ranges      = fitConfig$ranges,
      xrange      = fitConfig$xrange,
      convergence = optimised$convergence,
      value       = optimised$value,
      k           = k,
      version     = as.character(packageVersion("Tetmer")))
}



#' Fit a population genetic model to a k-mer spectrum non-interactively
#'
#' Runs Tetmer's optimisation machinery outside the Shiny app, allowing
#' scripted and batch analyses. Takes a spectrum object and parameter
#' bounds and returns a \code{tetmerFit} object directly.
#'
#' @param spec A \code{spectrum} object
#' @param model A string indicating the model type: \code{"d"} (diploid),
#'   \code{"tau"} (autotetraploid), \code{"tal"} (allotetraploid),
#'   \code{"traaa"} (autotriploid), \code{"traab"} (allotriploid),
#'   \code{"tse"} (segmental allotetraploid)
#' @param kcov Numeric vector of length 2: lower and upper bounds for
#'   monoploid k-mer coverage
#' @param vf Numeric vector of length 2: lower and upper bounds for
#'   variance factor (vf), where variance equals vf times the mean
#' @param log10theta Numeric vector of length 2: lower and upper bounds
#'   for log10 of theta per k-mer
#' @param log10Mbp Numeric vector of length 2: lower and upper bounds for
#'   log10 of haploid genome size in Mbp (e.g. \code{c(0, 3)})
#' @param xrange Numeric vector of length 2: lower and upper multiplicity
#'   bounds of the spectrum to use for fitting
#' @param diverg Numeric vector of length 2: lower and upper bounds for
#'   divergence time T in units of 2Ne (allopolyploid models only).
#'   Defaults to \code{c(0.1, 100)}.
#' @param pallo Numeric vector of length 2: lower and upper bounds for
#'   proportion allotetraploid (\code{"tse"} model only).
#'   Defaults to \code{c(0.01, 0.99)}.
#' @param maxAttempts Integer, maximum number of optimisation attempts
#'   with random restarts. Defaults to 5.
#'
#' @return A \code{tetmerFit} object, or \code{NULL} if all optimisation
#'   attempts failed
#' @export
#' @examples
#' \dontrun{
#' fit <- fitSpectrum(E030,
#'                   model  = "d",
#'                   kcov   = c(5, 100),
#'                   vf     = c(1, 100),
#'                   log10theta = c(-3, 0),
#'                   log10Mbp   = c(0, 3),
#'                   xrange = c(45, 200))
#' fit
#' spec <- addFit(E030, fit)
#' write.spectrum(spec, "E030_fitted.txt")
#' }
fitSpectrum <- function(spec,
                        model,
                        kcov,
                        vf,
                        log10theta,
                        log10Mbp,
                        xrange,
                        diverg      = c(0.1, 100),
                        pallo       = c(0.01, 0.99),
                        maxAttempts = 5) {
  fitConfig <- makeAutoFitConfig(
    model = model,
    kcov = kcov,
    vf = vf,
    log10theta = log10theta,
    log10Mbp = log10Mbp,
    xrange = xrange,
    diverg = diverg,
    pallo = pallo
  )

  probs   <- getProbs(fitConfig)
  factors <- getFactors(fitConfig)
  minFun  <- makeMinFun(fitConfig, probs, factors)
  sv      <- getStartingVals(fitConfig)

  result  <- doOptimisation(fitConfig, spec, minFun, sv,
                            maxAttempts = maxAttempts)

  if(is.null(result)){
    warning("fitSpectrum: optimisation failed after ",
            maxAttempts, " attempts -- try adjusting the parameter ranges.")
    return(NULL)
  }

  return(makeFitRecord(fitConfig, result, spec@k))
}




#' Get expected peak count for a model
#'
#' @param model Model code
#'
#' @return Integer number of expected peaks
#' @keywords internal
getPeakCount <- function(model){
  switch(model,
         d = 2,
         traaa = 3,
         traab = 3,
         tau = 4,
         tal = 4,
         tse = 4,
         NA_integer_)
}

#' Determine whether a model includes divergence time
#'
#' @param model Model code
#'
#' @return Logical scalar
#' @keywords internal
modelUsesDivergence <- function(model){
  model <- getModelCode(model)
  model %in% c("tal", "traab", "tse")
}

#' Determine whether a model includes p-allo
#'
#' @param model Model code
#'
#' @return Logical scalar
#' @keywords internal
modelUsesPallo <- function(model){
  identical(getModelCode(model), "tse")
}

#' Get the display label for a model
#'
#' @param model Model code
#'
#' @return Character scalar
#' @keywords internal
getModelLabel <- function(model){
  model <- getModelCode(model)
  switch(model,
         d = "DIPLOID",
         tau = "AUTOTETRAPLOID",
         traaa = "AUTOTRIPLOID",
         traab = "ALLOTRIPLOID",
         tal = "ALLOTETRAPLOID",
         tse = "SEG. ALLOTET.")
}

#' Get optimisation parameter names for a model
#'
#' @param model Model code
#'
#' @return Character vector
#' @keywords internal
getModelParameterNames <- function(model){
  model <- getModelCode(model)
  parNames <- c("cov", "vf", "theta", "haplSize")
  if(modelUsesDivergence(model)) parNames <- c(parNames, "diverg")
  if(modelUsesPallo(model)) parNames <- c(parNames, "pallo")
  parNames
}

#' Extract a model code from a supported config object
#'
#' @param x Model code or config-like object
#'
#' @return Character scalar
#' @keywords internal
getModelCode <- function(x){
  if(is.character(x) && length(x) == 1) return(x)
  if(is.list(x) && !is.null(x$model)) return(x$model)
  if(is.list(x) && !is.null(x$mod)) return(x$mod)
  stop("Could not determine model code.", call. = FALSE)
}

#' Create a normalized autofit config
#'
#' @param model Model code
#' @param kcov Numeric vector of length 2 for coverage bounds
#' @param vf Numeric vector of length 2 for variance factor bounds
#' @param log10theta Numeric vector of length 2 for log10 theta bounds
#' @param log10Mbp Numeric vector of length 2 for log10 haploid genome
#'   size bounds in Mbp
#' @param xrange Numeric vector of length 2 for fitting x-range
#' @param diverg Numeric vector of length 2 for divergence bounds
#' @param pallo Numeric vector of length 2 for p-allo bounds
#'
#' @return A named list representing autofit configuration
#' @keywords internal
makeAutoFitConfig <- function(model, kcov, vf, log10theta, log10Mbp, xrange,
                              diverg = c(0.1, 100), pallo = c(0.01, 0.99)){
  ranges <- list(
    kcov = kcov,
    vf = vf,
    log10theta = log10theta,
    log10Mbp = log10Mbp
  )
  if(modelUsesDivergence(model)) ranges$diverg <- diverg
  if(modelUsesPallo(model)) ranges$pallo <- pallo

  list(
    model = model,
    ranges = ranges,
    xrange = as.numeric(xrange)
  )
}

#' Normalize Shiny input or config into autofit config
#'
#' @param x Input from GUI or normalized autofit config
#'
#' @return A named list representing autofit configuration
#' @keywords internal
asAutoFitConfig <- function(x){
  if(is.list(x) && !is.null(x$model) && !is.null(x$ranges) && !is.null(x$xrange)){
    return(x)
  }

  if(is.list(x) && !is.null(x$mod)){
    return(makeAutoFitConfig(
      model = x$mod,
      kcov = x$akcov,
      vf = x$avf,
      log10theta = x$ath,
      log10Mbp = x$ayadj,
      xrange = x$axrange,
      diverg = x$adiv,
      pallo = x$apallo
    ))
  }

  stop("Could not convert object to autofit config.", call. = FALSE)
}

#' Get parameter range inputs for a model
#'
#' @param input Input from GUI or compatible list
#'
#' @return Named list of parameter ranges
#' @keywords internal
getModelRangeInputs <- function(input){
  asAutoFitConfig(input)$ranges
}

#' Extract model parameters from manual or optimised inputs
#'
#' @param input Input from GUI or compatible list
#' @param optimised Optimisation result; leave \code{NULL} for manual mode
#'
#' @return Named list of numeric parameter values
#' @keywords internal
getModelParameters <- function(input, optimised = NULL){
  if(is.null(optimised)){
    params <- list(
      cov = input$tkcov,
      vf = input$tvf,
      theta = input$tth,
      gs = input$tyadj
    )
    if(modelUsesDivergence(input$mod)) params$diverg <- input$tdiverg
    if(modelUsesPallo(input$mod)) params$pallo <- input$pallo
    return(params)
  }

  par <- optimised$par
  params <- list(
    cov = unname(par["cov"]),
    vf = unname(par["vf"]),
    theta = unname(par["theta"]),
    gs = unname(par["haplSize"])
  )
  if(modelUsesDivergence(input$mod)) params$diverg <- unname(par["diverg"])
  if(modelUsesPallo(input$mod)) params$pallo <- unname(par["pallo"])
  params
}

#' Build optimisation bounds and starting values
#'
#' @param input Input from GUI or compatible list
#'
#' @return A list containing \code{lower}, \code{upper}, and \code{start}
#' @keywords internal
buildOptimisationSpec <- function(input){
  fitConfig <- asAutoFitConfig(input)
  rangeInputs <- list(
    cov = fitConfig$ranges$kcov,
    vf = fitConfig$ranges$vf,
    theta = 10^fitConfig$ranges$log10theta,
    haplSize = 10^fitConfig$ranges$log10Mbp
  )
  if(modelUsesDivergence(fitConfig$model)) rangeInputs$diverg <- fitConfig$ranges$diverg
  if(modelUsesPallo(fitConfig$model)) rangeInputs$pallo <- fitConfig$ranges$pallo

  lower <- vapply(rangeInputs, function(x) x[1], numeric(1))
  upper <- vapply(rangeInputs, function(x) x[2], numeric(1))

  list(
    lower = lower,
    upper = upper,
    start = (lower + upper) / 2
  )
}

#' Evaluate the current model over an x-range
#'
#' @param input Input from GUI or compatible list
#' @param probs Expression for k-mer spectrum peak shapes
#' @param factors Expression for k-mer spectrum peak factors
#' @param xmin Integer, lower x limit
#' @param xmax Integer, upper x limit
#' @param optimised Optimisation result; leave \code{NULL} for manual mode
#'
#' @return Numeric vector of expected k-mer counts
#' @keywords internal
evalCurrentModel <- function(input, probs, factors, xmin, xmax, optimised = NULL){
  params <- getModelParameters(input, optimised = optimised)

  evalModel(probs, factors,
            xmin = xmin, xmax = xmax,
            kcov = params$cov, vf = params$vf,
            theta = params$theta, gs = params$gs,
            diverg = params$diverg, pallo = params$pallo)
}

#' Draw model peak markers onto the active plot
#'
#' @param cov Numeric, expected monoploid coverage
#' @param nPeaks Integer number of expected peaks
#' @param ypos Numeric y-position for the labels
#' @param col Character color used for lines and labels
#' @param cex Numeric text size for labels
#'
#' @return NULL
#' @keywords internal
drawPeakMarkers <- function(cov, nPeaks, ypos, col = "grey", cex = 4){
  if(is.null(cov) || length(cov) != 1 || !is.finite(cov) || cov <= 0) return(invisible(NULL))
  if(is.null(nPeaks) || length(nPeaks) != 1 || !is.finite(nPeaks) || nPeaks < 1) return(invisible(NULL))

  abline(v = seq_len(nPeaks) * cov, lty = 2, col = col)
  for(i in seq_len(nPeaks)){
    text(cov * i, ypos, paste0(i, "x"), col = col, cex = cex)
  }
  invisible(NULL)
}

addvertlines <- function(input, optimised = 0){
  nPeaks <- getPeakCount(input$mod)

  if(input$fitmod == "man"){
    cov   <- input$tkcov
    ypos  <- input$tymax * .7 * 1000
    drawPeakMarkers(cov = cov, nPeaks = nPeaks, ypos = ypos)
  }

  if(input$fitmod == "auto"){
    cov  <- as.numeric(optimised$par["cov"])
    ypos <- (10^input$ymax) * .7 * 1000000
    drawPeakMarkers(cov = cov, nPeaks = nPeaks, ypos = ypos)
  }
}

#' Plot annotated spectrum
#'
#' Wrapper around \code{plot.spectrum} adding some other elements.
#'
#' @param input Parameters from GUI
#' @param spec Spectrum object
#'
#' @keywords internal
#' @return NULL
plotSpecApp <- function(input, spec){

  if(input$fitmod=="man"){
    if(input$showData){
      plot(spec, xlim=c(0,input$txmax), ylim=c(0, input$tymax*1000),
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count",
           fitIndex = NA)
    } else {
      plot(spec, xlim=c(0,input$txmax), ylim=c(0, input$tymax*1000),
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count", type = 'n',
           fitIndex = NA)
    }
    legend("topright",
           col=c(1, 2),
           lwd=c(1, 2),
           lty=c(0, 1),
           pch=c(1, NA),
           legend=c("Data", "Fit")
    )
  }
  if(input$fitmod=="auto"){
    if(input$showData){
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, (10^input$ymax)*1000000),
           fitIndex = NA)
    } else {
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, (10^input$ymax)*1000000),
           type = 'n', fitIndex = NA)
    }
    abline(v=input$axrange)
    legend("topright",
           col=c(1, 2, 2, 4),
           lwd=c(1, 2, 2, 2),
           lty=c(0, 1, 2, 1),
           pch=c(1, NA, NA, NA),
           legend=c("Data", "Fit", "Extrap.", "Contam."))
  }

}



#' Plot the fit
#'
#' @param input From GUI
#' @param optimised From fit (only supplied in auto mode)
#' @param probs Expression for k-mer spectrum peak shapes
#' @param factors Expression for k-mer spectrum peak factors
#'
#' @keywords internal
#' @return NULL
pointsFit <- function(input, optimised = 0, probs, factors){
  if(input$fitmod == "man"){
    points(
      evalCurrentModel(input, probs, factors,
                       xmin = 1, xmax = input$txmax),
      col = "red", type = 'l', lty = 1, lwd = 2
    )
  }

  if(input$fitmod == "auto"){
    points(input$axrange[1]:input$axrange[2],
      evalCurrentModel(input, probs, factors,
                       xmin = input$axrange[1], xmax = input$axrange[2],
                       optimised = optimised),
      col = "red", type = 'l', lty = 1, lwd = 2
    )
  }
}

pointsExtrap <- function(input, optimised, probs, factors){
  points(1:input$axrange[1],
    evalCurrentModel(input, probs, factors,
                     xmin = 1, xmax = input$axrange[1],
                     optimised = optimised),
    col = "red", type = 'l', lty = 2, lwd = 2
  )
}

pointsContam <- function(input, optimised, spect, probs, factors){
  points(1:input$axrange[1],
    spect@data$count[1:input$axrange[1]] -
      evalCurrentModel(input, probs, factors,
                       xmin = 1, xmax = input$axrange[1],
                       optimised = optimised),
    type = 'l', col = 4, lwd = 2
  )
}

#' Generate text for Tetmer window
#'
#' @param input Input from the GUI (contains plotting range, model, etc.)
#' @param optimised Fitted values from \code{optim}.
#' @param spec The current \code{spectrum} object.
#'
#' @return A string of the fitted parameters produced by \code{shiny::renderText()}
#' to be displayed in the Tetmer window.
#' @keywords internal
textOut <- function(input, optimised, spec){

  k    <- spec@k
  hasK <- k > 0

  # Helper: per-nucleotide line, only shown when k is known
  perNuc <- function(val) {
    if(hasK) paste0("\n theta per nucleotide: ", round(val / k, 5)) else ""
  }
  divergPerNuc <- function(val) {
    if(hasK) paste0("\ndiverg per nucleotide: ", round(val / k, 4)) else ""
  }

  label <- getModelLabel(input$mod)

  if(input$fitmod == "man"){
    params <- getModelParameters(input)

    base <- paste0(
      label, " MODEL, MANUAL FIT",
      if(hasK) paste0("\n         k-mer length: ", k) else "",
      "\n  monoploid k-mer cov: ", params$cov,
      "\n      theta per k-mer: ", params$theta,
      perNuc(params$theta),
      "\n     non-rep GS (Mbp): ", params$gs,
      "\n variance factor (vf): ", params$vf
    )

    if(modelUsesDivergence(input$mod)){
      divergVal <- params$theta * params$diverg
      base <- paste0(base,
        "\n                    T: ", params$diverg,
        "\n     diverg per k-mer: ", round(divergVal, 4),
        divergPerNuc(divergVal)
      )
    }
    if(modelUsesPallo(input$mod)){
      base <- paste0(base, "\n       prop. allotet.: ", params$pallo)
    }

    return(base)
  }

  if(input$fitmod == "auto"){
    params <- getModelParameters(input, optimised = optimised)

    base <- paste0(
      label, " MODEL, AUTO FITTED",
      if(hasK) paste0("\n         k-mer length: ", k) else "",
      "\n  monoploid k-mer cov: ", round(params$cov, 3),
      "\n      theta per k-mer: ", round(params$theta, 4),
      perNuc(params$theta),
      "\n     non-rep GS (Mbp): ", round(params$gs, 1),
      "\n variance factor (vf): ", round(params$vf, 2)
    )

    if(modelUsesDivergence(input$mod)){
      divergVal <- params$theta * params$diverg
      base <- paste0(base,
        "\n                    T: ", round(params$diverg, 2),
        "\n     diverg per k-mer: ", round(divergVal, 4),
        divergPerNuc(divergVal)
      )
    }
    if(modelUsesPallo(input$mod)){
      base <- paste0(base,
        "\n       prop. allotet.: ", round(params$pallo, 2)
      )
    }

    # Starting ranges section
    ranges <- paste0(
      "\n\nSTARTING RANGES (MIN MAX)",
      "\n  monoploid k-mer cov: ", input$akcov[1], " ", input$akcov[2],
      "\nlog10 theta per k-mer: ", input$ath[1], " ", input$ath[2],
      "\nlog10 non-rep GS (Mbp): ", input$ayadj[1], " ", input$ayadj[2],
      "\n variance factor (vf): ", input$avf[1], " ", input$avf[2],
      "\n              x range: ", input$axrange[1], " ", input$axrange[2]
    )

    if(modelUsesDivergence(input$mod)){
      ranges <- paste0(ranges,
        "\n                    T: ", input$adiv[1], " ", input$adiv[2]
      )
    }
    if(modelUsesPallo(input$mod)){
      ranges <- paste0(ranges,
        "\n       prop. allotet.: ", input$apallo[1], " ", input$apallo[2]
      )
    }

    return(paste0(base, ranges))
  }
}

#' Generate styled text for the Tetmer window (HTML)
#'
#' Like \code{textOut()}, but returns an HTML \code{<pre>} block that can
#' highlight fitted parameters that land on the optimisation bounds.
#'
#' @param input Input from the GUI (contains plotting range, model, etc.)
#' @param optimised Fitted values from \code{optim}.
#' @param spec The current \code{spectrum} object.
#'
#' @return An HTML tag suitable for \code{shiny::renderUI()}.
#' @keywords internal
textOutHtml <- function(input, optimised, spec){
  k    <- spec@k
  hasK <- k > 0

  escapeHtml <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x <- gsub("'", "&#39;", x, fixed = TRUE)
    x
  }

  perNuc <- function(val) {
    if(hasK) paste0("\n theta per nucleotide: ", round(val / k, 5)) else ""
  }
  divergPerNuc <- function(val) {
    if(hasK) paste0("\ndiverg per nucleotide: ", round(val / k, 4)) else ""
  }

  bounds <- buildOptimisationSpec(input)

  boundStatus <- function(parName){
    if(is.null(optimised) || is.null(optimised$par) || is.null(names(optimised$par))) {
      return(list(on = FALSE, where = NULL))
    }
    if(!(parName %in% names(optimised$par))) return(list(on = FALSE, where = NULL))
    if(!(parName %in% names(bounds$lower)) || !(parName %in% names(bounds$upper))) {
      return(list(on = FALSE, where = NULL))
    }

    value <- unname(optimised$par[parName])
    lower <- unname(bounds$lower[parName])
    upper <- unname(bounds$upper[parName])

    if(!is.finite(value) || !is.finite(lower) || !is.finite(upper)) {
      return(list(on = FALSE, where = NULL))
    }

    span <- max(1, abs(upper - lower), abs(value), abs(lower), abs(upper))
    tol  <- 1e-8 * span

    onLower <- abs(value - lower) <= tol
    onUpper <- abs(value - upper) <= tol
    if(onLower) return(list(on = TRUE, where = "lower bound"))
    if(onUpper) return(list(on = TRUE, where = "upper bound"))
    list(on = FALSE, where = NULL)
  }

  fmtBounded <- function(parName, display){
    st <- boundStatus(parName)
    if(!isTRUE(st$on)) return(escapeHtml(display))

    paste0(
      "<span style=\"color:#b00020;font-weight:bold;\" title=\"at ", st$where, "\">",
      escapeHtml(display),
      "</span>"
    )
  }

  label  <- getModelLabel(input$mod)
  params <- getModelParameters(input, optimised = optimised)

  base <- paste0(
    escapeHtml(label), " MODEL, AUTO FITTED",
    if(hasK) paste0("\n         k-mer length: ", escapeHtml(as.character(k))) else "",
    "\n  monoploid k-mer cov: ", fmtBounded("cov", as.character(round(params$cov, 3))),
    "\n      theta per k-mer: ", fmtBounded("theta", as.character(round(params$theta, 4))),
    perNuc(params$theta),
    "\n     non-rep GS (Mbp): ", fmtBounded("haplSize", as.character(round(params$gs, 1))),
    "\n variance factor (vf): ", fmtBounded("vf", as.character(round(params$vf, 2)))
  )

  if(modelUsesDivergence(input$mod)){
    divergVal <- params$theta * params$diverg
    base <- paste0(base,
      "\n                    T: ", fmtBounded("diverg", as.character(round(params$diverg, 2))),
      "\n     diverg per k-mer: ", escapeHtml(as.character(round(divergVal, 4))),
      divergPerNuc(divergVal)
    )
  }
  if(modelUsesPallo(input$mod)){
    base <- paste0(base,
      "\n       prop. allotet.: ", fmtBounded("pallo", as.character(round(params$pallo, 2)))
    )
  }

  ranges <- paste0(
    "\n\nSTARTING RANGES (MIN MAX)",
    "\n  monoploid k-mer cov: ", escapeHtml(as.character(input$akcov[1])), " ", escapeHtml(as.character(input$akcov[2])),
    "\nlog10 theta per k-mer: ", escapeHtml(as.character(input$ath[1])), " ", escapeHtml(as.character(input$ath[2])),
    "\nlog10 non-rep GS (Mbp): ", escapeHtml(as.character(input$ayadj[1])), " ", escapeHtml(as.character(input$ayadj[2])),
    "\n variance factor (vf): ", escapeHtml(as.character(input$avf[1])), " ", escapeHtml(as.character(input$avf[2])),
    "\n              x range: ", escapeHtml(as.character(input$axrange[1])), " ", escapeHtml(as.character(input$axrange[2]))
  )

  if(modelUsesDivergence(input$mod)){
    ranges <- paste0(ranges,
      "\n                    T: ", escapeHtml(as.character(input$adiv[1])), " ", escapeHtml(as.character(input$adiv[2]))
    )
  }
  if(modelUsesPallo(input$mod)){
    ranges <- paste0(ranges,
      "\n       prop. allotet.: ", escapeHtml(as.character(input$apallo[1])), " ", escapeHtml(as.character(input$apallo[2]))
    )
  }

  note <- "\n\n(red = parameter at bound)"

  shiny::tags$pre(shiny::HTML(paste0(base, ranges, note)))
}




#' Generate a function to be minimised
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#' @param probs Expression for k-mer spectrum peak shapes
#' @param factors Expression for k-mer spectrum peak factors
#'
#' @return A function serving as input to \code{optim}
#' @keywords internal
makeMinFun <- function(input, probs, factors){
  fitConfig <- asAutoFitConfig(input)
  function(x, xlimits, spec) {
    if(is.null(names(x))){
      names(x) <- getModelParameterNames(fitConfig$model)
    }
    sum(
      (spec@data$count[xlimits[1]:xlimits[2]] -
         evalModel(probs, factors,
                   xmin = xlimits[1], xmax = xlimits[2],
                   kcov = x["cov"], vf = x["vf"],
                   theta = x["theta"], gs = x["haplSize"],
                   diverg = if("diverg" %in% names(x)) x["diverg"] else NULL,
                   pallo = if("pallo" %in% names(x)) x["pallo"] else NULL)) ^ 2
    )
  }
}

#' Starting values for optimisation
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#'
#' @return A named vector of starting values
#' @keywords internal
getStartingVals <- function(input){
  buildOptimisationSpec(input)$start
}

#' Optimise population genetic parameters
#'
#' Runs the L-BFGS-B optimisation algorithm with multiple random restarts
#' to find the best fit to the k-mer spectrum. On the first attempt, the
#' midpoint of the slider ranges is used as the starting value. On subsequent
#' attempts, starting values are randomly perturbed within the slider bounds.
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#' @param sp A \code{spectrum} object
#' @param minFun The objective function from \code{makeMinFun}
#' @param startingVals Named vector of starting values from \code{getStartingVals}
#' @param maxAttempts Integer, maximum number of optimisation attempts (default 5)
#'
#' @return An \code{optim} result list with named \code{$par}, or \code{NULL} if all attempts fail
#' @keywords internal
doOptimisation <- function(input, sp, minFun, startingVals, maxAttempts = 5){
  fitConfig <- asAutoFitConfig(input)
  optimisationSpec <- buildOptimisationSpec(input)
  lower <- optimisationSpec$lower
  upper <- optimisationSpec$upper
  parNames <- names(lower)

  bestResult <- NULL
  bestValue  <- Inf

  for(attempt in seq_len(maxAttempts)){

    # First attempt: midpoint starting values
    # Subsequent attempts: randomly perturb within bounds
    if(attempt == 1){
      sv <- startingVals
    } else {
      sv <- startingVals * runif(length(startingVals), 0.5, 1.5)
      sv <- pmax(pmin(sv, upper), lower)
    }

    result <- tryCatch({
      optim(sv, minFun,
            lower   = lower,
            upper   = upper,
            xlimits = fitConfig$xrange,
            spec    = sp,
            method  = "L-BFGS-B")
    }, error = function(e) NULL)

    # Keep best result across attempts
    if(!is.null(result) && result$value < bestValue){
      bestResult <- result
      bestValue  <- result$value
    }

    # Stop early if clean convergence
    if(!is.null(result) && result$convergence == 0) break
  }

  # Add names to parameters
  if(!is.null(bestResult)){
    names(bestResult$par) <- parNames
  }

  return(bestResult)
}

#' Factors for k-mer spectrum peaks
#'
#' @param input Input from GUI (contains model etc.)
#'
#' @return An expression of the factors by which the k-mer spectrum peaks are multiplied.
#' @keywords internal
getFactors <- function(input){
  model <- getModelCode(input)
  if(model=="tau")   return(factorAut)
  if(model=="tal")   return(factorAll)
  if(model=="tse")   return(factorTse)
  if(model=="traaa") return(factorTraaa)
  if(model=="traab") return(factorTraab)
  if(model=="d")     return(factorDip)
}

#' Evaluate the expected k-mer spectrum
#'
#' Core computation shared across pointsFit, pointsExtrap,
#' pointsContam and makeMinFun. Evaluates the model expressions
#' and returns the expected k-mer counts.
#'
#' @param probs Expression for peak shapes
#' @param factors Expression for peak factors
#' @param xmin Integer, lower x limit
#' @param xmax Integer, upper x limit
#' @param kcov Numeric, monoploid k-mer coverage
#' @param vf Numeric, variance factor (vf), where variance equals
#'   vf times the mean
#' @param theta Numeric, population-scaled mutation rate
#' @param gs Numeric, haploid genome size in millions
#' @param diverg Numeric, divergence time T (allopolyploids only)
#' @param pallo Numeric, proportion allotetraploid (tse model only)
#'
#' @return Numeric vector of expected k-mer counts
#' @keywords internal
evalModel <- function(probs, factors, xmin, xmax,
                      kcov, vf, theta, gs,
                      diverg = NULL, pallo = NULL) {

  probsEnv <- list(txmin = as.numeric(xmin),
                   txmax = as.numeric(xmax),
                   tkcov = as.numeric(kcov),
                   tvf = as.numeric(vf))

  factorsEnv <- list(tth = as.numeric(theta))
  if (!is.null(diverg)) factorsEnv$tdiverg <- as.numeric(diverg)
  if (!is.null(pallo))  factorsEnv$pal     <- as.numeric(pallo)

  colSums(eval(probs, envir = probsEnv) *
            eval(factors, envir = factorsEnv)) * as.numeric(gs) * 1000000
}

#' Raw k-mer spectrum peaks
#'
#' @param input Input from GUI (contains model etc.)
#'
#' @return An expression corresponding to the size and number of k-mer spectrum peaks
#'
#' @keywords internal
getProbs <- function(input){
  model <- getModelCode(input)
  if(model %in% c("tau", "tal", "tse")) return(probsTet)
  if(model %in% c("traaa", "traab"))    return(probsTri)
  if(model=="d")                        return(probsDip)
}


# UI things ---------------------------------------------------------------


#' Validate and normalise slider range overrides
#'
#' @param x A named list of slider range overrides
#'
#' @return A named list of numeric scalar slider range values
#' @keywords internal
normaliseSliderRanges <- function(x){
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x)))) {
    stop("Slider range overrides must be a named list.", call. = FALSE)
  }

  unknownNames <- setdiff(names(x), names(.defaultSliderRanges))
  if (length(unknownNames) > 0) {
    stop(
      paste0(
        "Unknown slider range name(s): ",
        paste(unknownNames, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  badLengths <- names(x)[lengths(x) != 1L]
  if (length(badLengths) > 0) {
    stop(
      paste0(
        "Slider range value(s) must have length 1: ",
        paste(badLengths, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coerced <- lapply(names(x), function(name) {
    numericValue <- suppressWarnings(as.numeric(x[[name]]))
    if (length(numericValue) != 1L || is.na(numericValue)) {
      stop(
        paste0("Slider range `", name, "` must be numeric."),
        call. = FALSE
      )
    }
    numericValue
  })

  names(coerced) <- names(x)
  coerced
}

#' Get slider ranges for Tetmer UI
#'
#' Returns the current slider range settings used to initialise the
#' Tetmer Shiny interface. These can be modified using
#' \code{setSliderRanges} before launching the app.
#'
#' @return A named list of slider range values
#' @export
#' @examples
#' sliderRanges()
sliderRanges <- function(){
  sliderOverrides <- getOption("Tetmer.sliderRanges", list())
  if (length(sliderOverrides) == 0) {
    return(.defaultSliderRanges)
  }

  utils::modifyList(
    .defaultSliderRanges,
    normaliseSliderRanges(sliderOverrides)
  )
}

#' Set slider ranges for Tetmer UI
#'
#' Updates the slider range settings used to initialise the Tetmer
#' Shiny interface for the current R session. Call this before
#' \code{tetmer()} to customise the parameter ranges for your data.
#'
#' @param x A named list of slider range overrides. Names must match
#'   those returned by \code{sliderRanges}. Values must be numeric.
#' @return NULL, invisibly
#' @export
#' @examples
#' \dontrun{
#' setSliderRanges(list(kcovMax = 500, ymax = 3))
#' }
setSliderRanges <- function(x){
  options(Tetmer.sliderRanges = normaliseSliderRanges(x))
  invisible(NULL)
}

#' Enable the segmental allotetraploid model
#'
#' Adds the segmental allotetraploid model (\code{"tse"}) to the list
#' of available models in the Tetmer UI. This model is hidden by default
#' as it is experimental. Call this function before launching the app
#' with \code{tetmer()}.
#'
#' @return NULL, invisibly
#' @importFrom utils assignInMyNamespace
#' @export
#' @examples
#' \dontrun{
#' allowSegTet()
#' tetmer(E028)
#' }
allowSegTet <- function(){
  utils::assignInMyNamespace("modelClasses",
    list("Diploid" = "d",
         "Triploid (aaa)" = "traaa",
         "Triploid (aab)" = "traab",
         "Tetraploid (aaaa)" = "tau",
         "Tetraploid (aabb)" = "tal",
         "Tetraploid (seg.)" = "tse"
    )
  )
}
