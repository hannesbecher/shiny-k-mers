
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
  agsl    = 6,
  agsh    = 9,
  akcovl  = 10,
  akcovh  = 100,
  avfl    = 1,
  avfh    = 10,
  athl    = -2,
  athh    = 0.6,
  adivl   = 0.1,
  adivh   = 100,
  axrangel = 45,
  axrangeh = 200,
  apallol  = 0.01,
  apalloh  = 0.99
)

.defaultSliderRanges <- list(gsMin=3,
                             gsMax=10,
                             kcovMin=5,
                             kcovMax=300,
                             vfMin=1,
                             vfMax=100,
                             thMin=-4,
                             thMax=1,
                             divMin=0.001,
                             divMax=100,
                             xrangeMin=1,
                             xrangeMax=500,
                             palloMin=0,
                             palloMax=1,
                             ymax=2)


modelClasses <-
  list("Diploid" = "d",
                     "Triploid (aaa)" = "traaa",
                     "Triploid (aab)" = "traab",
                     "Tetraploid (aaaa)" = "tau",
                     "Tetraploid (aabb)" = "tal"#,
                     #"Tetraploid (seg.)" = "tse"
)

# Suppress R CMD check notes for intentional package-level variables
utils::globalVariables(c("E028", "initialSpec"))

#' Run the Tetmer app server
#'
#' @param input Default argument for a \code{shiny} server function. Leave empty.
#' @param output Default argument for a \code{shiny} server function. Leave empty.
#' @param session Default argument for a \code{shiny} server function. Leave empty.
#' @param initialSpec An object of class \code{spectrum} used to initialize the reactive values of the Shiny application.
#'
#' @return NULL
#' @importFrom stats optim runif
#' @importFrom shiny reactiveValues renderPlot renderText observeEvent
#' @importFrom shiny fluidPage titlePanel fluidRow column plotOutput
#' @importFrom shiny verbatimTextOutput numericInput fileInput actionButton
#' @importFrom shiny icon wellPanel h4 checkboxInput radioButtons
#' @importFrom shiny conditionalPanel sliderInput shinyApp showNotification runApp stopApp
#' @importFrom shiny updateNumericInput
#' @importFrom graphics abline legend points text
tetServer <- function(input, output, session, initialSpec) {

  # Reactive values -- session-scoped state, replaces all <<- assignments
  rv <- reactiveValues(
    spec      = initialSpec,
    optimised = NULL,
    model     = NULL
  )

  # Helper: render the plot given current inputs and reactive state
  renderTetmerPlot <- function() {
    probs   <- getProbs(input)
    factors <- getFactors(input)

    if (input$fitmod == "man") {
      plotSpecApp(input, rv$spec)
      addvertlines(input)
      pointsFit(input, probs = probs, factors = factors)
      output$outText <- renderText(textOut(input, 0, rv$spec))
      return()
    }

    if (input$fitmod == "auto") {
      plotSpecApp(input, rv$spec)
      minFun       <- makeMinFun(input, probs = probs, factors = factors)
      startingVals <- getStartingVals(input)
      rv$optimised <- doOptimisation(input, rv$spec,
                                     minFun       = minFun,
                                     startingVals = startingVals)
      if(is.null(rv$optimised)){
        showNotification(
          "Optimisation failed after multiple attempts -- try adjusting the parameter ranges.",
          type = "error", duration = 10
        )
        return()
      }
      rv$model     <- input$mod
      addvertlines(input, rv$optimised)
      pointsFit(input, rv$optimised, probs = probs, factors = factors)
      pointsExtrap(input, rv$optimised, probs = probs, factors = factors)
      pointsContam(input, rv$optimised, rv$spec,
                   probs = probs, factors = factors)
      output$outText <- renderText(textOut(input, rv$optimised, rv$spec))
    }
  }

  output$plot <- renderPlot({
    renderTetmerPlot()
  })

  observeEvent(input$done, {
    # Build updated spectrum with fit stored (only if autofit was run)
    if(!is.null(rv$optimised)){
      fit     <- makeFitRecord(input, rv$optimised, rv$spec@k)
      updated <- addFit(rv$spec, fit)
    } else {
      updated <- rv$spec
    }
    stopApp(updated)
  })

  observeEvent(input$histFile, {
    fPath <- input$histFile$datapath
    fName <- input$histFile$name
    rv$spec <- prepareSpectrum(
      read.spectrum(fPath,
                    nam = strsplit(fName, "[.]")[[1]][[1]],
                    k = input$kVal)
    )
    updateNumericInput(session, "kVal", value = rv$spec@k)
  })

}


# User interface ####

makeUI <- function(spec){
  currentSliderRanges <- sliderRanges()
  fluidPage(titlePanel("Tetmer v2.3.2"),
                    "Fitting population parameters to k-mer spectra (by Hannes Becher)",
                    fluidRow(
                      column(8, plotOutput('plot')),
                      column(4,
                             verbatimTextOutput("outText"),

                      )
                    ),
            fluidRow(
              column(4,
                     numericInput("kVal", "Specify k-mer size before loading data", spec@k)
              ),
              column(4,
                     fileInput("histFile", "Choose or drop spectrum file")
              ),
              column(4,
                     fluidRow(
                       actionButton("done", "Done", icon=icon("check"),
                                    class="btn-success")
                     )
              )

            ),
                    fluidRow(
                      column(3,

                                              wellPanel(h4("1st: Select fitting mode and model"),
                                                        checkboxInput("showData", "Show data", value = TRUE),
                                                        radioButtons("fitmod", "Fitting mode",
                                                                     c("Manual" = "man",
                                                                       "Autofit" = "auto")
                                                        ),
                                                        radioButtons("mod", "Model", modelClasses

                                                        )
                                              )

                      ),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'man'",
                                              wellPanel(
                                                h4("2nd: Adjust plotting area, make all data peaks visible"),
                                                numericInput('txmax', 'Max multiplicity', .tetmerDefaults$txmax),
                                                numericInput('tymax', 'y axis max (x1000)', .tetmerDefaults$tymax)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'man'",
                                              wellPanel(
                                                h4("3rd: Param ranges"),
                                                numericInput('tkcov', 'Monoploid k-mer multiplicity', .tetmerDefaults$tkcov),
                                                numericInput('tvf', 'Variance factor (vf)', .tetmerDefaults$tvf,
                                                             min = 1, step = 0.1),
                                                numericInput('tth', 'theta', .tetmerDefaults$tth),
                                                numericInput('tyadj', 'Monoploid non-rep GS (Mbp)', .tetmerDefaults$tyadj)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'man') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              wellPanel(h4("4th: Only allopolyploids, adjust sub-genome split time"),
                                                        numericInput('tdiverg', 'T (in units of 2Ne)', .tetmerDefaults$tdiverg)
                                              )),
                             conditionalPanel(condition = "(input.fitmod == 'man') && (['tse'].includes(input.mod))",
                                              wellPanel(h4("5th: Only seg. allopolyploids, adjust p-allo"),
                                                        numericInput('pallo', 'p-allo', .tetmerDefaults$pallo)
                                              ))
                      ),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("2nd: Adjust the fitting area, make all data peaks visible"),
                                                        sliderInput("axrange", "x limits for fitting",
                                                                    min=currentSliderRanges$xrangeMin, max = currentSliderRanges$xrangeMax,
                                                                    value=c(.tetmerDefaults$axrangel, .tetmerDefaults$axrangeh)),
                                                        sliderInput("ymax", "y axis max (does not affect fit)",
                                                                    min=-2, max = currentSliderRanges$ymax,
                                                                    value=(-2 + currentSliderRanges$ymax)/2 + 1, step = (currentSliderRanges$ymax +2)/100 )
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("3rd: Param ranges"),

                                                        sliderInput('akcov', 'k-mer  multiplicity',
                                                                    min=currentSliderRanges$kcovMin, max = currentSliderRanges$kcovMax,
                                                                    value=c(.tetmerDefaults$akcovl, .tetmerDefaults$akcovh)),
                                                        sliderInput('avf', 'Variance factor (vf)',
                                                                    min=currentSliderRanges$vfMin, max = currentSliderRanges$vfMax,
                                                                    value=c(.tetmerDefaults$avfl, .tetmerDefaults$avfh), step = 0.1),
                                                        sliderInput('ath', "log10 of theta",
                                                                    min=currentSliderRanges$thMin, max = currentSliderRanges$thMax, step = 0.05,
                                                                    value=c(.tetmerDefaults$athl, .tetmerDefaults$athh)),
                                                        sliderInput('ayadj', 'Monoploid non-rep GS (Mbp)',
                                                                    min=currentSliderRanges$gsMin, max = currentSliderRanges$gsMax,
                                                                    value=c(.tetmerDefaults$agsl, .tetmerDefaults$agsh))
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'auto') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              wellPanel(h4("4th: Allopolyploids only, adjust sub-genome split time"),
                                                        sliderInput('adiv', 'T (in units of 2Ne)',
                                                                    min=currentSliderRanges$divMin, max=currentSliderRanges$divMax,
                                                                    value=c(.tetmerDefaults$adivl, .tetmerDefaults$adivh))
                                              )),
                             conditionalPanel(condition = "(input.fitmod == 'auto') && (['tse'].includes(input.mod))",
                                              wellPanel(h4("5th: Proportion of genome that is allopolyploid"),
                                                        sliderInput('apallo', 'p-allo',
                                                                    min=currentSliderRanges$palloMin, max=currentSliderRanges$palloMax,
                                                                    value=c(.tetmerDefaults$apallol, .tetmerDefaults$apalloh))
                                              ))

                      )
                    )

)
}



# function ####

#' Run the interactive Tetmer app
#'
#' @param sp A \code{spectrum} object as generated by \code{read.spectrum}.
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{result <- tetmer(E028)}
#' \dontrun{result <- tetmer(E030)}
tetmer <- function(sp=E028){
  spec   <- prepareSpectrum(sp)
  server <- function(input, output, session){
    tetServer(input, output, session, initialSpec = spec)
  }
  # runApp returns the value passed to stopApp() -- i.e. the updated spectrum
  result <- shiny::runApp(shinyApp(ui = makeUI(spec), server = server))
  invisible(result)
}

#' A stored Tetmer fit result
#'
#' Stores the result of a single Tetmer model fit, including the model type,
#' fitted parameters, parameter ranges used in optimisation, the x range
#' of the spectrum used for fitting, convergence information, and the
#' version of Tetmer used.
#'
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
  model       = "character",
  par         = "numeric",
  ranges      = "list",
  xrange      = "numeric",
  convergence = "numeric",
  value       = "numeric",
  k           = "numeric",
  version     = "character"
))

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
})

#' Create a fit record from the current app state
#'
#' Constructs a \code{tetmerFit} object from the current Shiny input
#' and optimisation result.
#'
#' @param input Input from Shiny GUI
#' @param optimised Result from \code{doOptimisation}
#' @param k Numeric k-mer length
#'
#' @return A \code{tetmerFit} object
#' @export
#' @importFrom methods new
#' @importFrom utils packageVersion
makeFitRecord <- function(input, optimised, k){
  ranges <- list(
    kcov  = input$akcov,
    vf  = input$avf,
    theta = input$ath,
    gs    = input$ayadj
  )
  if(input$mod %in% c("tal", "traab", "tse")){
    ranges$diverg <- input$adiv
  }
  if(input$mod == "tse"){
    ranges$pallo <- input$apallo
  }
  new("tetmerFit",
      model       = input$mod,
      par         = optimised$par,
      ranges      = ranges,
      xrange      = as.numeric(input$axrange),
      convergence = optimised$convergence,
      value       = optimised$value,
      k           = k,
      version     = as.character(packageVersion("Tetmer")))
}

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

#' Read in a k-mer spectrum
#'
#' Uses \code{read.table()} to read in a k-mer spectrum.
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
  sp <- read.table(f, ...)
  sp <- sp[sp[,1] <= cropAt,]
  names(sp) = c("mult", "count")
  spc <- new("spectrum",
             name=nam,
             data=sp,
             k=k,
             fits=list()
  )
  if(no0==FALSE) {
    return(prepareSpectrum(spc))
  } else {
    return(spc)
  }

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
#' @param theta Numeric vector of length 2: lower and upper bounds for
#'   log10 of theta per k-mer
#' @param gs Numeric vector of length 2: lower and upper bounds for
#'   log10 of haploid genome size in Mbp (e.g. \code{c(6, 9)})
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
#'                   theta  = c(-3, 0),
#'                   gs     = c(6, 9),
#'                   xrange = c(45, 200))
#' fit
#' spec <- addFit(E030, fit)
#' write.spectrum(spec, "E030_fitted.txt")
#' }
fitSpectrum <- function(spec,
                        model,
                        kcov,
                        vf,
                        theta,
                        gs,
                        xrange,
                        diverg      = c(0.1, 100),
                        pallo       = c(0.01, 0.99),
                        maxAttempts = 5) {

  # Build a minimal input list mimicking the Shiny input structure
  # so that existing helper functions can be reused without modification
  input <- list(
    fitmod  = "auto",
    mod     = model,
    akcov   = kcov,
    avf   = vf,
    ath     = theta,
    ayadj   = gs,
    axrange = xrange,
    adiv    = diverg,
    apallo  = pallo
  )

  probs   <- getProbs(input)
  factors <- getFactors(input)
  minFun  <- makeMinFun(input, probs, factors)
  sv      <- getStartingVals(input)

  result  <- doOptimisation(input, spec, minFun, sv,
                            maxAttempts = maxAttempts)

  if(is.null(result)){
    warning("fitSpectrum: optimisation failed after ",
            maxAttempts, " attempts -- try adjusting the parameter ranges.")
    return(NULL)
  }

  return(makeFitRecord(input, result, spec@k))
}


#' Plot a k-mer spectrum
#'
#' Uses the plot function to generate a scatter plot of a k-mer spectrum.
#'
#' It may be useful to set \code{log="xy"} or to limit the plotting range
#' using \code{xlim} or \code{ylim}.
#'
#' @param x An object of class \code{spectrum}
#' @param main (optional) A string passed to \code{plot()}
#' @param xlab (optional) A string passed to \code{plot()}
#' @param ylab (optional) A string passed to \code{plot()}
#' @param ... other keyword arguments to be passed to \code{plot()}
#'
#' @return NULL
#' @export
#' @examples
#' plot(E030, log="xy")
#' plot(E030, xlim=c(0,200), ylim=c(0,10000000))
plot.spectrum <- function(x,
                          main,
                          xlab,
                          ylab,
                          ...){
  if(missing(main)) main=x@name
  if(missing(xlab)) xlab="K-mer multiplicity (coverage)"
  if(missing(ylab)) ylab="K-mer count"
  plot(count ~ mult, data = x@data, main=main, xlab=xlab, ylab=ylab, ...)
}


#' Write a k-mer spectrum to a text file
#'
#' Writes a \code{spectrum} object to a plain text file in a human-readable
#' format. The file uses \code{#} comment lines to store metadata (name, k,
#' and any stored fits) followed by two columns: multiplicity and count.
#' The format is compatible with \code{read.spectrum} and can be inspected
#' with any text editor or scrolled through using \code{less}.
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
  writeLines(paste0("# Tetmer spectrum file -- written by Tetmer v",
                    packageVersion("Tetmer")), con)
  writeLines(paste0("# name: ", x@name), con)
  writeLines(paste0("# k: ",    x@k),    con)

  # -- Stored fits --
  nFits <- length(x@fits)
  writeLines(paste0("# fits: ", nFits), con)
  if(nFits > 0){
    for(i in seq_len(nFits)){
      fit <- x@fits[[i]]
      writeLines(paste0("# fit.", i, ".model: ",       fit@model),       con)
      writeLines(paste0("# fit.", i, ".k: ",           fit@k),           con)
      writeLines(paste0("# fit.", i, ".version: ",     fit@version),     con)
      writeLines(paste0("# fit.", i, ".convergence: ", fit@convergence), con)
      writeLines(paste0("# fit.", i, ".value: ",       fit@value),       con)
      writeLines(paste0("# fit.", i, ".xrange: ",
                        paste(fit@xrange, collapse = " ")),               con)
      for(nm in names(fit@par)){
        writeLines(paste0("# fit.", i, ".par.", nm, ": ", fit@par[nm]), con)
      }
      for(nm in names(fit@ranges)){
        vals <- fit@ranges[[nm]]
        writeLines(paste0("# fit.", i, ".range.", nm, ": ",
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

#' Add vertical lines to k-mer spectrum plot
#'
#' Using \code{input}, this function plots dashed vertical lines at multiples
#' of \code{n} corresponding to the sample's ploidy.
#'
#' @param input Input object from GUI
#' @param optimised The fit obtained from optim.
#' @keywords internal
#' @return NULL
addvertlines <- function(input, optimised = 0){

  # Determine number of peaks based on model
  nPeaks <- switch(input$mod,
    d     = 2,
    traaa = 3,
    traab = 3,
    tau   = 4,
    tal   = 4,
    tse   = 4
  )

  if(input$fitmod == "man"){
    cov   <- input$tkcov
    ypos  <- input$tymax * .7 * 1000
    abline(v = seq_len(nPeaks) * cov, lty = 2, col = "grey")
    for(i in seq_len(nPeaks)){
      text(cov * i, ypos, paste0(i, "x"), col = "grey", cex = 4)
    }
  }

  if(input$fitmod == "auto"){
    cov  <- as.numeric(optimised$par["cov"])
    ypos <- (10^input$ymax) * .7 * 1000000
    abline(v = seq_len(nPeaks) * cov, lty = 2, col = "grey")
    for(i in seq_len(nPeaks)){
      text(cov * i, ypos, paste0(i, "x"), col = "grey", cex = 4)
    }
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
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count")
    } else {
      plot(spec, xlim=c(0,input$txmax), ylim=c(0, input$tymax*1000),
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count", type = 'n')
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
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, (10^input$ymax)*1000000))
    } else {
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, (10^input$ymax)*1000000),
           type = 'n')
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

  diverg <- if(input$mod %in% c("tal","traab","tse")) {
    if(input$fitmod == "man") input$tdiverg else optimised$par["diverg"]
  } else NULL

  pallo <- if(input$mod == "tse") {
    if(input$fitmod == "man") input$pallo else optimised$par["pallo"]
  } else NULL

  if(input$fitmod == "man"){
    points(
      evalModel(probs, factors,
                xmin = 1, xmax = input$txmax,
                kcov = input$tkcov, vf = input$tvf,
                theta = input$tth, gs = input$tyadj,
                diverg = diverg, pallo = pallo),
      col = "red", type = 'l', lty = 1, lwd = 2
    )
  }

  if(input$fitmod == "auto"){
    points(input$axrange[1]:input$axrange[2],
      evalModel(probs, factors,
                xmin = input$axrange[1], xmax = input$axrange[2],
                kcov = optimised$par["cov"], vf = optimised$par["vf"],
                theta = optimised$par["theta"], gs = optimised$par["haplSize"],
                diverg = diverg, pallo = pallo),
      col = "red", type = 'l', lty = 1, lwd = 2
    )
  }
}

pointsExtrap <- function(input, optimised, probs, factors){

  diverg <- if(input$mod %in% c("tal", "traab","tse")) optimised$par["diverg"] else NULL
  pallo  <- if(input$mod =="tse") optimised$par["pallo"] else NULL

  points(1:input$axrange[1],
    evalModel(probs, factors,
              xmin = 1, xmax = input$axrange[1],
              kcov = optimised$par["cov"], vf = optimised$par["vf"],
              theta = optimised$par["theta"], gs = optimised$par["haplSize"],
              diverg = diverg, pallo = pallo),
    col = "red", type = 'l', lty = 2, lwd = 2
  )
}

pointsContam <- function(input, optimised, spect, probs, factors){

  diverg <- if(input$mod %in% c("tal","traab","tse")) optimised$par["diverg"] else NULL
  pallo  <- if(input$mod == "tse") optimised$par["pallo"] else NULL

  points(1:input$axrange[1],
    spect@data$count[1:input$axrange[1]] -
      evalModel(probs, factors,
                xmin = 1, xmax = input$axrange[1],
                kcov = optimised$par["cov"], vf = optimised$par["vf"],
                theta = optimised$par["theta"], gs = optimised$par["haplSize"],
                diverg = diverg, pallo = pallo),
    type = 'l', col = 4, lwd = 2
  )
}

#' Generate text for Tetmer window
#'
#' @param input Input from the GUI (contains plotting range, model, etc.)
#' @param optimised Fitted values from \code{optim}.
#' @param spec The current \code{spectrum} object.
#'
#' @return A string of the fitted parameters produced by \code{renderText()}
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

  # Model label lookup
  label <- switch(input$mod,
    d     = "DIPLOID",
    tau   = "AUTOTETRAPLOID",
    traaa = "AUTOTRIPLOID",
    traab = "ALLOTRIPLOID",
    tal   = "ALLOTETRAPLOID",
    tse   = "SEG. ALLOTET."
  )

  if(input$fitmod == "man"){

    kcov   <- input$tkcov
    vf   <- input$tvf
    theta  <- input$tth
    gs     <- input$tyadj
    diverg <- input$tdiverg
    pallo  <- input$pallo

    base <- paste0(
      label, " MODEL, MANUAL FIT",
      if(hasK) paste0("\n         k-mer length: ", k) else "",
      "\n  monoploid k-mer cov: ", kcov,
      "\n      theta per k-mer: ", theta,
      perNuc(theta),
      "\n     non-rep GS (Mbp): ", gs,
      "\n variance factor (vf): ", vf
    )

    if(input$mod %in% c("traab", "tal", "tse")){
      divergVal <- theta * diverg
      base <- paste0(base,
        "\n                    T: ", diverg,
        "\n     diverg per k-mer: ", round(divergVal, 4),
        divergPerNuc(divergVal)
      )
    }
    if(input$mod == "tse"){
      base <- paste0(base, "\n       prop. allotet.: ", pallo)
    }

    return(base)
  }

  if(input$fitmod == "auto"){

    # Extract named parameters -- no magic numbers
    kcov   <- optimised$par["cov"]
    vf   <- optimised$par["vf"]
    theta  <- optimised$par["theta"]
    gs     <- optimised$par["haplSize"]
    diverg <- optimised$par["diverg"]
    pallo  <- optimised$par["pallo"]

    base <- paste0(
      label, " MODEL, AUTO FITTED",
      if(hasK) paste0("\n         k-mer length: ", k) else "",
      "\n  monoploid k-mer cov: ", round(kcov, 3),
      "\n      theta per k-mer: ", round(theta, 4),
      perNuc(theta),
      "\n     non-rep GS (Mbp): ", round(gs, 1),
      "\n variance factor (vf): ", round(vf, 2)
    )

    if(input$mod %in% c("traab", "tal", "tse")){
      divergVal <- theta * diverg
      base <- paste0(base,
        "\n                    T: ", round(diverg, 2),
        "\n     diverg per k-mer: ", round(divergVal, 4),
        divergPerNuc(divergVal)
      )
    }
    if(input$mod == "tse"){
      base <- paste0(base,
        "\n       prop. allotet.: ", round(pallo, 2)
      )
    }

    # Starting ranges section
    ranges <- paste0(
      "\n\nSTARTING RANGES (MIN MAX)",
      "\n  monoploid k-mer cov: ", input$akcov[1], " ", input$akcov[2],
      "\nlog10 theta per k-mer: ", input$ath[1], " ", input$ath[2],
      "\n     non-rep GS (Mbp): ", input$ayadj[1], " ", input$ayadj[2],
      "\n variance factor (vf): ", input$avf[1], " ", input$avf[2],
      "\n              x range: ", input$axrange[1], " ", input$axrange[2]
    )

    if(input$mod %in% c("traab", "tal", "tse")){
      ranges <- paste0(ranges,
        "\n                    T: ", input$adiv[1], " ", input$adiv[2]
      )
    }
    if(input$mod == "tse"){
      ranges <- paste0(ranges,
        "\n       prop. allotet.: ", input$apallo[1], " ", input$apallo[2]
      )
    }

    return(paste0(base, ranges))
  }
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
  sp@data <- merge(allMults, sp@data, on="mult", all.x=TRUE)
  sp@data[is.na(sp@data[, 2]), 2] <- 0
  # Initialise fits slot if not present (e.g. old spectrum objects)
  if(!.hasSlot(sp, "fits")) sp@fits <- list()
  return(sp)
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

  diverg_idx <- if(input$mod %in% c("tal", "traab","tse")) 5 else NULL
  pallo_idx  <- if(input$mod == "tse") 6 else NULL

  function(x, xlimits, spec) {
    diverg <- if(!is.null(diverg_idx)) x[diverg_idx] else NULL
    pallo  <- if(!is.null(pallo_idx))  x[pallo_idx]  else NULL
    sum(
      (spec@data$count[xlimits[1]:xlimits[2]] -
         evalModel(probs, factors,
                   xmin = xlimits[1], xmax = xlimits[2],
                   kcov = x[1], vf = x[2],
                   theta = x[3], gs = x[4],
                   diverg = diverg, pallo = pallo)) ^ 2
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
  if(input$mod %in% c("d", "tau", "traaa")){
    return(c(
      cov      = (input$akcov[1]+input$akcov[2])/2,
      vf     = (input$avf[1]+input$avf[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^((input$ayadj[1]-6))+10^((input$ayadj[2]-6)))/2
    ))
  }
  if(input$mod %in% c("tal", "traab")){
    return(c(
      cov      = (input$akcov[1]+input$akcov[2])/2,
      vf     = (input$avf[1]+input$avf[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^(input$ayadj[1]-6)+10^(input$ayadj[2]-6))/2,
      diverg   = (input$adiv[1]+input$adiv[2])/2
    ))
  }
  if(input$mod == "tse"){
    return(c(
      cov      = (input$akcov[1]+input$akcov[2])/2,
      vf     = (input$avf[1]+input$avf[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^(input$ayadj[1]-6)+10^(input$ayadj[2]-6))/2,
      diverg   = (input$adiv[1]+input$adiv[2])/2,
      pallo    = (input$apallo[1]+input$apallo[2])/2
    ))
  }
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

  # Define bounds based on model
  if(input$mod %in% c("tal", "traab")){
    lower <- c(input$akcov[1], input$avf[1], 10^input$ath[1], 10^(input$ayadj[1]-6), input$adiv[1])
    upper <- c(input$akcov[2], input$avf[2], 10^input$ath[2], 10^(input$ayadj[2]-6), input$adiv[2])
    parNames <- c("cov", "vf", "theta", "haplSize", "diverg")
  } else if(input$mod == "tse"){
    lower <- c(input$akcov[1], input$avf[1], 10^input$ath[1], 10^(input$ayadj[1]-6), input$adiv[1], input$apallo[1])
    upper <- c(input$akcov[2], input$avf[2], 10^input$ath[2], 10^(input$ayadj[2]-6), input$adiv[2], input$apallo[2])
    parNames <- c("cov", "vf", "theta", "haplSize", "diverg", "pallo")
  } else {
    lower <- c(input$akcov[1], input$avf[1], 10^input$ath[1], 10^(input$ayadj[1]-6))
    upper <- c(input$akcov[2], input$avf[2], 10^input$ath[2], 10^(input$ayadj[2]-6))
    parNames <- c("cov", "vf", "theta", "haplSize")
  }

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
            xlimits = c(input$axrange[1], input$axrange[2]),
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
  if(input$mod=="tau")   return(factorAut)
  if(input$mod=="tal")   return(factorAll)
  if(input$mod=="tse")   return(factorTse)
  if(input$mod=="traaa") return(factorTraaa)
  if(input$mod=="traab") return(factorTraab)
  if(input$mod=="d")     return(factorDip)
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
  if(input$mod %in% c("tau", "tal", "tse")) return(probsTet)
  if(input$mod %in% c("traaa", "traab"))    return(probsTri)
  if(input$mod=="d")                        return(probsDip)
}

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

#' @keywords internal
makeExpectedSpectrum <- function(params, modelType, nam="", k=0){
  input <- params
  input$mod <- modelType
  probs   <- getProbs(input)
  factors <- getFactors(input)

  optimised <- params$optimised

  kcov   <- optimised$par["cov"]
  vf   <- optimised$par["vf"]
  theta  <- optimised$par["theta"]
  gs     <- optimised$par["haplSize"]
  diverg <- optimised$par["diverg"]
  pallo  <- optimised$par["pallo"]

  if(modelType %in% c("d", "tau", "traaa")){
    points(input$axrange[1]:input$axrange[2],
           evalModel(probs, factors,
                     xmin = input$axrange[1], xmax = input$axrange[2],
                     kcov = kcov, vf = vf, theta = theta, gs = gs),
           col="red", type='l', lty=1, lwd=2
    )
  }
  if(modelType %in% c("tal", "traab")){
    points(input$axrange[1]:input$axrange[2],
           evalModel(probs, factors,
                     xmin = input$axrange[1], xmax = input$axrange[2],
                     kcov = kcov, vf = vf, theta = theta, gs = gs,
                     diverg = diverg),
           col="red", type='l', lty=1, lwd=2
    )
  }
  if(modelType == "tse"){
    points(input$axrange[1]:input$axrange[2],
           evalModel(probs, factors,
                     xmin = input$axrange[1], xmax = input$axrange[2],
                     kcov = kcov, vf = vf, theta = theta, gs = gs,
                     diverg = diverg, pallo = pallo),
           col="red", type='l', lty=1, lwd=2
    )
  }

  return(new("spectrum",
             name=nam,
             data=data.frame(mult=integer(0), count=integer(0)),
             k=k,
             fits=list()))
}
