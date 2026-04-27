
# initial parameters for manual fitting ####
txmax <- 200
txmin <- 5
tymax <- 10000
tkcov <- 15
tbias <- -1.8
tth <- 0.04
tyadj <- 200
tdiverg <- 30
pallo <- 0

# starting ranges for auto fitting ####

agsl <- 6
agsh <- 9
akcovl <- 10
akcovh <- 100
abiasl <- -5
abiash <- -3
athl <- -2
athh <- 0.6
adivl <- 0.1
adivh <- 100
axrangel <- 45
axrangeh <- 200
apallol <- 0.01
apalloh <- 0.99

.sliderRanges <- list(gsMin=3,
                     gsMax=10,
                     kcovMin=5,
                     kcovMax=300,
                     biasMin=-6,
                     biasMax=2,
                     thMin=-4,
                     thMax=1,
                     divMin=0.001,
                     divMax=100,
                     xrangeMin=0,
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
utils::globalVariables(c(".spec", "E028"))

#' Run the Tetmer app server
#'
#' @param input Default argument for a \code{shiny} server function. Leave empty.
#' @param output Default argument for a \code{shiny} server function. Leave empty.
#'
#' @return NULL
#' @importFrom stats optim
#' @importFrom shiny reactiveValues renderPlot renderText observeEvent
#' @importFrom shiny fluidPage titlePanel fluidRow column plotOutput
#' @importFrom shiny verbatimTextOutput numericInput fileInput actionButton
#' @importFrom shiny icon wellPanel h4 checkboxInput radioButtons
#' @importFrom shiny conditionalPanel sliderInput shinyApp
#' @importFrom graphics abline legend points text
#' @importFrom utils assignInMyNamespace
tet.server <- function(input, output) {

  # Reactive values — session-scoped state, replaces all <<- assignments
  rv <- reactiveValues(
    spec      = .spec,
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

  observeEvent(input$saveFit, {
    outF <- paste0(getwd(), "/", rv$spec@name, ".fit.txt")
    print(paste("Saving fit to", outF))
    writeLines(
      ifelse(input$fitmod == "man",
             textOut(input, 0, rv$spec),
             textOut(input, rv$optimised, rv$spec)),
      con = outF
    )
  })

  observeEvent(input$saveFitAs, {
    outF <- file.choose(new = TRUE)
    if (!endsWith(outF, "txt")) outF <- paste0(outF, ".txt")
    print(paste("Saving fit to", outF))
    writeLines(
      ifelse(input$fitmod == "man",
             textOut(input, 0, rv$spec),
             textOut(input, rv$optimised, rv$spec)),
      con = outF
    )
  })

  observeEvent(input$histFile, {
    fPath <- input$histFile$datapath
    fName <- input$histFile$name
    rv$spec <- prepare.spectrum(
      read.spectrum(fPath,
                    nam = strsplit(fName, "[.]")[[1]][[1]],
                    k = input$kVal)
    )
  })

}


# User interface ####

makeUI <- function(spec){
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
                       actionButton("saveFit", "Save fit", icon=icon("save")),
                       actionButton("saveFitAs", "Save fit as...", icon=icon("save"))
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
                                                numericInput('txmax', 'Max multiplicity', txmax),
                                                numericInput('tymax', 'y axis max (x1000)', tymax)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'man'",
                                              wellPanel(
                                                h4("3rd: Param ranges"),
                                                numericInput('tkcov', 'Monoploid k-mer multiplicity', tkcov),
                                                numericInput('tbias', 'Peak width', tbias),
                                                numericInput('tth', 'theta', tth),
                                                numericInput('tyadj', 'Monoploid non-rep GS (Mbp)', tyadj)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'man') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              wellPanel(h4("4th: Only allopolyploids, adjust sub-genome split time"),
                                                        numericInput('tdiverg', 'T (in units of 2Ne)', tdiverg)
                                              )),
                             conditionalPanel(condition = "(input.fitmod == 'man') && (['tse'].includes(input.mod))",
                                              wellPanel(h4("5th: Only seg. allopolyploids, adjust p-allo"),
                                                        numericInput('pallo', 'p-allo', pallo)
                                              ))
                      ),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("2nd: Adjust the fitting area, make all data peaks visible"),
                                                        sliderInput("axrange", "x limits for fitting",
                                                                    min=.sliderRanges$xrangeMin, max = .sliderRanges$xrangeMax,
                                                                    value=c(axrangel,axrangeh)),
                                                        sliderInput("ymax", "y axis max (does not affect fit)",
                                                                    min=-2, max = .sliderRanges$ymax,
                                                                    value=1, step = 0.1)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("3rd: Param ranges"),

                                                        sliderInput('akcov', 'k-mer  multiplicity',
                                                                    min=.sliderRanges$kcovMin, max = .sliderRanges$kcovMax,
                                                                    value=c(akcovl, akcovh)),
                                                        sliderInput('abias', 'Peak width',
                                                                    min=.sliderRanges$biasMin, max = .sliderRanges$biasMax,
                                                                    value=c(abiasl, abiash), step = 0.1),
                                                        sliderInput('ath', "log10 of theta",
                                                                    min=.sliderRanges$thMin, max = .sliderRanges$thMax, step = 0.05,
                                                                    value=c(athl, athh)),
                                                        sliderInput('ayadj', 'Monoploid non-rep GS (Mbp)',
                                                                    min=.sliderRanges$gsMin, max = .sliderRanges$gsMax,
                                                                    value=c(agsl, agsh))
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'auto') && (['tal', 'traab', 'tse'].includes(input.mod))",
                                              wellPanel(h4("4th: Allopolyploids only, adjust sub-genome split time"),
                                                        sliderInput('adiv', 'T (in units of 2Ne)',
                                                                    min=.sliderRanges$divMin, max=.sliderRanges$divMax,
                                                                    value=c(adivl, adivh))
                                              )),
                             conditionalPanel(condition = "(input.fitmod == 'auto') && (['tse'].includes(input.mod))",
                                              wellPanel(h4("5th: Proportion of genome that is allopolyploid"),
                                                        sliderInput('apallo', 'p-allo',
                                                                    min=.sliderRanges$palloMin, max=.sliderRanges$palloMax,
                                                                    value=c(apallol, apalloh))
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
#' @examples \dontrun{tetmer(E028)}
#' \dontrun{tetmer(E030)}
tetmer <- function(sp=E028){
  assignInMyNamespace(".spec", prepare.spectrum(sp))
  shinyApp(ui = makeUI(.spec), server = tet.server)
}

#' A named k-mer spectrum class
#'
#' @slot name A string indicating the name of the spectrum (or sample)
#'
#' @slot data A \code{data.frame} with numeric columns \code{mult} and \code{count}.
#' @slot k A \code{numeric} indicating the k-mer length.
#'
#' @export
setClass("spectrum", slots=list(name="character", data="data.frame", k="numeric"))

# Initialise .spec at package level — set properly by tetmer() before app launch
.spec <- methods::new("spectrum",
                      name = "",
                      data = data.frame(mult = integer(0), count = integer(0)),
                      k    = 0L)

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
#'   \code{prepare.spectrum} is run on the data to insert missing
#'   multiplicity values. If \code{TRUE}, this step is skipped.
#' @param ... keyword arguments to be passed to \code{read.table}
#' @return A \code{spectrum} object
#' @export
#' @importFrom methods new
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
             k=k
  )
  if(no0==FALSE) {
    return(prepare.spectrum(spc))
  } else {
    return(spc)
  }

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


#' Add vertical lines to k-mer spectrum plot
#'
#' Using \code{input}, this function plots dashed vertical lines at multiples
#' of \code{n} corresponding to the sample's ploidy.
#'
#' @param input Input object from GUI
#' @param optimised The fit obtained from optim.
#' @keywords internal
#' @return NULL
addvertlines <- function(input, optimised=0){
  if(input$fitmod=="man"){
    if(input$mod=="d") {
      abline(v=c(1,2)*input$tkcov, lty=2,col="grey")
      text(input$tkcov[1], input$tymax*.7*1000, "1x", col="grey", cex=4)
      text(input$tkcov[1]*2, input$tymax*.7*1000, "2x", col="grey", cex=4)
    }
    if(input$mod %in% c("traaa", "traab")){
      abline(v=c(1,2,3)*input$tkcov, lty=2,col="grey")
      text(input$tkcov[1], input$tymax*.7*1000, "1x", col="grey", cex=4)
      text(input$tkcov[1]*2, input$tymax*.7*1000, "2x", col="grey", cex=4)
      text(input$tkcov[1]*3, input$tymax*.7*1000, "3x", col="grey", cex=4)
    }
    if(input$mod %in% c("tau", "tal", "tse")){
      abline(v=c(1,2,3,4)*input$tkcov, lty=2,col="grey")
      text(input$tkcov[1], input$tymax*.7*1000, "1x", col="grey", cex=4)
      text(input$tkcov[1]*2, input$tymax*.7*1000, "2x", col="grey", cex=4)
      text(input$tkcov[1]*3, input$tymax*.7*1000, "3x", col="grey", cex=4)
      text(input$tkcov[1]*4, input$tymax*.7*1000, "4x", col="grey", cex=4)
    }
  } #if manual mode
  if(input$fitmod=="auto"){
    if(input$mod=="d"){
      abline(v=c(1,2)*optimised$par[1], lty=2,col="grey")
      text(optimised$par[1], (10^input$ymax)*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, (10^input$ymax)*.7*1000000, "2x", col="grey", cex=4)
    }
    if(input$mod%in%c("traaa", "traab")){
      abline(v=c(1,2,3)*optimised$par[1], lty=2,col="grey")
      text(optimised$par[1], (10^input$ymax)*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, (10^input$ymax)*.7*1000000, "2x", col="grey", cex=4)
      text(optimised$par[1]*3, (10^input$ymax)*.7*1000000, "3x", col="grey", cex=4)
    }
    if(input$mod%in%c("tau", "tal", "tse")){
      abline(v=c(1,2,3,4)*optimised$par[1], lty=2,col="grey")
      text(optimised$par[1], (10^input$ymax)*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, (10^input$ymax)*.7*1000000, "2x", col="grey", cex=4)
      text(optimised$par[1]*3, (10^input$ymax)*.7*1000000, "3x", col="grey", cex=4)
      text(optimised$par[1]*4, (10^input$ymax)*.7*1000000, "4x", col="grey", cex=4)
    }
  } # if auto mode
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
pointsFit <- function(input, optimised=0, probs, factors){

  if(input$fitmod=="man"){
    if(input$mod %in% c("d", "tau", "traaa")){
      points(
        colSums(eval(probs, envir = list(txmin=1, txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                  eval(factors, envir=list(tth=input$tth)))*input$tyadj*1000000,
        col="red", type='l', lty=1, lwd=2
      )
    }
    if(input$mod %in% c("tal", "traab")){
      points(
        colSums(eval(probs, envir = list(txmin=1, txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                  eval(factors, envir=list(tth=input$tth, tdiverg=input$tdiverg)))*input$tyadj*1000000,
        col="red", type='l', lty=1, lwd=2
      )
    }
    if(input$mod == "tse"){
      points(
        colSums(eval(probs, envir = list(txmin=1, txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                  eval(factors, envir=list(tth=input$tth, tdiverg=input$tdiverg, pal=input$pallo)))*input$tyadj*1000000,
        col="red", type='l', lty=1, lwd=2
      )
    }
  }
  if(input$fitmod=="auto"){
    if(input$mod %in% c("d", "tau", "traaa")){
      points(input$axrange[1]:input$axrange[2],
             colSums(eval(probs, envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                              tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3])))*optimised$par[4]*1000000,
             col="red", type='l', lty=1, lwd=2
      )
    }
    if(input$mod %in% c("tal", "traab")){
      points(input$axrange[1]:input$axrange[2],
             colSums(eval(probs,
                          envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                       tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5]))
             )*optimised$par[4]*1000000,
             col="red", type='l', lty=1, lwd=2
      )
    }
    if(input$mod == "tse"){
      points(input$axrange[1]:input$axrange[2],
             colSums(eval(probs,
                          envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                       tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5], pal=optimised$par[6]))
             )*optimised$par[4]*1000000,
             col="red", type='l', lty=1, lwd=2
      )
    }
  }
}


# only run in auto mode
pointsExtrap <- function(input, optimised, probs, factors){
  if(input$mod %in% c("d", "tau", "traaa")){
    points(1:input$axrange[1],
           colSums(eval(probs, envir = list(txmin=1, txmax=input$axrange[1],
                                            tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3])))*optimised$par[4]*1000000,
           col="red", type='l', lty=2, lwd=2
    )
  }
  if(input$mod %in% c("tal", "traab")){
    points(1:input$axrange[1],
           colSums(eval(probs,
                        envir = list(txmin=1, txmax=input$axrange[1],
                                     tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5]))
           )*optimised$par[4]*1000000,
           col="red", type='l', lty=2, lwd=2
    )
  }
  if(input$mod == "tse"){
    points(1:input$axrange[1],
           colSums(eval(probs,
                        envir = list(txmin=1, txmax=input$axrange[1],
                                     tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5], pal=optimised$par[6]))
           )*optimised$par[4]*1000000,
           col="red", type='l', lty=2, lwd=2
    )
  }
}

pointsContam <- function(input, optimised, spect, probs, factors){
  if(input$mod %in% c("d", "tau", "traaa")){
    points(1:input$axrange[1],
           spect@data$count[1:input$axrange[1]] -
             colSums(eval(probs,
                          envir = list(txmin=1, txmax=input$axrange[1],
                                       tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3]))
             ) * optimised$par[4]*1000000,
           type='l', col=4, lwd=2
    )
  }
  if(input$mod %in% c("tal", "traab")){
    points(1:input$axrange[1],
           spect@data$count[1:input$axrange[1]] -
             colSums(eval(probs,
                          envir = list(txmin=1, txmax=input$axrange[1],
                                       tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5]))
             ) * optimised$par[4]*1000000,
           type='l', col=4, lwd=2
    )
  }
  if(input$mod == "tse"){
    points(1:input$axrange[1],
           spect@data$count[1:input$axrange[1]] -
             colSums(eval(probs,
                          envir = list(txmin=1, txmax=input$axrange[1],
                                       tkcov=optimised$par[1], tbias=optimised$par[2])) *
                       eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5], pal=optimised$par[6]))
             ) * optimised$par[4]*1000000,
           type='l', col=4, lwd=2
    )
  }
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
  if(spec@k > 0){
    if(input$fitmod == "man"){
      if(input$mod=="d"){
        return(paste(
          "DIPLOID MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      }
      if(input$mod=="tau"){
        return(paste(
          "AUTOTETRAPLOID MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      }
      if(input$mod=="traaa"){
        return(paste(
          "AUTOTRIPLOID MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      }
      if(input$mod=="traab"){
        return(paste(
          "ALLOTRIPLOID MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n                    T:", input$tdiverg,
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias,
          "\n     diverg per k-mer:", round(input$tth*input$tdiverg, 4),
          "\ndiverg per nucleotide:", round(input$tth*input$tdiverg/spec@k, 4)
        ))
      }
      if(input$mod=="tal"){
        return(paste(
          "ALLOTETRAPLOID MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n                    T:", input$tdiverg,
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias,
          "\n     diverg per k-mer:", round(input$tth*input$tdiverg, 4),
          "\ndiverg per nucleotide:", round(input$tth*input$tdiverg/spec@k, 4)
        ))
      }
      if(input$mod=="tse"){
        return(paste(
          "SEG. ALLOTET. MODEL, MANUAL FIT",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", input$tkcov,
          "\n      theta per k-mer:", input$tth,
          "\n theta per nucleotide:", round(input$tth/spec@k, 4),
          "\n                    T:", input$tdiverg,
          "\n     non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias,
          "\n     diverg per k-mer:", round(input$tth*input$tdiverg, 4),
          "\ndiverg per nucleotide:", round(input$tth*input$tdiverg/spec@k, 4),
          "\n       prop. allotet.:", input$pallo
        ))
      }
    } # if man
    if(input$fitmod == "auto"){
      if(input$mod=="tal"){
<<<<<<< HEAD
        return(paste(
          "ALLOTETRAPLOID MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n                    T:", round(optimised$par[5],2),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
          "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/spec@k, 4),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2]
        ))
      }
      if(input$mod=="tse"){
        return(paste(
          "SEG. ALLOTET. MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n                    T:", round(optimised$par[5],2),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
          "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/spec@k, 4),
          "\n       prop. allotet.:", round(optimised$par[6], 2),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2],
          "\n       prop. allotet.:", input$apallo[1], input$apallo[2]
        ))
      }
      if(input$mod=="traab"){
        return(paste(
          "ALLOTRIPLOID MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n                    T:", round(optimised$par[5],2),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
          "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/spec@k, 4),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2]
        ))
      }
      if(input$mod=="tau"){
        return(paste(
          "AUTOTETRAPLOID MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2]
        ))
      }
      if(input$mod=="traaa"){
        return(paste(
          "AUTOTRIPLOID MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2]
        ))
      }
      if(input$mod=="d"){
        return(paste(
          "DIPLOID MODEL, AUTO FITTED",
          "\n         k-mer length:", spec@k,
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n theta per nucleotide:", round(optimised$par[3]/spec@k,5),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1], input$axrange[2]
        ))
      }
=======
        return(
          paste(
            "ALLOTETRAPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n  monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n                    T:", round(optimised$par[5],2),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
            "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/.spec@k, 4),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n                    T:", input$adiv[1], input$adiv[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))

      } # if allotet
      if(input$mod=="tse"){
        return(
          paste(
            "ALLOTETRAPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n  monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n                    T:", round(optimised$par[5],2),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
            "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/.spec@k, 4),
            "\n       prop. allotet.:", round(optimised$par[6], 2),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n                    T:", input$adiv[1], input$adiv[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2],
            "\n       prop. allotet.:", input$apallo[1], input$apallo[2]
          ))

      } # if allotet
      if(input$mod=="traab"){
        return(
          paste(
            "ALLOTRIPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n  monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n                    T:", round(optimised$par[5],2),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n diverg per k-mer:", round(optimised$par[3]*optimised$par[5], 4),
            "\ndiverg per nucleotide:", round(optimised$par[3]*optimised$par[5]/.spec@k, 4),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n                    T:", input$adiv[1], input$adiv[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))

      } # if allotrip
      if(input$mod=="tau"){
        return(
          paste(
            "AUTOTETRAPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n  monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))

      } # if tau
      if(input$mod=="traaa"){
        return(
          paste(
            "AUTOTRIPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))

      } # if traab
      if(input$mod=="d"){
        return(
          paste(
            "DIPLOID MODEL, AUTO FITTED",
            "\n         k-mer length:", .spec@k,
            "\n  monoploid k-mer cov:", round(optimised$par[1],1),
            "\n      theta per k-mer:", round(optimised$par[3],4),
            "\n theta per nucleotide:", round(optimised$par[3] / .spec@k,5),
            "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],3),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
            "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))

      } #if d
>>>>>>> cf2f9a70a95d4e4486e77e9b04a41d10ced47caa
    } # if auto fit
  } # if k > 0
  # k == 0: no per-nucleotide estimates
  if(input$fitmod == "man"){
    if(input$mod=="d"){
      return(paste(
        "DIPLOID MODEL, MANUAL FIT",
        "\n  monoploid k-mer cov:", input$tkcov,
        "\n      theta per k-mer:", input$tth,
        "\n     non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias
      ))
    }
    if(input$mod=="tau"){
      return(paste(
        "AUTOTETRAPLOID MODEL, MANUAL FIT",
        "\n  monoploid k-mer cov:", input$tkcov,
        "\n      theta per k-mer:", input$tth,
        "\n     non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias
      ))
    }
    if(input$mod=="traaa"){
      return(paste(
        "AUTOTRIPLOID MODEL, MANUAL FIT",
        "\n  monoploid k-mer cov:", input$tkcov,
        "\n      theta per k-mer:", input$tth,
        "\n     non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias
      ))
    }
    if(input$mod=="traab"){
      return(paste(
        "ALLOTRIPLOID MODEL, MANUAL FIT",
        "\n  monoploid k-mer cov:", input$tkcov,
        "\n      theta per k-mer:", input$tth,
        "\n                    T:", input$tdiverg,
        "\n     non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias,
        "\n     diverg per k-mer:", round(input$tth*input$tdiverg, 4)
      ))
    }
    if(input$mod=="tal"){
      return(paste(
        "ALLOTETRAPLOID MODEL, MANUAL FIT",
        "\n  monoploid k-mer cov:", input$tkcov,
        "\n      theta per k-mer:", input$tth,
        "\n                    T:", input$tdiverg,
        "\n     non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias,
        "\n     diverg per k-mer:", round(input$tth*input$tdiverg, 4)
      ))
    }
  } # if man
  if(input$fitmod == "auto"){
    if(input$mod=="tal"){
<<<<<<< HEAD
      return(paste(
        "ALLOTETRAPLOID MODEL, AUTO FITTED",
        "\n  monoploid k-mer cov:", round(optimised$par[1],1),
        "\n      theta per k-mer:", round(optimised$par[3],4),
        "\n                    T:", round(optimised$par[5],2),
        "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
        "\n    bias (peak width):", round(optimised$par[2],1),
        "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
        "\n\nSTARTING RANGES (MIN MAX)",
        "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
        "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
        "\n                    T:", input$adiv[1], input$adiv[2],
        "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
        "\n    bias (peak width):", input$abias[1], input$abias[2],
        "\n              x range:", input$axrange[1], input$axrange[2]
      ))
    }
    if(input$mod=="traab"){
      return(paste(
        "ALLOTRIPLOID MODEL, AUTO FITTED",
        "\n  monoploid k-mer cov:", round(optimised$par[1],1),
        "\n      theta per k-mer:", round(optimised$par[3],4),
        "\n                    T:", round(optimised$par[5],2),
        "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
        "\n    bias (peak width):", round(optimised$par[2],1),
        "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
        "\n\nSTARTING RANGES (MIN MAX)",
        "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
        "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
        "\n                    T:", input$adiv[1], input$adiv[2],
        "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
        "\n    bias (peak width):", input$abias[1], input$abias[2],
        "\n              x range:", input$axrange[1], input$axrange[2]
      ))
    }
    if(input$mod=="tau"){
      return(paste(
        "AUTOTETRAPLOID MODEL, AUTO FITTED",
        "\n  monoploid k-mer cov:", round(optimised$par[1],1),
        "\n      theta per k-mer:", round(optimised$par[3],4),
        "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
        "\n    bias (peak width):", round(optimised$par[2],1),
        "\n\nSTARTING RANGES (MIN MAX)",
        "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
        "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
        "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
        "\n    bias (peak width):", input$abias[1], input$abias[2],
        "\n              x range:", input$axrange[1], input$axrange[2]
      ))
    }
    if(input$mod=="traaa"){
      return(paste(
        "AUTOTRIPLOID MODEL, AUTO FITTED",
        "\n  monoploid k-mer cov:", round(optimised$par[1],1),
        "\n      theta per k-mer:", round(optimised$par[3],4),
        "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
        "\n    bias (peak width):", round(optimised$par[2],1),
        "\n\nSTARTING RANGES (MIN MAX)",
        "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
        "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
        "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
        "\n    bias (peak width):", input$abias[1], input$abias[2],
        "\n              x range:", input$axrange[1], input$axrange[2]
      ))
    }
    if(input$mod=="d"){
      return(paste(
        "DIPLOID MODEL, AUTO FITTED",
        "\n  monoploid k-mer cov:", round(optimised$par[1],1),
        "\n      theta per k-mer:", round(optimised$par[3],4),
        "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
        "\n    bias (peak width):", round(optimised$par[2],1),
        "\n\nSTARTING RANGES (MIN MAX)",
        "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
        "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
        "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
        "\n    bias (peak width):", input$abias[1], input$abias[2],
        "\n              x range:", input$axrange[1], input$axrange[2]
      ))
    }
=======
      return(
        paste(
          "ALLOTETRAPLOID MODEL, AUTO FITTED",
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n                    T:", round(optimised$par[5],2),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],3),
          "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))

    } # if allotet
    if(input$mod=="traab"){
      return(
        paste(
          "ALLOTRIPLOID MODEL, AUTO FITTED",
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n                    T:", round(optimised$par[5],2),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],3),
          "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))

    } # if allotrip
    if(input$mod=="tau"){
      return(
        paste(
          "AUTOTETRAPLOID MODEL, AUTO FITTED",
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))

    } # if tau
    if(input$mod=="traaa"){
      return(
        paste(
          "AUTOTRIPLOID MODEL, AUTO FITTED",
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))

    } # if traab
    if(input$mod=="d"){
      return(
        paste(
          "DIPLOID MODEL, AUTO FITTED",
          "\n  monoploid k-mer cov:", round(optimised$par[1],1),
          "\n      theta per k-mer:", round(optimised$par[3],4),
          "\n     non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n  monoploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 theta per k-mer:", input$ath[1], input$ath[2],
          "\n     non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))

    } #if d
>>>>>>> cf2f9a70a95d4e4486e77e9b04a41d10ced47caa
  } # if auto fit
}



#' Prepare spectrum for plotting and fitting
#'
#' @param spe A \code{spectrum} object
#'
#' @keywords internal
#'
#' @return A \code{spectrum} object.
prepare.spectrum <- function(spe){
  sp <- spe
  maxMult <- sp@data[nrow(sp@data), 1]
  if(maxMult < 500) maxMult <- 500
  allMults <- data.frame(mult=1:maxMult)
  sp@data <- merge(allMults, sp@data, on="mult", all.x=TRUE)
  sp@data[is.na(sp@data[, 2]), 2] <- 0
  return(sp)
}


#' Generate a function to be minimised
#'
#' This is then minimised by \code{optim}.
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#' @param probs Expression for k-mer spectrum peak shapes
#' @param factors Expression for k-mer spectrum peak factors
#'
#' @return A function serving as input to \code{optim}
#'
#' @keywords internal
makeMinFun <- function(input, probs, factors){

  if(input$mod=="d"){
    return(function(x, xlimits, spec) {
      sum(
        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3])))*x[4]*1000000) ^2
      )
    })
  }
  if(input$mod %in% c("tal", "traab")){
    return(function(x, xlimits, spec) {
      sum(
        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3], tdiverg=x[5])))*x[4]*1000000) ^2
      )
    })
  }
  if(input$mod == "tse"){
    return(function(x, xlimits, spec) {
      sum(
        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3], tdiverg=x[5], pal=x[6])))*x[4]*1000000) ^2
      )
    })
  }
  if(input$mod %in% c("tau", "traaa")){
    return(function(x, xlimits, spec) {
      sum(
        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3])))*x[4]*1000000) ^2
      )
    })
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
      bias     = (input$abias[1]+input$abias[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^((input$ayadj[1]-6))+10^((input$ayadj[2]-6)))/2
    ))
  }
  if(input$mod %in% c("tal", "traab")){
    return(c(
      cov      = (input$akcov[1]+input$akcov[2])/2,
      bias     = (input$abias[1]+input$abias[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^(input$ayadj[1]-6)+10^(input$ayadj[2]-6))/2,
      diverg   = (input$adiv[1]+input$adiv[2])/2
    ))
  }
  if(input$mod == "tse"){
    return(c(
      cov      = (input$akcov[1]+input$akcov[2])/2,
      bias     = (input$abias[1]+input$abias[2])/2,
      theta    = (10^input$ath[1]+10^input$ath[2])/2,
      haplSize = (10^(input$ayadj[1]-6)+10^(input$ayadj[2]-6))/2,
      diverg   = (input$adiv[1]+input$adiv[2])/2,
      pallo    = (input$apallo[1]+input$apallo[2])/2
    ))
  }
}

doOptimisation <- function(input, sp, minFun, startingVals){
  if(input$mod %in% c("tal", "traab")){
    return(
      optim(startingVals,
            minFun,
            lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], 10^(input$ayadj[1]-6), input$adiv[1]),
            upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], 10^(input$ayadj[2]-6), input$adiv[2]),
            xlimits=c(input$axrange[1],input$axrange[2]),
            spec=sp,
            method = "L-BFGS-B"
      )
    )
  }
  if(input$mod == "tse"){
    return(
      optim(startingVals,
            minFun,
            lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], 10^(input$ayadj[1]-6), input$adiv[1], input$apallo[1]),
            upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], 10^(input$ayadj[2]-6), input$adiv[2], input$apallo[2]),
            xlimits=c(input$axrange[1],input$axrange[2]),
            spec=sp,
            method = "L-BFGS-B"
      )
    )
  }
  if(input$mod %in% c("d", "tau", "traaa")){
    return(
      optim(startingVals,
            minFun,
            lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], 10^(input$ayadj[1]-6)),
            upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], 10^(input$ayadj[2]-6)),
            xlimits=c(input$axrange[1],input$axrange[2]),
            spec = sp,
            method = "L-BFGS-B"
      )
    )
  }
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
  return(.sliderRanges)
}

#' Set slider ranges for Tetmer UI
#'
#' Updates the slider range settings used to initialise the Tetmer
#' Shiny interface. Call this before \code{tetmer()} to customise
#' the parameter ranges for your data.
#'
#' @param x A named list of slider range values in the same format
#'   as returned by \code{sliderRanges}
#' @return NULL, invisibly
#' @export
#' @examples
#' \dontrun{
#' ranges <- sliderRanges()
#' ranges$kcovMax <- 500
#' setSliderRanges(ranges)
#' }
setSliderRanges <- function(x){
  assignInMyNamespace(".sliderRanges", x)
}

#' Enable the segregating allotetraploid model
#'
#' Adds the segregating allotetraploid model (\code{"tse"}) to the list
#' of available models in the Tetmer UI. This model is hidden by default
#' as it is experimental. Call this function before launching the app
#' with \code{tetmer()}.
#'
#' @return NULL, invisibly
#' @export
#' @examples
#' \dontrun{
#' allowSegTet()
#' tetmer(E028)
#' }
allowSegTet <- function(){
  assignInMyNamespace("modelClasses",
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

  if(modelType %in% c("d", "tau", "traaa")){
    points(input$axrange[1]:input$axrange[2],
           colSums(eval(probs, envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                            tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3])))*optimised$par[4]*1000000,
           col="red", type='l', lty=1, lwd=2
    )
  }
  if(modelType %in% c("tal", "traab")){
    points(input$axrange[1]:input$axrange[2],
           colSums(eval(probs,
                        envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                     tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5]))
           )*optimised$par[4]*1000000,
           col="red", type='l', lty=1, lwd=2
    )
  }
  if(modelType == "tse"){
    points(input$axrange[1]:input$axrange[2],
           colSums(eval(probs,
                        envir = list(txmin=input$axrange[1], txmax=input$axrange[2],
                                     tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3], tdiverg=optimised$par[5], pal=optimised$par[6]))
           )*optimised$par[4]*1000000,
           col="red", type='l', lty=1, lwd=2
    )
  }

  return(new("spectrum",
             name=nam,
             data=data.frame(mult=integer(0), count=integer(0)),
             k=k))
}
