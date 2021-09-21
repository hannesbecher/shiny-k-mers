

# initial parameters for manual fitting ####
txmax <- 200
txmin <- 5
tymax <- 10000
tkcov <- 15
tbias <- 0.5
tth <- 0.04
tyadj <- 200
tdiverg <- 30

# starting ranges for auto fitting ####

agsl <- 1
agsh <- 2000
akcovl <- 10
akcovh <- 100
abiasl <- 0.1
abiash <- 3
athl <- -2
athh <- 0.6
adivl <- 0.1
adivh <- 100
axrangel <- 5
axrangeh <- 200


#' Run the Tetmer app
#'
#' @param input Default argument fo ar `shiny` server function. Leave empty.
#' @param output Default argument fo ar `shiny` server function. Leave empty.
#'
#' @return NULL
#' @importFrom stats optim
tet.server <- function(input, output) {



  output$plot <- renderPlot({

    probs <<- getProbs(input)
    factors <<- getFactors(input)

    # read spectrum
    spec <- prepare.spectrum(sp)

    if(input$fitmod=="man"){


      plotSpecApp(input, spec)
      addvertlines(input) # vertical lines, multiples of n
      pointsFit(input)
      output$outText <- textOut(input)
      return()

    }

    if(input$fitmod == "auto"){

      plotSpecApp(input, spec)
      minFun <<- makeMinFun(input)
      startingVals <<- getStartingVals(input)
      optimised <- doOptimisation(input)
      addvertlines(input, optimised)
      pointsFit(input, optimised)
      output$outText <- textOut(input, optimised)

    }# if auto

  })
}


# The user interface
tet.ui <- fluidPage(titlePanel("Tetmer"),
                    "Fitting paramters to k-mer spectra (by Hannes Becher)",
                    fluidRow(
                      column(8, plotOutput('plot')),
                      column(4,
                             verbatimTextOutput("outText")
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
                                       radioButtons("mod", "Model",
                                                    c("Tetraploid (aabb)" = "tal",
                                                      "Tetraploid (aaaa)" = "tau",
                                                      "Diploid" = "d",
                                                      "Triploid (aaa)" = "traaa",
                                                      "Triploid (aab)" = "traab"
                                                    )
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
                                                h4("3rd: Adjust parameters, make the fit match the data"),
                                                numericInput('tkcov', 'k-mer  multiplicity', tkcov),
                                                numericInput('tbias', 'Peak width', tbias),
                                                numericInput('tth', 'theta', tth),
                                                numericInput('tyadj', 'Haploid non-rep GS (Mbp)', tyadj)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'man') && (['tal', 'traab'].includes(input.mod))",
                                              wellPanel(h4("4th: Only allopolyploids, adjust sub-genome split time"),
                                                        numericInput('tdiverg', 'T (in units of 2Ne)', tdiverg)
                                              )
                             )),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("2nd: Adjust the fitting area, make all data peaks visible"),
                                                        sliderInput("axrange", "x limits for fitting",
                                                                    min=0, max = 500,
                                                                    value=c(axrangel,axrangeh)),
                                                        sliderInput("ymax", "y axis max (does not affect fit)",
                                                                    min=0, max = 100,
                                                                    value=10)
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.fitmod == 'auto'",
                                              wellPanel(h4("3rd: adjust parameter ranges for fitting"),

                                                        sliderInput('akcov', 'k-mer  multiplicity',
                                                                    min=5, max = 200,
                                                                    value=c(akcovl, akcovh)),
                                                        sliderInput('abias', 'Peak width',
                                                                    min=0.01, max = 4,
                                                                    value=c(abiasl, abiash)),
                                                        sliderInput('ath', "log10 of theta",
                                                                    min=-4, max = 1, step = 0.05,
                                                                    value=c(athl, athh)),
                                                        sliderInput('ayadj', 'Haploid non-rep GS (Mbp)',
                                                                    min=1, max =2000,
                                                                    value=c(agsl, agsh))
                                              ))),
                      column(3,
                             conditionalPanel(condition = "(input.fitmod == 'auto') && (['tal', 'traab'].includes(input.mod))",
                                              wellPanel(h4("4th: Allopolyploids only, adjst sub-genome split time"),
                                                        sliderInput('adiv', 'T (in units of 2Ne)',
                                                                    min=0.001, max=100,
                                                                    value=c(adivl, adivh))
                                              )
                             ))
                    )

)


#' Run the interactive Tetmer app
#'
#' @param sp A `spectrum` object as generate by `read.spectrum`.
#'
#' @return NULL
#' @export
#'
#' @examples \dontrun{tetmer(E028)}
#' \dontrun{tetmer(E030)}
tetmer <- function(sp){
  sp <<- sp
  shinyApp(ui = tet.ui, server = tet.server)
}

#' A named k-mer spectrum class
#'
#' @slot name A string indicating the name of the spectrum (or sample)
#'
#' @slot data A `data.frame` with numeric cloumns mult and count.
#' @slot k A `numeric` indiciating the k-mer length.
#'
#' @export
setClass("spectrum", slots=list(name="character", data="data.frame", k="numeric"))



#' Read in a k-mer spectrum
#'
#' Uses `read.table()` to read in a k-mer spectrum.
#'
#'
#' @param f A string indicating the path to th spectrum file
#' @param nam (otional) A string indicating the name of the spectrum
#' @param k (otional) A numeric indicating the k-mer length
#' @param ... keyword arguments to be passed to `read.table()`
#' @return A two-element list comprising a names string and a two-column data frame of the
#' k-mer spectrum
#' @export
#' @importFrom methods new
#' @importFrom utils read.table
#' @examples
#' testdf <- data.frame(mult=1:100, count = 100:1)
#' tf <- tempfile()
#' write.table(testdf, tf)
#' sp <- read.spectrum(tf, nam = "Test spectrum")
#' unlink(tf)
#' plot(sp)
read.spectrum <- function(f,
                          nam="MySpectrum",
                          k=0,
                          ...){
  sp = read.table(f, ...)
  names(sp) = c("mult", "count")
  return(new("spectrum",
             name=nam,
             data=sp,
             k=k
  )
  )
}



#' Plot a k-mer spectrum
#'
#' Uses the plot function to generate a scatter plot of a k-mer spectrum.
#'
#' It may be useful to set `log="xy"` or to limit the plotting range using `xlims` or `ylims`.
#'
#' @param x An object of class spectrum
#' @param main (optional) A string passed to `plot()`
#' @param xlab (optional) A two-element numeric to `plot()`
#' @param ylab (optional) A two-element numeric passed to `plot()`
#' @param ... other keyword arguments to be passed to `plot()`
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
#' Using `input`, this function plots dashed vertical lines at multiples of `n`
#' corresponding to the samples's ploidy.
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
    if(input$mod %in% c("tau", "tal")){
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
      text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
    }
    if(input$mod%in%c("traaa", "traab")){
      abline(v=c(1,2,3)*optimised$par[1], lty=2,col="grey")
      text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
      text(optimised$par[1]*3, input$ymax*.7*1000000, "3x", col="grey", cex=4)

    }
    if(input$mod%in%c("tau", "tal")){
      abline(v=c(1,2,3,4)*optimised$par[1], lty=2,col="grey")
      text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
      text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
      text(optimised$par[1]*3, input$ymax*.7*1000000, "3x", col="grey", cex=4)
      text(optimised$par[1]*4, input$ymax*.7*1000000, "4x", col="grey", cex=4)
    }
  } # if auto mode
}

#' Plot annotated spectrum
#'
#' Wrapper around `plot.spectrum` adding some other bits.
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

  }
  if(input$fitmod=="auto"){
    if(input$showData){
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, input$ymax*1000000))
    } else {
      plot(spec, xlim=c(0,input$axrange[2]), ylim=c(0, input$ymax*1000000),
           type = 'n')
    }
    abline(v=input$axrange)

  }
  legend("topright", col=c(1,2,2), lwd=c(2, 2, 1), lty=c(0,1, 2), pch=c(1,NA, NA), legend=c("Data","Fit", "Extrapolation"))


}



#' Plot the fit
#'
#'
#' @param input From GUI
#' @param optimised From fit (only supplied in auto mode)
#'
#' @keywords internal
#'
#' @return NULL
pointsFit <- function(input, optimised=0){

  if(input$fitmod=="man"){
    if(input$mod %in% c("d", "tau", "traaa")){
    points(
      colSums(eval(probs, envir = list(txmin= 1, txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                eval(factors, envir=list(tth=input$tth)))*input$tyadj*1000000,
      col="red", type = 'l', lty=1, lwd=2
    )
    }
    if(input$mod %in% c("tal", "traab")){
      points(
        colSums(eval(probs, envir = list(txmin= 1, txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                  eval(factors, envir=list(tth=input$tth, tdiverg=input$tdiverg)))*input$tyadj*1000000,
        col="red", type = 'l', lty=1, lwd=2
      )
    }
  }
  if(input$fitmod=="auto"){
    if(input$mod %in% c("d", "tau", "traaa")){
    points(input$axrange[1]:input$axrange[2],
           colSums(eval(probs, envir = list(txmin=input$axrange[1], txmax=input$axrange[2], tkcov=optimised$par[1], tbias=optimised$par[2])) *
                     eval(factors, envir=list(tth=optimised$par[3])))*optimised$par[4]*1000000,
           col="red", type = 'l', lty=1, lwd=2
    )
    }
    if(input$mod %in% c("tal", "traab")){
      points(input$axrange[1]:input$axrange[2],
             colSums(eval(probs,
                          envir = list(txmin=input$axrange[1],
                                       txmax=input$axrange[2],
                                       tkcov=optimised$par[1],
                                       tbias=optimised$par[2]
                                       )) *
                       eval(factors,
                            envir=list(tth=optimised$par[3],
                                       tdiverg=optimised$par[5]))
                     )*optimised$par[4]*1000000,
             col="red", type = 'l', lty=1, lwd=2

      )
      points(0:input$axrange[1],
             colSums(eval(probs,
                          envir = list(txmin=0,
                                       txmax=input$axrange[1],
                                       tkcov=optimised$par[1],
                                       tbias=optimised$par[2]
                          )) *
                       eval(factors,
                            envir=list(tth=optimised$par[3],
                                       tdiverg=optimised$par[5]))
             )*optimised$par[4]*1000000,
             col="red", type = 'l', lty=2, lwd=1
      )


    }
  }
}


#' Generate text for Tetmer window
#'
#' @param input Input from the GUI (contains plotting range, model, etc.)
#' @param optimised Fitted values from `optim`.
#'
#' @return A string of the fitted parameters produced by \code{renderText()}
#' to be displayed in the Tetmer window.
#' @keywords internal
textOut <- function(input, optimised=0){
  if(sp@k > 0){
    if(input$fitmod == "man"){
      if(input$mod=="d"){
        return(renderText(paste(
          "DIPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
        )
      }
      if(input$mod=="tau"){
        return(
          renderText(paste(
            "AUTOTETRAPLOID MODEL, MANUAL FIT",
            "\n    haploid k-mer cov:", input$tkcov,
            "\n      per k-mer theta:", input$tth,
            "\nhapl non-rep GS (Mbp):", input$tyadj,
            "\n    bias (peak width):", input$tbias
          ))
        )
      }
      if(input$mod=="traaa"){
        return(
          renderText(paste(
            "AUTOTRIPLOID MODEL, MANUAL FIT",
            "\n    haploid k-mer cov:", input$tkcov,
            "\n      per k-mer theta:", input$tth,
            "\nhapl non-rep GS (Mbp):", input$tyadj,
            "\n    bias (peak width):", input$tbias
          ))
        )
      }
      if(input$mod=="traab"){
        return(
          renderText(paste(
            "ALLOTRIPLOID MODEL, MANUAL FIT",
            "\n    haploid k-mer cov:", input$tkcov,
            "\n      per k-mer theta:", input$tth,
            "\n                    T:", input$tdiverg,
            "\nhapl non-rep GS (Mbp):", input$tyadj,
            "\n    bias (peak width):", input$tbias
          ))
        )
      }
      if(input$mod=="tal"){
        return(
          renderText(paste(
            "ALLOTETRAPLOID MODEL, MANUAL FIT",
            "\n    haploid k-mer cov:", input$tkcov,
            "\n      per k-mer theta:", input$tth,
            "\n                    T:", input$tdiverg,
            "\nhapl non-rep GS (Mbp):", input$tyadj,
            "\n    bias (peak width):", input$tbias
          ))
        )
      }
    } # if man
    if(input$fitmod == "auto"){
      if(input$mod=="tal"){
        return(
          renderText(paste(
            "ALLOTETRAPLOID MODEL, AUTO FITTED",
            "\n    haploid k-mer cov:", round(optimised$par[1],1),
#            "\n      per k-mer theta:", round(optimised$par[3],4),
            "\n per nucleotide theta:", round(optimised$par[3] / sp@k,5),
            "\n                    T:", round(optimised$par[5],2),
            "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],1),
            "\nper nucleotide diverg:", round(optimised$par[3]*optimised$par[5]/sp@k, 4),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n         k-mer length:", sp@k,
            "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
            "\n                    T:", input$adiv[1], input$adiv[2],
            "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))
        )
      } # if allotet
      if(input$mod=="traab"){
        return(
          renderText(paste(
            "ALLOTRIPLOID MODEL, AUTO FITTED",
            "\n    haploid k-mer cov:", round(optimised$par[1],1),
#            "\n      per k-mer theta:", round(optimised$par[3],4),
            "\n per nucleotide theta:", round(optimised$par[3] / sp@k,5),
            "\n                    T:", round(optimised$par[5],2),
            "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],1),
#            "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
            "\nper nucleotide diverg:", round(optimised$par[3]*optimised$par[5], 4),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n         k-mer length:", sp@k,
            "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
            "\n                    T:", input$adiv[1], input$adiv[2],
            "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))
        )
      } # if allotrip
      if(input$mod=="tau"){
        return(
          renderText(paste(
            "AUTOTETRAPLOID MODLE, AUTO FITTED",
            "\n    haploid k-mer cov:", round(optimised$par[1],1),
#            "\n      per k-mer theta:", round(optimised$par[3],4),
            "\n per nucleotide theta:", round(optimised$par[3] / sp@k,5),
            "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],1),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n         k-mer length:", sp@k,
            "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
            "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))
        )
      } # if tau
      if(input$mod=="traaa"){
        return(
          renderText(paste(
            "AUTOTRIPLOID MODLE, AUTO FITTED",
            "\n    haploid k-mer cov:", round(optimised$par[1],1),
#            "\n      per k-mer theta:", round(optimised$par[3],4),
            "\n per nucleotide theta:", round(optimised$par[3] / sp@k,5),
            "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],1),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n         k-mer length:", sp@k,
            "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
            "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))
        )
      } # if traab
      if(input$mod=="d"){
        return(
          renderText(paste(
            "DIPLOID MODEL, AUTO FITTED",
            "\n    haploid k-mer cov:", round(optimised$par[1],1),
#            "\n      per k-mer theta:", round(optimised$par[3],4),
            "\n per nucleotide theta:", round(optimised$par[3] / sp@k,5),
            "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
            "\n    bias (peak width):", round(optimised$par[2],1),
            "\n\nSTARTING RANGES (MIN MAX)",
            "\n         k-mer length:", sp@k,
            "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
            "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
            "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
            "\n    bias (peak width):", input$abias[1], input$abias[2],
            "\n              x range:", input$axrange[1],input$axrange[2]
          ))
        )
      } #if d
    } # if auto fit
  } # if k > 0
  if(input$fitmod == "man"){
    if(input$mod=="d"){
      return(renderText(paste(
        "DIPLOID MODEL, MANUAL FIT",
        "\n    haploid k-mer cov:", input$tkcov,
        "\n      per k-mer theta:", input$tth,
        "\nhapl non-rep GS (Mbp):", input$tyadj,
        "\n    bias (peak width):", input$tbias
      ))
      )
    }
    if(input$mod=="tau"){
      return(
        renderText(paste(
          "AUTOTETRAPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      )
    }
    if(input$mod=="traaa"){
      return(
        renderText(paste(
          "AUTOTRIPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      )
    }
    if(input$mod=="traab"){
      return(
        renderText(paste(
          "ALLOTRIPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\n                    T:", input$tdiverg,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      )
    }
    if(input$mod=="tal"){
      return(
        renderText(paste(
          "ALLOTETRAPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\n                    T:", input$tdiverg,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
      )
    }
  } # if man
  if(input$fitmod == "auto"){
    if(input$mod=="tal"){
      return(
        renderText(paste(
          "ALLOTETRAPLOID MODEL, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
          "\n                    T:", round(optimised$par[5],2),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
      )
    } # if allotet
    if(input$mod=="traab"){
      return(
        renderText(paste(
          "ALLOTRIPLOID MODEL, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
          "\n                    T:", round(optimised$par[5],2),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[5], 3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
      )
    } # if allotrip
    if(input$mod=="tau"){
      return(
        renderText(paste(
          "AUTOTETRAPLOID MODLE, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
      )
    } # if tau
    if(input$mod=="traaa"){
      return(
        renderText(paste(
          "AUTOTRIPLOID MODLE, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
      )
    } # if traab
    if(input$mod=="d"){
      return(
        renderText(paste(
          "DIPLOID MODEL, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
#          "\n      per nucleotide theta:", round(optimised$par[3] / sp@k,4),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[4],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
      )
    } #if d
  } # if auto fit
}



#' Prepare spectrum for plotting and fitting
#'
#' @param spe A `spectrum` object
#'
#' @keywords internal
#'
#' @return A `spectrum` object.
prepare.spectrum <- function(spe){
  sp <- spe
  # cut off first row if 0-based spectrum

  if(sp@data[1,1] == 0) {
    # print("Cutting off multiplicity zero.")
    sp@data <- sp@data[2:nrow(sp@data),]
  }
  # add lines if spectrum starts with multiplicity > 1
  offs = sp@data[1,1] - 1
  if(offs > 0){
    # print("Padding with zeros to fill in lower multiplicities.")
    prepDF = data.frame(mult=1:offs, count=rep(0, offs))
    sp@data <- rbind(prepDF, sp@data)
  }
  return(sp)
}


#' Generate a function to be minimised
#'
#' This is then mimimised by \code{optim}.
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#'
#' @return A function servoing as input to \code{optim}
#'
#' @keywords internal
makeMinFun <- function(input){

  if(input$mod=="d"){
    # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - GS in Mbp

    return(function(x, xlimits, spec) {
      sum(

        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3])))*x[4]*1000000) ^2

      )
    })
  }
  if(input$mod %in% c("tal", "traab")){
    # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - GS in Mbp, 5 - divergence
    return(function(x, xlimits, spec) {
      sum(

        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3], tdiverg=x[5])))*x[4]*1000000) ^2

      )
    }
    )
  }
  if(input$mod %in% c("tau", "traaa")){
    # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - GS in Mbp
    return(function(x,xlimits, spec) {
      sum(

        (spec@data$count[xlimits[1]:xlimits[2]] -
           colSums(eval(probs, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                     eval(factors, envir=list(tth=x[3])))*x[4]*1000000) ^2

      )
    }
    )
  }
}

#' Starting values for optimisation
#'
#' @param input Input from GUI (contains model, parameter ranges, etc.)
#'
#' @return A named vector of starting vlues
#' @keywords internal
getStartingVals <- function(input){
  if(input$mod%in% c("d", "tau", "traaa")){
    return(
      # kmer coverage, bias, theta, genome size
      c(
        cov=(input$akcov[1]+input$akcov[2])/2,
        bias=(input$abias[1]+input$abias[2])/2,
        theta=(10^input$ath[1]+10^input$ath[2])/2,
        haplSize=(input$ayadj[1]+input$ayadj[2])/2
      )
    )
  }
  if(input$mod %in% c("tal", "traab")){
    return(
      # kmer coverage, bias, theta, genome size, divergence
      c(
        cov=(input$akcov[1]+input$akcov[2])/2,
        bias=(input$abias[1]+input$abias[2])/2,
        theta=(10^input$ath[1]+10^input$ath[2])/2,
        haplSize=(input$ayadj[1]+input$ayadj[2])/2,
        diverg=(input$adiv[1]+input$adiv[2])/2
      )
    )
  }
}

doOptimisation <- function(input){
  if(input$mod %in% c("tal", "traab")){

    return(
      optim(startingVals, # starting values (vector of)
                       minFun,
                       lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], input$ayadj[1], input$adiv[1]),
                       upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], input$ayadj[2], input$adiv[2]),
                       xlimits=c(input$axrange[1],input$axrange[2]),
                       spec=sp,
                       method = "L-BFGS-B"
    )
    )

  } # if tal
  if(input$mod %in% c("d", "tau", "traaa")){
    return(
      optim(startingVals, # starting values (vector of)
                       minFun,
                       lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], input$ayadj[1]),
                       upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], input$ayadj[2]),
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
  if(input$mod=="tau"){
    return(factorAut)
  }
  if(input$mod=="tal"){
    return(factorAll)
  }
  if(input$mod=="traaa"){
    return(factorTraaa)
  }
  if(input$mod=="traab"){
    return(factorTraab)
  }
  if(input$mod=="d"){
    return(factorDip)
  }
}

#' Raw k-mer spectrum peaks
#'
#' @param input Input from GUI (contains model etc.)
#'
#' @return An expression corresponding to the size and number of k-mer spectrum peaks
#'
#' @keywords internal
getProbs <- function(input){
  if(input$mod %in% c("tau", "tal")){
    return(probsTet)
  }
  if(input$mod %in% c("traaa", "traab")){
    return(probsTri)
  }
  if(input$mod=="d"){
    return(probsDip)
  }
}

