library(shiny)


# initial parameters for manual fitting
txmax <- 100
txmin <- 5
tymax <- 10000
tkcov <- 15
tbias <- 0.5
tth <- 0.04
tyadj <- 200
tdiverg <- 30

# starting ranges for auto fitting

agsl <- 20
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



# Define the UI
tet.ui <- fluidPage(titlePanel("Tetmer"),
                    "Fitting paramters to allotetraploid k-mer spectra (by Hannes Becher)",
                    fluidRow(
                      column(8, plotOutput('plot')),
                      column(4,
                             verbatimTextOutput("outText")
                      )
                    ),
                    fluidRow(
                      column(3,
                             wellPanel(h4("1st: Select fitting mode and model"),
                                       radioButtons("fitmod", "Fitting mode",
                                                    c("Manual" = "man",
                                                      "Autofit" = "auto")
                                       ),
                                       radioButtons("mod", "Model",
                                                    c("Tetraploid, allo" = "tal",
                                                      "Tetraploid, auto" = "tau",
                                                      "Diploid" = "d"
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
                             conditionalPanel(condition = "input.mod == 'tal' && input.fitmod == 'man'",
                                              wellPanel(h4("4th: Only allotetraploids, adjust sub-genome split time"),
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
                                                                    min=10, max =2000,
                                                                    value=c(agsl, agsh))
                                              ))),
                      column(3,
                             conditionalPanel(condition = "input.mod == 'tal' && input.fitmod == 'auto'",
                                              wellPanel(h4("4rd: Allotetraploids only, adjst sub-genome split time"),
                                                        sliderInput('adiv', 'T (in units of 2Ne)',
                                                                    min=0.001, max=100,
                                                                    value=c(adivl, adivh))
                                              )
                             ))
                    )
                    
)

tet.server <- function(input, output) {
  
  output$plot <- renderPlot({
    
    
    # read spectrum
    counts <- read.table(specPath, skip=6, col.names = c("Muliplicity","Count"))
    
    probsDip <- expression(rbind(
      dnbinom(txmin:txmax, tkcov/tbias*1, mu = tkcov * 1),
      dnbinom(txmin:txmax, tkcov/tbias*2, mu = tkcov * 2)
    )
    )
    
    probsTet <- expression(rbind(
      dnbinom(txmin:txmax, tkcov/tbias*1, mu = tkcov * 1),
      dnbinom(txmin:txmax, tkcov/tbias*2, mu = tkcov * 2),
      dnbinom(txmin:txmax, tkcov/tbias*3, mu = tkcov * 3),
      dnbinom(txmin:txmax, tkcov/tbias*4, mu = tkcov * 4)
    ))
    
    factorDip <- expression(
      c(
        2*tth/(1+tth),
        1/(1+tth)
      )
    )
    factorAut <- expression(
      c(4*tth/(3+tth),
        6*tth/(6+5*tth+tth^2),
        8*tth/(6+11*tth+6*tth^2+tth^3),
        6/(6+11*tth+6*tth^2+tth^3))
    )
    factorAll <- expression(
      c(
        4*exp(-3*tth*tdiverg)*tth*(
          -2*exp(tdiverg*(tth-2)) +
            exp(3*tdiverg*tth) * (2+tth)^2 * (3 + tth) -
            2*exp(2*tth*tdiverg)*(3+4*tth+tth^2)
        ),
        2*exp(-1/2*tdiverg*(4+5*tth))*(
          6*exp(tdiverg*tth/2)*tth+
            exp(1/2*tdiverg*(4+5*tth)) * (2+tth)^2 * (3 + tth) +
            2*exp(2*tdiverg+3/2*tdiverg*tth)*
            (-6-8*tth+tth^2+tth^3)
        ),
        8*exp(-tdiverg*(3+4*tth))*tth*(
          -exp(tdiverg+2*tdiverg*tth)+
            exp(3+tdiverg*(1+tth))*(3+tth)
        ),
        2*exp(-2*tdiverg*(1+tth))*(
          tth + 2*exp(tdiverg*(2+tth))*(3+tth)
        )
      )/ (1+tth) / (2+tth)^2 / (3+tth)
      
    )
    
    if(input$fitmod=="man"){
      
      # plot spectrum
      plot(counts[,1],counts[,2], xlim=c(0,input$txmax), ylim=c(0, input$tymax*1000),
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count", main = specPath)
      legend("topleft", col=c(1,2), lwd=2, lty=c(0,2), pch=c(1,NA), legend=c("Data","Fit"))
      
      if(input$mod =="d") {
        abline(v=c(1,2)*input$tkcov, lty=2,col="grey")
        text(input$tkcov[1], input$tymax*.7*1000, "1x", col="grey", cex=4)
        text(input$tkcov[1]*2, input$tymax*.7*1000, "2x", col="grey", cex=4)
      } else {
        abline(v=c(1,2,3,4)*input$tkcov, lty=2,col="grey")
        text(input$tkcov[1], input$tymax*.7*1000, "1x", col="grey", cex=4)
        text(input$tkcov[1]*2, input$tymax*.7*1000, "2x", col="grey", cex=4)
        text(input$tkcov[1]*3, input$tymax*.7*1000, "3x", col="grey", cex=4)
        text(input$tkcov[1]*4, input$tymax*.7*1000, "4x", col="grey", cex=4)
      }
      
      
      if(input$mod=="d"){
        points(
          colSums(eval(probsDip, envir = list(txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                    eval(factorDip, envir=list(tth=input$tth)))*input$tyadj*1000000,
          col="red", type = 'l', lty=1, lwd=2
        )
        output$outText <- renderText(paste(
          "DIPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
        return()
      }
      if(input$mod=="tau"){
        points(
          colSums(eval(probsTet, envir = list(txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                    eval(factorAut, envir=list(tth=input$tth)))*input$tyadj*1000000,
          col="red", type = 'l', lty=1, lwd=2
        )
        output$outText <- renderText(paste(
          "AUTOTETRAPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
        return()
      }
      if(input$mod=="tal"){
        points(
          colSums(eval(probsTet, envir = list(txmax=input$txmax, tkcov=input$tkcov, tbias=input$tbias)) *
                    eval(factorAll, envir=list(tth=input$tth, diverg=input$tdiverg)))*input$tyadj*1000000,
          col="red", type = 'l', lty=1, lwd=2
        )
        output$outText <- renderText(paste(
          "ALLOTETRAPLOID MODEL, MANUAL FIT",
          "\n    haploid k-mer cov:", input$tkcov,
          "\n      per k-mer theta:", input$tth,
          "\n                    T:", input$tdiverg,
          "\nhapl non-rep GS (Mbp):", input$tyadj,
          "\n    bias (peak width):", input$tbias
        ))
        return()
      }
    }
    
    if(input$fitmod == "auto"){
      
      plot(counts[,1],counts[,2], xlim=c(0,input$axrange[2]), ylim=c(0, input$ymax*1000000),
           xlab="K-mer multiplicity (coverage)", ylab="K-mer count", main = specPath)
      
      
      
      # show limits
      abline(v=input$axrange)        
      legend("topleft", col=c(1,2), lwd=2, lty=c(0,2), pch=c(1, NA), legend=c("Data","Fit"))
      
      if(input$mod=="tal"){
        
        # function to be minimised
        # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - divergence, 5 - GS in Mbp
        minAllo <- function(x,xlimits) {
          sum(
            (counts[xlimits[1]:xlimits[2],2] -
               colSums(eval(probsTet, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                         eval(factorAll, envir=list(tth=x[3], tdiverg=x[4])))*x[5]*1000000) ^2)
        }
        
        # function to plot the fit
        pointsFit <- function(x,xlimits){
          points(xlimits[1]:xlimits[2],
                 colSums(eval(probsTet, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                           eval(factorAll, envir=list(tth=x[3], tdiverg=x[4])))*x[5]*1000000,
                 col="red", type = 'l', lty=1, lwd=2
          )
        }
        
        # kmer coverage, bias, theta, divergence, genome size
        startingVals <- c(
          cov=(input$akcov[1]+input$akcov[2])/2,
          bias=(input$abias[1]+input$abias[2])/2,
          theta=(10^input$ath[1]+10^input$ath[2])/2,
          diverg=(input$adiv[1]+input$adiv[2])/2,
          haplSize=(input$ayadj[1]+input$ayadj[2])/2
        )
        optimised <- optim(startingVals, # starting values (vector of)
                           minAllo,
                           lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], input$adiv[1], input$ayadj[1]),
                           upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], input$adiv[2], input$ayadj[2]),
                           xlimits=c(input$axrange[1],input$axrange[2]),
                           method = "L-BFGS-B"
        )
        print(optimised$par)
        abline(v=c(1,2,3,4)*optimised$par[1], lty=2,col="grey")
        text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
        text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
        text(optimised$par[1]*3, input$ymax*.7*1000000, "3x", col="grey", cex=4)
        text(optimised$par[1]*4, input$ymax*.7*1000000, "4x", col="grey", cex=4)
        
        pointsFit(optimised$par,input$axrange)
        
        
        output$outText <- renderText(paste(
          "ALLOTETRAPLOID MODEL, AUTO FITTED",
          "\n    haploid k-mer cov:", round(optimised$par[1],1),
          "\n      per k-mer theta:", round(optimised$par[3],4),
          "\n                    T:", round(optimised$par[4],2),
          "\nhapl non-rep GS (Mbp):", round(optimised$par[5],1),
          "\n    bias (peak width):", round(optimised$par[2],1),
          "\n     per k-mer diverg:", round(optimised$par[3]*optimised$par[4], 3),
          "\n\nSTARTING RANGES (MIN MAX)",
          "\n    haploid k-mer cov:", input$akcov[1], input$akcov[2],
          "\nlog10 per k-mer theta:", input$ath[1], input$ath[2],
          "\n                    T:", input$adiv[1], input$adiv[2],
          "\nhapl non-rep GS (Mbp):", input$ayadj[1], input$ayadj[2],
          "\n    bias (peak width):", input$abias[1], input$abias[2],
          "\n              x range:", input$axrange[1],input$axrange[2]
        ))
        return()
        
      } # if tal
      if(input$mod=="tau"){
        # function to be minimised
        # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - GS in Mbp
        minAuto <- function(x,xlimits) {
          sum(
            (counts[xlimits[1]:xlimits[2],2] -
               colSums(eval(probsTet, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                         eval(factorAut, envir=list(tth=x[3])))*x[4]*1000000) ^2)
        }
        
        # function to plot the fit
        pointsFit <- function(x,xlimits){
          points(xlimits[1]:xlimits[2],
                 colSums(eval(probsTet, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                           eval(factorAut, envir=list(tth=x[3])))*x[4]*1000000,
                 col="red", type = 'l', lty=1, lwd=2
          )
        }
        
        
        # kmer coverage, bias, theta, divergence, genome size
        startingVals <- c(
          cov=(input$akcov[1]+input$akcov[2])/2,
          bias=(input$abias[1]+input$abias[2])/2,
          theta=(10^input$ath[1]+10^input$ath[2])/2,
          haplSize=(input$ayadj[1]+input$ayadj[2])/2
        )
        optimised <- optim(startingVals, # starting values (vector of)
                           minAuto,
                           lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], input$ayadj[1]),
                           upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], input$ayadj[2]),
                           xlimits=c(input$axrange[1],input$axrange[2]),
                           method = "L-BFGS-B"
                           #      method = "Nelder-Mead"
        )
        abline(v=c(1,2,3,4)*optimised$par[1], lty=2,col="grey")
        text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
        text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
        text(optimised$par[1]*3, input$ymax*.7*1000000, "3x", col="grey", cex=4)
        text(optimised$par[1]*4, input$ymax*.7*1000000, "4x", col="grey", cex=4)
        
        pointsFit(optimised$par,input$axrange)
        output$outText <- renderText(paste(
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
        return()
      }
      if(input$mod=="d"){
        # function to be minimised
        # parameter vector: 1 - k-mer cov, 2 - peak width, 3 - theta, 4 - GS in Mbp
        minDip <- function(x,xlimits) {
          sum(
            (counts[xlimits[1]:xlimits[2],2] -
               colSums(eval(probsDip, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                         eval(factorDip, envir=list(tth=x[3])))*x[4]*1000000) ^2)
        }
        
        # function to plot the fit
        pointsFit <- function(x,xlimits){
          points(xlimits[1]:xlimits[2],
                 colSums(eval(probsDip, envir = list(txmin=xlimits[1], txmax=xlimits[2], tkcov=x[1], tbias=x[2])) *
                           eval(factorDip, envir=list(tth=x[3])))*x[4]*1000000,
                 col="red", type = 'l', lty=1, lwd=2
          )
        }
        
        
        # kmer coverage, bias, theta, divergence, genome size
        startingVals <- c(
          cov=(input$akcov[1]+input$akcov[2])/2,
          bias=(input$abias[1]+input$abias[2])/2,
          theta=(10^input$ath[1]+10^input$ath[2])/2,
          haplSize=(input$ayadj[1]+input$ayadj[2])/2
        )
        print("Diploid starting vals:")
        print(startingVals)
        optimised <- optim(startingVals, # starting values (vector of)
                           minDip,
                           lower=c(input$akcov[1], input$abias[1], 10^input$ath[1], input$ayadj[1]),
                           upper=c(input$akcov[2], input$abias[2], 10^input$ath[2], input$ayadj[2]),
                           xlimits=c(input$axrange[1],input$axrange[2]),
                           method = "L-BFGS-B"
                           #      method = "Nelder-Mead"
        )
        print("Diploid results:")
        print(optimised$par)
        abline(v=c(1,2)*optimised$par[1], lty=2,col="grey")
        text(optimised$par[1], input$ymax*.7*1000000, "1x", col="grey", cex=4)
        text(optimised$par[1]*2, input$ymax*.7*1000000, "2x", col="grey", cex=4)
        
        
        pointsFit(optimised$par,input$axrange)
        
        
        output$outText <- renderText(paste(
          "DIPLOID MODEL, AUTO FITTED",
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
        
        
        
        return()
      }
    }# if auto
    
  })
}

# Run all the code above ####

# Set path to directory with k-mer spectra:
setwd("~/git_repos/shiny-k-mers/data/")

 
# Set "specPath" to the name of th espectrum file you want to fit:

# # diplods
# specPath <- 'AN'
# specPath <- 'RO'
# specPath <- 'RI' 
# specPath <- 'VI'
# 
# # tetraploids
# specPath <- 'A0'
# specPath <- 'A1'
# specPath <- 'A2'
# specPath <- 'A3'
# specPath <- 'F1'
# specPath <- 'F2'
# specPath <- 'F3'
# specPath <- 'F4'
# specPath <- 'M0'
# specPath <- 'M1'
# specPath <- 'M2'
# specPath <- 'M3'
# specPath <- 'X1'
# specPath <- 'X2'
  

# Return a Shiny app object
shinyApp(ui = tet.ui, server = tet.server)

