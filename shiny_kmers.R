library(shiny)

# Sepcify path to matrix file here:
#f <- "~/Dropbox/Euphrasia/kmer_model/Polyploid/comp0123-main.mx"
f <- "comp0167-main.mx"
# Then, execute everything below at once:


# Global variables can go here
xmax <- 200
ymax <- 10000000
kcov <- 40
bias <- 0.1
th <- 0.1
yadj <- 1000000
diverg <- 1




# Define the UI
ui <- fluidPage(titlePanel("K-mer explorer (by Hannes)"),
                fluidRow(
                  plotOutput('plot')
                ),
                fluidRow(
                  column(3,
                         wellPanel(
                           h4("1) make the fat black line fit the plotting area."),
                           numericInput('xmax', 'Max multiplicity', xmax),
                           numericInput('ymax', 'Max k-mer coverage', ymax)
                         )),
                  column(3,
                         wellPanel(
                           h4("2) make the black lines match."),
                           numericInput('kcov', 'Sequencing coverage', kcov),
                           numericInput('bias', 'Amount of noise', bias),
                           numericInput('th', 'theta', th),
                           numericInput('yadj', 'Adjust y-scale of expectation', yadj)
                         )),
                  column(3,
                         wellPanel(
                           h4("3) adjust the divergence to make the coloured lines match."),
                           numericInput('diverg', 'Divergence', diverg)
                         )),
                  column(3,
                         wellPanel(
                           h4("Other useful things"),
                           checkboxInput("transp.mat", label = "Swap samples", value = FALSE)
                         )
                  )
                )
                
)



# Define the server code
server <- function(input, output) {
  output$plot <- renderPlot({
    if(input$transp.mat)
      mat <- t(read.table(f, sep=" "))
    else
      mat <- read.table(f, sep=" ")
    t1 <- rowSums(mat)[2:1001]
    s1 <- rowSums(mat[,2:1001])[2:1001]
    p1 <- mat[2:1001, 1]
    
    t2 <- colSums(mat)[2:1001]
    s2 <- colSums(mat[2:1001,])[2:1001]
    p2 <- mat[1, 2:1001]
    
    plot(t1, xlim=c(0,input$xmax), ylim=c(0, input$ymax), type = "l", lwd=2,
         xlab="K-mer multiplicity", ylab="K-mer coverage")
    points(s1, col = 2, type="l")
    points(p1, col = 4, type="l")
    legend("topleft", col=c(0,1, 2, 4, 0, 1, 2, 4), ncol = 2, lwd=c(1,2, 1, 1, 1, 2, 1, 1),
           lty=c(1, 1, 1, 1, 2, 2, 2, 2), legend=c("Observed:","Total","shared","private",
                                                   "Fit:", "Total","shared","private")
           )
    distros <- t(sapply(1:2, function(peak){
      dnbinom(1:(input$xmax), input$kcov/input$bias*peak, mu = input$kcov * peak)
    }
    ))
    factort <- sapply(1:input$xmax, function(x) c(2-2/(1+input$th), 1/(1+input$th)))
    points(colSums(distros*factort)*input$yadj, col = 1, type = "l", lty=2, lwd=2)
    
    factorp <- sapply(1:input$xmax, function(x) {
      c(
        (2*input$th*(24+30*input$th+10*input$th^2+input$th^3))/
          ((1+input$th)*(2+input$th)*(6+input$th)^2),
        
        (input$th*(12+6*input$th+input$th^2))/
          ((1+input$th)*(2+input$th)*(6+input$th)^2)
      )
    }
    )
    split.priv = c(
      (2*exp(-input$diverg*(2+input$th))*input$th*(-12*(-2+input$th) -4*exp(input$diverg+input$diverg*input$th/2)*input$th^2*(6+input$th) +2*exp(2*input$diverg)*(6+input$th)^2 +exp(input$diverg*(2+input$th))*(6+input$th)^2*(-4+input$th^2))) /
        ((-2 + input$th)*(1 + input$th)*(2 + input$th)*(6 + input$th)^2),
      (exp(-input$diverg*(2+input$th))*(12*(-2+input$th)*input$th -8*exp(input$diverg+input$diverg*input$th/2)*input$th^2*(6+input$th) +4*exp(2*input$diverg)*(6+input$th)^2 +exp(input$diverg*(2+input$th))*(6+input$th)^2*(-4+input$th^2))) /
        ((-2 + input$th)*(1 + input$th)*(2 + input$th)*(6 + input$th)^2)
    )
        points(colSums(distros*split.priv)*input$yadj, col = 4, type = "l", lty=2, lwd=1)
    
    factors <- sapply(1:input$xmax, function(x){
      c(
        (4*input$th*(24+15*input$th+2*input$th^2))/
          ((1+input$th)*(2+input$th)*(6+input$th)^2),
        
        (8*(3+input$th^2))/
          ((1+input$th)*(2+input$th)*(6+input$th)^2)
      )
    }
    )
    split.shared = c(
      (4*exp(-input$diverg*(2+input$th))*input$th*(6*(-2+input$th) +2*exp(input$diverg+input$diverg*input$th/2)*input$th^2*(6+input$th) -exp(2*input$diverg)*(6+input$th)^2)) /
        ((-2 + input$th)*(1 + input$th)*(2 + input$th)*(6 + input$th)^2),
      -(4*exp(-input$diverg*(2+input$th))*(3*(-2+input$th)*input$th -2*exp(input$diverg+input$diverg*input$th/2)*input$th^2*(6+input$th) +exp(2*input$diverg)*(6+input$th)^2)) /
        ((-2 + input$th)*(1 + input$th)*(2 + input$th)*(6 + input$th)^2)
    )
    
    
    points(colSums(distros*split.shared)*input$yadj, col = 2, type = "l", lty=2, lwd=1)

  })
  
  
  
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)
