install.packages("/Users/hbecher/git_repos/shiny-k-mers/Tetmer_2.3.2.tar.gz",
                 repo=NULL)


roxygen2::roxygenise(clean=TRUE)
library(Tetmer)
browseVignettes("Tetmer")



mySpec <- tetmer(E028)
mySpec <- tetmer(mySpec)
plot(mySpec, fitIndex = 1)
plot(mySpec, fitIndex = 2)

?plot.spectrum
plot(mySpec, residuals=TRUE)

plot(mySpec, residuals=TRUE, fitIndex = 2)

plot(mySpec, residuals="overlay", col=4, colFit=2, colResid=3)

plot(mySpec, residuals="overlay", fitIndex = 2)

plot(E028, residuals=TRUE)
plot(E028, residuals="overlay")

plot(1:10, type="l", col="orange", lwd=2)
?plot.spectrum
?expand.grid


getwd()
?write.spectrum
write.spectrum(mySpec, file="~/temp/mySpec.tsp")
write.spectrum(E030, file="~/temp/E030.tsp")
plot(E030, log="xy")



ee030 <- read.spectrum("~/temp/E030.tsp")
plot(ee030, log="xy")


msp <- read.spectrum("~/temp/mySpec.tsp")
plot(mySpec)
plot(msp)
msp@fits
mySpec@fits


# Non-interactive ---------------------------------------------------------




?fitSpectrum
mySpec@fits
fit2 <- fitSpectrum(E028,
                    model  = "tal",
                    kcov   = c(35, 45),
                    vf     = c(1, 5),
                    log10theta  = c(-3, 0),
                    log10Mbp =  c(0, 3),
                    xrange = c(30, 200))
fit2@par
plot(fit2)
predict(fit2)
msp2 <- addFit(E028, fit2)
plot(msp2, residuals = "overlay")
