library(usethis)
E030 <- read.spectrum("../data/E030full.hist.no0", "E. anglica, E030", 21)
E028 <- read.spectrum("../data/A0", "E. arctica, E028", 27)
# later # E028 <- temer(E028) # did interactive fit and overwrote object
use_data(E030, overwrite = T)
use_data(E028, overwrite = T)
