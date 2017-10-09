# the first simple example
library(RSA)
fold.change <- read.csv(system.file("extdata","foldchange.csv", package = "RSA"))
head(RSA(fold.change))

# the secondary example
