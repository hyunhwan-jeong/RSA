# the first simple example
library(RSA)
fold.change <- read.csv(system.file("extdata","foldchange.csv", package = "RSA"))
head(fold.change)
head(RSA(fold.change))

# the secondary example, Hit-Select format UMUC3 dataset
library(tidyverse)
read.count <- read.table(system.file("extdata","readcount.tsv", package = "RSA"), sep="\t", header=T, stringsAsFactors = F)
head(read.count)

ctrl.count <- read.count %>% select(starts_with("CONTROL_")) %>%
  apply(1, median, na.rm = TRUE)
test.count <- read.count %>% select(starts_with("TEST_")) %>%
  apply(1, median, na.rm = TRUE)
FC <- test.count / ctrl.count
df.RSA <- data.frame(Gene_ID=read.count$Gene,
                     Well_ID=read.count$guide.RNA_ID,
                     Score=FC)
head(df.RSA)

ret.RSA <- RSA(df.RSA)
ret.RSA %>% group_by(Gene_ID) %>% summarise(score = mean(LogP)) %>%
  select(gene = Gene_ID, score = score)
