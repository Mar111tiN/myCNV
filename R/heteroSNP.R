library(tidyverse)
library(glue)

outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/heteroSNP/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)



read_tsv(glue(outpath,"01_A.snp.csv")) %>% 
  filter(Chr == "chr7") %>% 
  filter(Depth > 50) %>% 
  filter(VAF > 0.04 & VAF < 0.95) %>% 
  ggplot(aes(ExonPos, VAF)) +
  geom_point()

read_tsv(glue(outpath,"03_A.snp.csv")) %>% 
  filter(VAF > 0.05 & VAF < 0.98) %>%
  ggplot(aes(VAF)) +
  geom_histogram(bins=200)


read_tsv(glue(outpath,"03_A.snp.csv"))
