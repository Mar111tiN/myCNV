library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)


color.pal <- pal_uchicago("dark")(8)
chrom.cols <- colorRampPalette(color.pal, 23)

outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)

read_tsv(glue(outpath,"heteroSNP/01_A.snp"))

read_tsv(glue(outpath,"heteroSNP/01_A.snp")) %>% 
  filter(Depth > 20 & EBscore > 0.5) %>% 
  ggplot() +
  geom_point(
    aes(
      x=FullExonPos,
      y=VAF, 
      size=EBscore,
      color=Chr
      )) +
      scale_color_manual(values = chrom.cols) +
      theme_bw(base_size = 20)
  )
