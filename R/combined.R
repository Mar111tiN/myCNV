library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)


color.pal <- pal_uchicago("dark")(8)
color.pal

chrom.cols <- colorRampPalette(color.pal)(23)

chrom.cols

outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)

read_tsv(glue(outpath,"heteroSNP/01_A.snp"))

read_tsv(glue(outpath,"heteroSNP/01_A.snp")) %>% 
  filter(Chr=="chr17") %>% 
  ggplot(aes(
    x=FullExonPos,
    y=VAF, 
    alpha=1 - EBscore,
    color=Chr
  )) +
  geom_point(
    size=.5
  ) +
  geom_line(
    aes(FullExonPos, scoreL),
    color="red"
  ) +
  geom_line(
    aes(FullExonPos, scoreR),
    color="blue"
  ) +
  scale_color_manual(values = chrom.cols) +
  theme_bw(base_size = 20) +
  geom_line(
    aes(FullExonPos, scoreDiff),
    color="green"
  )


# look at all the outliers
read_tsv(glue(outpath,"heteroSNP/01_A.snp")) %>% 
  # filter(scoreL > 0.7 | scoreR > 0.7) %>% 
  filter(Chr == "chr1") %>% 
  ggplot(aes(FullExonPos, VAF, color=Chr)) +
  geom_point()





geom_line(
  aes(FullExonPos, stdL),
  color="red"
) +
  geom_line(
    aes(FullExonPos, stdR),
    color="blue"
  ) +


