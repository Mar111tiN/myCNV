library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)

# set the colors
color.pal <- pal_uchicago("dark")(8)
# expand color palette
chrom.cols <- colorRampPalette(color.pal)(23)


outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)

glimpse(read_tsv(glue(outpath,"heteroSNP/01_A.snp")))

plot.rolling.SNP <- function(path, col, chrom) {
  plot <- read_tsv(glue(outpath, path)) %>% 
    filter(Chr == chrom) %>% 
    ggplot(aes(
      x=FullExonPos,
      y=VAF, 
      color=Chr
    )) +
    geom_point(
      size=.25
    ) +
    geom_line(
      aes(FullExonPos, !!sym(glue(col, "L"))),
      color="red",
      alpha=.25
    ) +
    geom_line(
      aes(FullExonPos, !!sym(glue(col, "R"))),
      color="blue",
      alpha=.25
    ) + 
    geom_line(
      aes(FullExonPos, !!sym(glue(col, "Diff"))),
      color="green"
    ) +
    geom_line(
      aes(FullExonPos, !!sym(col)),
      color="orange",
      size=1.5
    ) + guides(
      color=FALSE
    ) +
    scale_color_manual(values = chrom.cols) +
    theme_bw(base_size = 20)
  
  return(plot)
}
  

'deltaVAFstd'
'absVAFsum'
'VAFstd'

(plot <- plot.rolling.SNP("heteroSNP/01_A.snp",'absVAFsum', "chr7"))


###### rolling coverage info

glimpse(read_tsv(glue(outpath,"covDif/01_A.cov")))

plot.rolling.cov <- function(sample, col, chrom) {
  snp.data <- read_tsv(glue(outpath, "heteroSNP/", sample, ".snp")) %>% 
    filter(Chr == chrom)
  cov.data <- read_tsv(glue(outpath, "covDif/", sample, ".cov")) %>% 
    filter(Chr == chrom) %>% 
    filter(Coverage > 100)
  
  plot <- cov.data %>% 
    ggplot(aes(
      x=FullExonPos,
      y=log2ratio
    )) +
    geom_point(color="green") + 
    geom_line(
      aes(FullExonPos, !!sym(glue(col, "Diff"))),
      color="green"
    ) +
    geom_line(
      aes(FullExonPos, !!sym(col)),
      color="orange",
      size=1.5
    ) + guides(
      color=FALSE
    ) +
    geom_point(
      data=snp.data,
      aes(
        x=FullExonPos,
        y=VAF
      )
    ) +
    scale_color_manual(values = chrom.cols) +
    theme_bw(base_size = 20)
  return(plot)
}

(plot <- plot.rolling.cov("01_A",'log2ratiomean', "chr7"))







# look at all the snps
read_tsv(glue(outpath,"heteroSNP/01_A.snp")) %>% 
  ggplot(aes(offCenter, log2ratio, color=Chr)) +
  geom_point(
    size=0.2,
    alpha=0.4
  ) + 
  scale_y_continuous(limits=c(-2,2.5))





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


