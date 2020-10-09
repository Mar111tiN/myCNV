library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)

# set the colors
color.pal <- pal_uchicago("dark")(8)
# expand color palette
chrom.cols <- colorRampPalette(color.pal)(23)


outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/cluster/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)

glimpse(read_tsv(glue(outpath,"01_A.cnv")))

read_tsv(glue(outpath,"01_A.cnv")) %>% 
  drop_a

# plot function for all important data points
plot.rolling.cov <- function(sample, chrom, cols=list('log2ratiomean'='blue'), plot.diff=FALSE) {
  data <- read_tsv(glue(outpath, sample, ".cnv"))
  # draw first plot:-)
  plot <- data %>% 
    ggplot() +
    geom_point(
      aes(FullExonPos, VAF),
      color="black",
      alpha=.8,
      size=.8) +
    scale_color_manual(values = chrom.cols) +
    theme_bw(base_size = 20)
  
  for (col in names(cols)) {
    here.data <- data %>% drop_na(!!sym(col))
    plot <- plot +
      geom_line(
        data =  here.data,
        aes(FullExonPos, !!sym(col)),
        color=toString(cols[col]),
        size=1.5
      )
    if (plot.diff) {
      plot <- plot +
      geom_line(
        data= here.data,
        aes(FullExonPos, !!sym(glue(col,"Diff"))),
        color=toString(cols[col]),
        size=0.5,
        alpha=0.5
      )
    }
  }
  return(plot)
}

column.code <- list(
  'log2ratiomean'='blue',
  'absVAFsum'='green',
  'deltaVAFstd'='brown',
  'VAFstd'='red'
  )

code <- list(
  'log2ratiomean'='blue'
)
  
  
plot.rolling.cov("01_A_staggered", chrom="chrX", cols=code, plot.diff=TRUE)


glimpse(read_tsv(glue(outpath,"01_A_merged.cnv")))

'deltaVAFstd'
'absVAFsum'
'VAFstd'
'log2ratiomean'

read_tsv(glue(outpath,"01_A_merged.cnv")) %>% 
  ggplot(aes(absVAFsum,VAFstd, color=Chr)) +
   geom_point(
     size=0.5,
     alpha=.5
   ) +
  scale_color_manual(values = chrom.cols) +
  theme_bw(base_size = 20)
  
