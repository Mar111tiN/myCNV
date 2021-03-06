library(tidyverse)
library(glue)

## SET PATHS
outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
data.type <- "covDif"

## APPLY PATHS
outpath <- glue(outpath, data.type, "/")
plotpath <- glue(plotpath, data.type, "/")
list.files(outpath)
outpath
# INSPECT DATA
data <- read_tsv(glue(outpath, "01_A.cov"))
glimpse(data)
#
### show the coverage distribution
plot.cov <- function (data, zoom=c(0,Inf), ymax=200) {
    plot <- data %>% 
      ggplot(aes(FullExonPos, Coverage)) +
        geom_line(size=2) + 
        geom_line(
          aes(FullExonPos, PONmeanCovmean),
          size=0.2,
          alpha=.4) +
        geom_ribbon(
          aes(
            x=FullExonPos,
            y=PONmeanCovmean,
            ymin=PONmeanCovmean - PONstd,
            ymax = PONmeanCovmean + PONstd
            ),
          alpha=0.6
        ) + 
      scale_y_continuous(limits=c(0,ymax))
  
  return(plot)
}
data <- read_tsv(glue(outpath, "01_A.cov"))

data %>% 
  filter(Chr == "chr5") %>% 
  plot.cov()




(plot <- plot.cov("01_A.chr7.covDif"))

################ RUN THE FILTER


(plot <- plot.mean.cov("PON_coverage.chr7.removed.csv", zoom=c(2000,28000)))


ggsave(glue(plotpath, "PON_coverage_filter.ZOOM.png"), plot=plot, width=12, height=6)
(plot <- plot.mean.cov("PON_coverage_removed.csv", zoom=c(14000,16000)))
ggsave(glue(plotpath, "PON_coverage_removed.ZOOM.png"), plot=plot, width=12, height=6)

(plot <- plot.mean.cov("PON_coverage_filter.csv"))
ggsave(glue(plotpath, "PON_coverage_filter.png"), plot=plot, width=12, height=6)
(plot <- plot.mean.cov("PON_coverage_removed.csv"))
ggsave(glue(plotpath, "PON_coverage_removed.png"), plot=plot, width=12, height=6)
















#################
##### DIF #######
diff_file <- "003_A.bam.diff.csv"
df <- read_tsv(glue(outpath, diff_file, sep='/'))
glimpse(df)

read_tsv(glue(outpath, "003_A.bam.diff.csv")) %>% 
  filter(normCov > 300) %>% 
  ggplot(aes(ExonPos, rollingPloidy)) +
  geom_line() +
  scale_x_continuous(limits=c(0,2423534)) +
  scale_y_continuous(limits=c(0,3))






# make plot function
plot.diff <- function (coverage.file) {
  return(read_tsv(glue(outpath, coverage.file)) %>% 
           ggplot(aes(ExonPos, ploidy)) +
           geom_line() +
           scale_x_continuous(limits=c(0,2423534)) +
           scale_y_continuous(limits=c(0,3))
  )
}


plot.diff("003_A.bam.diff.csv")



plot.cov <- function (df) {
  return(df %>% 
  ggplot(aes(ExonPos, rploidy)) +
    geom_line() +
    scale_x_continuous(limits=c(0,200)) +
    scale_y_continuous(limits=c(0,2000))
  )
}


read_tsv(glue(outpath, diff_file, sep='/')) %>% 
  plot.cov

read_tsv("bedCov/all_norm_tidy.bedCov") %>% 
  plot.cov
