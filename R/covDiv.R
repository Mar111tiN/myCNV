library(tidyverse)
library(glue)

## SET PATHS
outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
data.type <- "PONcoverage"

## APPLY PATHS
outpath <- glue(outpath, data.type, "/")
plotpath <- glue(plotpath, data.type, "/")
list.files(outpath)

# INSPECT DATA
data <- read_tsv(glue(outpath, "PON_coverage_mean.csv"))

#
### show the coverage distribution
plot.cov <- function (coverage.file, sample.list="all", zoom=c(0,2423534), ymax=500) {
  
  # load the data and tidy
  data <- read_tsv(glue(outpath, coverage.file)) %>% 
    select(c(1,2,3) | ends_with("_B")) %>% 
    gather(ends_with("_B"), key=sample, value=Coverage)
  
  if (sample.list != "all") {
    data <- data %>% filter(sample %in% sample.list)
  }
  if ((zoom[2] - zoom[1]) > 10000) {
    plot <- data %>% 
      ggplot(aes(ExonPos, Coverage, color=sample)) +
      geom_point(size=.1, alpha=.1)
  } else {
    plot <- data %>% 
      ggplot(aes(ExonPos, Coverage, color=sample)) +
      geom_line()
  }
  
  # further edit the plot (scales etc)
  plot <- plot + 
    scale_x_continuous(limits=zoom) +
    scale_y_continuous(limits=c(0,ymax))
  return(plot)
}

# full coverage on three samples
(pon.cov <- plot.cov("PON_coverage.csv", sample.list=c("001_B", "002_B", "003_B", "004_B"), ymax=1000))
 ggsave(glue(plotpath, "PONcoverage.png"), plot=pon.cov, width=12, height=6)
       
# normalized coverage on three samples       
(norm.cov <- plot.cov("PON_coverage_normtidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B")))
 (norm.cov <- plot.cov("PON_coverage_normtidy.csv", zoom=c(500000,505000)))
ggsave(glue(plotpath, "PONnormalizedCoverage.png"), plot=norm.cov, width=12, height=6)


(zoom.cov.plot <- zoom.cov("PON_coverage_tidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B", "013_B", "015_B", "008_B")))
ggsave(glue(plotpath, "ZOOMcoverage.png"), plot=zoom.cov.plot, width=12, height=6)

(zoom.norm.cov.plot <- zoom.cov("PON_coverage_normtidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B", "013_B", "015_B", "008_B")))
(zoom.norm.cov.plot <- zoom.cov("PON_coverage_normtidy.csv"))
ggsave(glue(plotpath, "ZOOMnormalizedCoverage.png"), plot=zoom.norm.cov.plot, width=12, height=6)


################ RUN THE MEANS
# PON coverage as a mean of all samples

# INSPECT DATA
data <- read_tsv(glue(outpath, "PON_coverage_mean.csv"))

data %>% 
  select(c(1,2,3, "meanCov", "medianCov", "std"))

# now the combined plot function
plot.mean.cov <- function (coverage.file, sample.list="all", zoom=c(0,2423534), ymax=500) {
  
  # load the data
  all.data <- read_tsv(glue(outpath, coverage.file))
  
  # split data into avg and sample data
  data.samples <- all.data  %>% 
    select(c(1,2,3) | ends_with("_B")) %>% 
    gather(ends_with("_B"), key=sample, value=Coverage)
  
  data.mean <- all.data %>% 
    select(c(1,2,3, "meanCov", "medianCov", "std"))
  
  
  if (sample.list != "all") {
    data.samples <- data.samples %>% filter(sample %in% sample.list)
  }
  if ((zoom[2] - zoom[1]) > 10000) {
    plot <- data.samples %>% 
      ggplot() +
      geom_point(
        aes(ExonPos, Coverage, color=sample),
        size=.1,
        alpha=.1
        )
  } else {
    plot <- data.samples %>% 
      ggplot() +
      geom_line(
        aes(ExonPos, Coverage, color=sample),
        alpha=.4
      )
  }
  
  # further edit the plot (scales etc)
  plot <- plot + 
    scale_x_continuous(limits=zoom) +
    scale_y_continuous(limits=c(0,ymax))
  
  # add the mean values
  plot <- plot + 
    geom_line(
      data=data.mean,
      aes(ExonPos, meanCov),
      size=0.2,
      alpha=.4) +
    geom_ribbon(
      data=data.mean,
      aes(
        x=ExonPos,
        y=meanCov,
        ymin=meanCov - std,
        ymax = meanCov + std
        ),
      alpha=0.6
      )
  return(plot)
}


(plot <- plot.mean.cov("PON_coverage_mean.csv"))

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
