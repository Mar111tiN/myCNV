library(tidyverse)
library(glue)

outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)
read_tsv(glue(outpath, "PON_coverage_normtidy.csv"))

### show the coverage distribution
plot.cov <- function (coverage.file, sample.list="all", zoom=c(0,2423534), ymax=500) {
  data <- read_tsv(glue(outpath, coverage.file))
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
(pon.cov <- plot.cov("PON_coverage_tidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B"), ymax=1000))
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

# first I need a function extracting the mean, median and std into separate df
get.mean.table <- function(data) {
  mean.cov <- data %>% 
    filter(sample == "meanCov") %>% 
    rename(meanCoverage = Coverage) %>% 
    select(ExonPos, meanCoverage)
  
  median.cov <- data %>% 
    filter(sample == "medianCov") %>% 
    rename(medianCoverage = Coverage) %>% 
    select(ExonPos, medianCoverage)
  
  std.cov <- data %>% 
    filter(sample == "std") %>% 
    rename(std = Coverage) %>% 
    select(ExonPos, std)
  
  
  combined.tibble <- mean.cov %>% 
    full_join(median.cov) %>% 
    full_join(std.cov)
  return(combined.tibble)
}

# test get.mean.table

(data <- read_tsv(glue(outpath, "PON_coverage_mean.csv")) %>% 
    get.mean.table())

# now the combined plot function
plot.mean.cov <- function (coverage.file, sample.list="all", zoom=c(0,2423534), ymax=500) {
  
  # load the data
  all.data <- read_tsv(glue(outpath, coverage.file))
  
  # split data into avg and sample data
  data.samples <- all.data %>% 
    filter(grepl("_B", sample))
  print(data.samples)
  data.mean <- all.data %>% 
    get.mean.table
  
  
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
      aes(ExonPos, meanCoverage),
      size=0.2,
      alpha=.4) +
    geom_ribbon(
      data=data.mean,
      aes(
        x=ExonPos,
        y=meanCoverage,
        ymin=meanCoverage - std,
        ymax = meanCoverage + std
        ),
      alpha=0.6
      )
  return(plot)
}


(plot <- plot.mean.cov("PON_coverage_mean.csv"))

################ RUN THE FILTER


(plot <- plot.mean.cov("PON_coverage_filter.csv", zoom=c(10000,20000)))




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


difCov <- 

(difCov <- read_tsv("bedCov/divCov_tidy.bedCov"))
difCov%>% 
  filter(sample %in% c("016_A","014_A", "021_A")) %>% 
  ggplot(aes(`#`, Coverage, color=sample)) +
  geom_point(size=.5, alpha=.1) +
  scale_y_continuous(limits=c(0,5)) +
  geom_smooth()

unique(difCov$sample)
