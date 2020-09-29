library(tidyverse)
library(glue)

outpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/output/"
plotpath <- "/Users/martinscience/Dropbox/Icke/Work/somVar/tooldata/myCNVdata/plot/"
list.files(outpath)
read_tsv(glue(outpath, "PON_coverage_tidy.csv"))

### show the coverage distribution

plot.cov <- function (coverage.file, sample.list) {
  return(read_tsv(glue(outpath, coverage.file)) %>% 
           filter(sample %in% sample.list) %>% 
           ggplot(aes(ExonPos, Coverage, color=sample)) +
           geom_line() +
           scale_x_continuous(limits=c(0,2423534)) +
           scale_y_continuous(limits=c(0,2000))
  )
}

# full coverage on three samples
(pon.cov <- plot.cov("PON_coverage_tidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B")))
ggsave(glue(plotpath, "PONcoverage.png"), plot=pon.cov, width=12, height=6)
       
# normalized coverage on three samples       
(norm.cov <- plot.cov("PON_coverage_normtidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B")))
ggsave(glue(plotpath, "PONnormalizedCoverage.png"), plot=norm.cov, width=12, height=6)


## zoom in
zoom.cov <- function (coverage.file, sample.list) {
  return(read_tsv(glue(outpath, coverage.file)) %>% 
           filter(sample %in% sample.list) %>% 
           ggplot(aes(ExonPos, Coverage, color=sample)) +
           geom_line() +
           scale_x_continuous(limits=c(500000,505000)) +
           scale_y_continuous(limits=c(0,800))
  )
}

(zoom.cov.plot <- zoom.cov("PON_coverage_tidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B", "013_B", "015_B", "008_B")))
ggsave(glue(plotpath, "ZOOMcoverage.png"), plot=zoom.cov.plot, width=12, height=6)

(zoom.norm.cov.plot <- zoom.cov("PON_coverage_normtidy.csv", sample.list=c("001_B", "002_B", "003_B", "004_B", "013_B", "015_B", "008_B")))
ggsave(glue(plotpath, "ZOOMnormalizedCoverage.png"), plot=zoom.norm.cov.plot, width=12, height=6)


# PON coverage as a mean of all samples

(Pon.cov.plot <- read_tsv(glue(outpath, "PON_coverage.chr7.csv")) %>% 
  ggplot(aes(ExonPos, meanCov)) +
  geom_line() +
  scale_x_continuous(limits=c(500000,505000)) +
  scale_y_continuous(limits=c(0,800)) + 
  geom_ribbon(aes(ymin=meanCov - std, ymax = meanCov + std), alpha=0.6))
ggsave(glue(plotpath, "ZOOMPONcoverage.png"), plot=Pon.cov.plot, width=12, height=6)



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
