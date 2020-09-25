library(tidyverse)

list.files('bedCov')

plot.cov <- function (df) {
  return(df %>% 
  ggplot(aes(`#`, Coverage, color=sample)) +
    geom_line() +
    scale_x_continuous(limits=c(0,200)) +
    scale_y_continuous(limits=c(0,2000))
  )
}


read_tsv("bedCov/all_tidy.bedCov") %>% 
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
