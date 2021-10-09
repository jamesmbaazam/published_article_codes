### Preamble
suppressWarnings({
  suppressMessages({
    require(tidyr)
    require(gridExtra)
    require(ggplot2)
    require(gridExtra)
    require(RColorBrewer)
    require(plyr)
    require(dplyr)
    require(cowplot)
    require(reshape)
    require(reshape2)
    require(metR)
    require(grid)
    require(scales)
    #require(viridis)
  })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_sims <- args[1]
output_path <- args[2]


# plot correlations between outcomes
out = read.csv(input_sims)

out_plot = cast(out, d+tot.tests + sens+Ds+sim ~ objective, value = "median")

ggplot(data = out_plot %>% filter(sens == 1))+
  geom_point(aes(x = tot.cases, y = tot.hosp), size = 1.5, alpha = 0.8)+
  geom_point(aes(x = tot.cases, y= death), size = 1.5, alpha = 0.8, color = "blue")+
  scale_x_continuous(expand = c(0,0), name = "Total infections")+
  scale_y_continuous(expand = c(0,0), name = "Outcome")+
  theme_classic()
ggsave(file.path(output_path),width = 5, height = 5, units = "in")

