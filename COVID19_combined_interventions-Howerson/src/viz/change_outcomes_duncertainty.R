### Preamble
suppressWarnings({
  suppressMessages({
    require(tidyr)
    require(gridExtra)
    require(ggplot2)
    require(RColorBrewer)
    require(plyr)
    require(dplyr)
    require(cowplot)
    require(reshape)
    require(reshape2)
    require(metR)
    require(grid)
    require(scales)
  })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_contour_sims <- args[1]
input_viz_tools <- args[2]
output_path <- args[3]


source(input_viz_tools)


# immediate
out = read.csv(input_contour_sims)
out[which(round(out$Ds,8) == round(1/24,8)),"Ds"] = 1/24 # adjust for truncation when saving

out_compare2 = cast(out[,c("tot.tests", "Ds", "objective", "contour_FLAG", "median")], tot.tests + Ds  + objective ~ contour_FLAG, value = "median")
colnames(out_compare2)[colnames(out_compare2) %in% c("0","1")] = c("LessD", "CorrectD")

out_compare2$delta = with(out_compare2, LessD - CorrectD)
out_compare2$pctchange = with(out_compare2, delta/CorrectD)
out_compare2$objective = as.factor(out_compare2$objective)
out_compare2$objective = factor(out_compare2$objective, levels(out_compare2$objective)[c(3,7,6,1:2,4:5)])

# values for text
# out_compare2 %>% filter(objective == "tot.cases") %>% summarise(mean = mean(pctchange, na.rm = TRUE))
# out_compare2 %>% filter(objective == "tot.cases", tot.tests == 2000, Ds %in% c(1, 7))
# out_compare2 %>% filter(objective == "tot.cases", tot.tests %in% c(100,2000), Ds == 2)


p <- ggplot(data = out_compare2 %>% filter(objective == "tot.cases"))+
  geom_point(aes(x = Ds, y = pctchange, fill = as.factor(tot.tests/10000)), size = 3.5, shape = 21, color = "black")+
  scale_fill_brewer(palette = "YlGnBu", name = "Tests/day", labels = c("1%", "5%","20%","50%"))+
  scale_x_continuous(limits = c(0,NA),name = "Average test delay\n(days from test to isolation)")+
  scale_y_continuous(labels = percent_format(accuracy = 1),limits = c(0,NA), name = "Percent increase in infections\n(due to additional lifting of NPIs)")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        plot.margin = margin(10, 10, 10, 10),
        panel.spacing = unit(1.15, "cm"),
        strip.background = element_blank())
legend <- get_legend(p)

p <- p + annotation_custom(legend, xmin = 6, xmax = 8, ymin = 0.005, ymax = 0.075)+
  theme(legend.position = "none")
ggsave(file.path(output_path),p,width = 4, height =4, units = "in")
