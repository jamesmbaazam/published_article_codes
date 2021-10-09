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

input_ICs <- args[1]
output_path <- args[2]


# plot all simulated ICs, and those chosen

inits = readRDS(input_ICs)
max_time = dim(inits[["full"]])[3]


inits_new_cases_full = inits[["full"]][,"new_cases",]
inits_new_cases_full= melt(inits_new_cases_full)
colnames(inits_new_cases_full) = c("sims","time","new_cases")
inits_new_cases_full$select = ifelse(inits_new_cases_full$sims == inits[["mid"]], "Mid", "none")
inits_new_cases_full$select_flag = ifelse(inits_new_cases_full$select %in% c("Low", "Mid","High"),TRUE, FALSE)

text = inits_new_cases_full %>% 
  group_by(sims, time, select, select_flag) %>%
  filter(time == max_time, select_flag == TRUE)
text$lab = paste0(round(text$new_cases/100,2),"%")


p1a = ggplot(data = inits_new_cases_full, aes(x = time, y = new_cases, group = sims, color = select, size = select_flag, alpha = select_flag))+
  geom_vline(aes(xintercept = 15), linetype = "dashed")+
  geom_line(data = inits_new_cases_full %>% filter(select_flag == FALSE))+
  geom_line(data = inits_new_cases_full %>% filter(select_flag == TRUE))+
  geom_text(data = text, aes(x = time + 1, label = lab), size = 3.75, hjust = 0)+
  guides(alpha = FALSE, size = FALSE)+
  scale_alpha_manual(values =c(0.2,1))+
  scale_color_manual(values = c("black","grey"))+
  scale_size_manual(values = c(0.7,1))+
  scale_x_continuous(expand = c(0,0), name = "Days since initial seeding", limits = c(0,max_time + 17))+
  scale_y_continuous(expand = c(0,0), name = "Daily incidence")+
  theme_classic()+
  theme(legend.position = "none")


inits_prev_full = inits[["full"]][,"I1",] + inits[["full"]][,"I2",] + inits[["full"]][,"I3",] + inits[["full"]][,"I4",] + inits[["full"]][,"A",]
inits_prev_full= melt(inits_prev_full)
colnames(inits_prev_full) = c("sims","time","prevalence")
inits_prev_full$select = ifelse(inits_prev_full$sims == inits[["mid"]], "Mid", "none")
inits_prev_full$select_flag = ifelse(inits_prev_full$select %in% c("Low", "Mid","High"),TRUE, FALSE)

text = inits_prev_full %>% 
  group_by(sims, time, select, select_flag) %>%
  filter(time == max_time, select_flag == TRUE)
text$lab = paste0(round(text$prevalence/100,2),"%")

p1b = ggplot(data = inits_prev_full, aes(x = time, y = prevalence, group = sims, color = select, size = select_flag, alpha = select_flag))+
  geom_vline(aes(xintercept = 15), linetype = "dashed")+
  geom_line(data = inits_prev_full %>% filter(select_flag == FALSE))+
  geom_line(data = inits_prev_full %>% filter(select_flag == TRUE))+
  geom_text(data = text, aes(x = time + 1, label = lab), size = 3.75, hjust = 0)+
  guides(alpha = FALSE, size = FALSE, color = FALSE)+
  scale_alpha_manual(values =c(0.2,1))+
  scale_color_manual(values = c("black","grey"))+
  scale_size_manual(values = c(0.7,1))+
  scale_x_continuous(expand = c(0,0), name = "Days since initial seeding", limits = c(0,max_time + 17))+
  scale_y_continuous(expand = c(0,0), name = "Daily prevalence")+
  theme_classic()+
  theme(legend.position = "none")


# for recovered individuals
inits_recov_full = inits[["full"]][,"R",]
inits_recov_full= melt(inits_recov_full)
colnames(inits_recov_full) = c("sims","time","recoveries")
inits_recov_full$select = ifelse(inits_recov_full$sims == inits[["mid"]], "Mid", "none")
inits_recov_full$select_flag = ifelse(inits_recov_full$select %in% c("Low", "Mid","High"),TRUE, FALSE)

text = inits_recov_full %>% 
  group_by(sims, time, select, select_flag) %>%
  filter(time == max_time, select_flag == TRUE)
text$lab = paste0(round(text$recoveries/100,1),"%")

p2 = ggplot(data = inits_recov_full, aes(x = time, y = recoveries, group = sims, color = select, size = select_flag, alpha = select_flag))+
  geom_vline(aes(xintercept = 15), linetype = "dashed")+
  geom_line(data = inits_recov_full %>% filter(select_flag == FALSE))+
  geom_line(data = inits_recov_full %>% filter(select_flag == TRUE))+
  geom_text(data = text, aes(x = time + 1, label = lab), size = 3.75, hjust = 0)+
  scale_alpha_manual(values =c(0.2,1))+
  scale_color_manual(values = c("black","grey"))+
  scale_size_manual(values = c(0.7,1))+
  scale_x_continuous(expand = c(0,0), name = "Days since initial seeding", limits = c(0,max_time + 17))+
  scale_y_continuous(expand = c(0,0), name = "Cumulative recoveries")+
  theme_classic()+
  theme(legend.position = "none")


p = plot_grid(p1a+theme(legend.position = "none"), 
              p1b+theme(legend.position = "none"), 
              p2+theme(legend.position = "none"), 
              align = "vh",labels = c("A","B", "C"),nrow =1)
ggsave(file.path(output_path),p,width = 14, height = 3.5, units = "in")


