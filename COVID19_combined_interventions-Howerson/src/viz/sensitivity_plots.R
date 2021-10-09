### Preamble
suppressWarnings({
  suppressMessages({
    require(ggplot2)
    require(tidyr)
    require(gridExtra)
    require(gridExtra)
    require(RColorBrewer)
    require(plyr)
    require(dplyr)
    require(cowplot)
    require(reshape)
    require(reshape2)
    require(grid)
    require(scales)
    })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_model_sims <- args[1]
input_sens_folder <- args[2] # should contain sensitivity_IC_high_simulations.csv, 
                                           #  sensitivity_IC_higher_simulations.csv
                                           #  sensitivity_asymp_simulations.csv
                                           #  sensitivity_asymp_rho_simulations.csv
input_viz_tools <- args[3]
output_path <- args[4]

source(input_viz_tools)

# visualize sensitivity analyses

# read in files
out_tmp = read.csv(paste0(input_sens_folder,"/sensitivity_IC_high_simulations.csv"))
out_tmp$sens_type = "high"
out_tmp$sens_type2 = "Prior immunity: 10%"
out = out_tmp

out_tmp = read.csv(paste0(input_sens_folder,"/sensitivity_IC_higher_simulations.csv"))
out_tmp$sens_type = "higher"
out_tmp$sens_type2 = "Prior immunity: 20%"
out = rbind(out, out_tmp)

out_tmp = read.csv(paste0(input_sens_folder,"/sensitivity_asymp_simulations.csv"))
out_tmp$sens_type = paste("p",out_tmp$p, sep = "_")
out_tmp$sens_type2 = paste0("Asymptomatic: ",out_tmp$p*100, "%")
out = rbind(out, out_tmp[,colnames(out_tmp) != "p"])

out_tmp = read.csv(paste0(input_sens_folder,"/sensitivity_asymp_rho_simulations.csv"))
out_tmp$sens_type = paste("rho", out_tmp$rho, sep="_")
out_tmp$sens_type2 = paste("Relative infectiousness:",out_tmp$rho)
out = rbind(out, out_tmp[,colnames(out_tmp) != "rho"])

out_tmp = read.csv(input_model_sims)
out_tmp = out_tmp %>% filter(sens == 1)
out_tmp$sens_type = "original"
out_tmp$sens_type2 = "Base"
out = rbind(out, out_tmp)



p1 = ggplot(data = out %>% filter(objective == "tot.cases", d %in% seq(0,0.3,0.1), sens_type %in% c("original", "high", "higher")))+
  geom_contour(aes(x = Ds, y = tot.tests/10000, z = median, color = as.factor(sens_type2)), breaks = c(250))+
  facet_grid(cols = vars(d), labeller = labeller(d = labs_d))+
  scale_color_manual(values = c("black",brewer.pal(5,"Dark2")[1:2]), name = "")+
  scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)", limits = c(0,8))+
  scale_y_continuous(expand = c(0,0), name = "Test administration\n(% population tested daily)", label = percent_format(accuracy = 1), limits = c(0,0.5))+
  theme_bw()+
  theme(legend.position = "bottom")

  
p2 = ggplot(data = out %>% filter(objective == "tot.cases", d %in% seq(0,0.3,0.1), sens_type %in% c("original", "p_0.2", "p_0.6")))+
  geom_contour(aes(x = Ds, y = tot.tests/10000, z = median, color = as.factor(sens_type2)), breaks = c(250))+
  facet_grid(cols = vars(d), labeller = labeller(d = labs_d))+
  scale_color_manual(values = c(brewer.pal(5,"Dark2")[3:4],"black"), name = "")+
  scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)", limits = c(0,8))+
  scale_y_continuous(expand = c(0,0), name = "Test administration\n(% population tested daily)", label = percent_format(accuracy = 1), limits = c(0,0.5))+
  theme_bw()+
  theme(legend.position = "bottom")


p3 = ggplot(data = out %>% filter(objective == "tot.cases", d %in% seq(0,0.3,0.1), sens_type %in% c("original", "rho_0.5")))+
  geom_contour(aes(x = Ds, y = tot.tests/10000, z = median, color = as.factor(sens_type2)), breaks = c(250))+
  facet_grid(cols = vars(d), labeller = labeller(d = labs_d))+
  scale_color_manual(values = c("black",brewer.pal(5,"Dark2")[5]), name = "")+
  scale_x_continuous(expand = c(0,0), name = "Average test delay\n(days from test to isolation)", limits = c(0,8))+
  scale_y_continuous(expand = c(0,0), name = "Test administration\n(% population tested daily)", label = percent_format(accuracy = 1), limits = c(0,0.5))+
  theme_bw()+
  theme(legend.position = "bottom")


p =plot_grid(p1,p2,p3,ncol = 1, align = "h", labels = c('A','B','C'))

ggsave(file = file.path(output_path),p, width = 9, height =9, units = "in")


