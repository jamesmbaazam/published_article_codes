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

input_contour_sims_delay <- args[1]
input_contour_sims_immed <- args[2]
input_IC <- args[3]
input_viz_tools <- args[4]
input_fig_to_make <- args[5] # if "both" include both models on asymp fig, if "immed" only include immediate, if "diffiso" plot differences
input_example_flag <- args[6] # if "TRUE", highlight example from text, if "FALSE" do not include
output_path <- args[7]

source(input_viz_tools)
IC = read.csv(input_IC)

# read in both
out_i = read.csv(input_contour_sims_immed)
out_i$symp_delay = FALSE
out_d = read.csv(input_contour_sims_delay)
out_d$symp_delay = TRUE
out = rbind(out_d, out_i)
out[which(round(out$Ds,8) == round(1/24,8)),"Ds"] = 1/24 # adjust for truncation when saving

## plot number of isolations
plot_iso = out %>% filter(objective %in% c("sympiso", "asympiso", "tot.cases"), contour_FLAG == 1)
plot_iso = cast(plot_iso, tot.tests + d + Ds + symp_delay~ objective, value = "median")
plot_iso$totiso = with(plot_iso, asympiso+sympiso)
plot_iso$pctasymp = with(plot_iso, asympiso/totiso)
plot_iso$d_new = ifelse(plot_iso$symp_delay == TRUE, -1*plot_iso$d, plot_iso$d)

plot_iso$denom = plot_iso$tot.cases + with(IC, I1+I2+I3+I4+A+WI1+WI2+WI3+WI4+WA)
plot_iso$tot_iso_frac = with(plot_iso, (asympiso + sympiso)/denom)
plot_iso$asymp_iso_frac = with(plot_iso, (asympiso)/denom)

line = data.frame(x = seq(0,max(plot_iso$tot_iso_frac),0.01))

plot_iso_sub = rbind(plot_iso %>% filter(round(d,6) %in% round(seq(0.1,0.3,by = 0.1),6)))

if(input_fig_to_make == "both"){
  # both models together
  p = ggplot(data = plot_iso)+
    geom_line(data = line, aes(x = x, y = x), color = "grey")+
    geom_line(data = line, aes(x = x, y = x/2), color = "grey")+
    geom_line(data = line, aes(x = x, y = x/4), color = "grey")+
    geom_line(data = line, aes(x = x, y = 3*x/4), color = "grey")+
    geom_line(aes(x = tot_iso_frac, y = asymp_iso_frac, group = paste(tot.tests, symp_delay), alpha = symp_delay))+
    geom_path(data = plot_iso_sub, aes(x = tot_iso_frac, y = asymp_iso_frac, color = d_new, group = paste(d,symp_delay), alpha = symp_delay), linetype = "dashed")+
    geom_point(aes(x = tot_iso_frac, y = asymp_iso_frac, size = Ds,  fill = symp_delay,  color = d_new, alpha = symp_delay))+
    scale_alpha_manual(values = c(0.8,0.8))+
    scale_color_gradientn(colours = c("grey20", "grey40", "grey60", "grey80","grey95", brewer.pal(4, "Reds")), limits = c(-max(abs(plot_iso$d_new)),max(abs(plot_iso$d_new))))+
    scale_size_continuous(breaks = c(1/24,1,2,5), labels = c("Rapid","1 day", "2 days", "5 days"), name = c("Delay"))+
    scale_x_continuous(expand = c(0,0), name = "Total isolations\n(as % of all infections)", limits = c(0, (max(line)+0.025)), labels = percent_format(accuracy = 1))+
    scale_y_continuous(expand = c(0,0), name = "Isolations before reporting\n(as % of all infections)", labels = percent_format(accuracy = 1))+
    theme_classic()+
    theme(legend.position = "none")
  
  # move legends
  l_grey = get_legend(ggplot(data = plot_iso %>% filter(symp_delay== TRUE))+
                        geom_point(aes(x = sympiso + asympiso, y = asympiso,  color = -d_new))+
                        scale_color_gradientn(colours = c("grey95", "grey80", "grey60", "grey40","grey20"), 
                                              limits = c(0,max(abs(plot_iso$d_new))), 
                                              name = "", labels = percent_format(accuracy = 1))+
                        theme(legend.text = element_blank()))
  l_red = get_legend(ggplot(data = plot_iso %>% filter(symp_delay== FALSE))+
                       geom_point(aes(x = sympiso + asympiso, y = asympiso,  color = d_new))+
                       scale_color_gradientn(colours = c("grey95", brewer.pal(4, "Reds")), 
                                             limits = c(0,max(abs(plot_iso$d_new))), 
                                             name = "", labels = percent_format(accuracy = 1)))
  l_size = get_legend(ggplot(data = plot_iso %>% filter(symp_delay== TRUE))+
                        geom_point(aes(x = sympiso + asympiso, y = asympiso, size = Ds))+
                        scale_size_continuous(breaks = c(1/24,1,2,5), labels = c("Rapid","1 day", "2 days", "5 days"), name = c("Delay"))+
                        theme_classic())
  p = p +
    theme(legend.position = "none")+
    annotation_custom(grob=l_grey,ymin = max(line$x)/1.5, ymax=Inf, xmin=0.025, xmax=0.09) +
    annotation_custom(grob=l_red,ymin = max(line$x)/1.5, ymax=Inf, xmin=0.09, xmax=0.155) +
    annotation_custom(grob=l_size,ymin = max(line$x)/4, ymax=max(line$x)/1.5, xmin=0.025, xmax=0.155) +
    geom_text(x = 0.025, y = 0.69, label = "NPI Intensity", hjust = 0, size = 3.5)
  ggsave(file.path(output_path),p, width = 6, height = 5)
  
}

if(input_fig_to_make == "immed"){
  # immediate only
  text = plot_iso_sub %>% filter(symp_delay == FALSE, tot.tests == 5000) 
  text$x = with(text,tot_iso_frac+0.025)
  text$y = with(text,asymp_iso_frac-0.025)
  text$lab=paste0(text$d*100,"%")
  
  text2 = data.frame(x = rep(max(line$x),4)+max(line$x)/100,
                     y = c(max(line$x),max(line$x)/2,max(line$x)/4,3*max(line$x)/4),
                     lab = c("100% of\nisolations\nbefore\nreporting","50%","25%","75%"))
  
  p = ggplot(data = plot_iso %>% filter(symp_delay == FALSE))+
    geom_line(data = line, aes(x = x, y = x), color = "grey")+
    geom_line(data = line, aes(x = x, y = x/2), color = "grey")+
    geom_line(data = line, aes(x = x, y = x/4), color = "grey")+
    geom_line(data = line, aes(x = x, y = 3*x/4), color = "grey")+
    geom_line(aes(x = tot_iso_frac, y = asymp_iso_frac, group = tot.tests))+
    geom_text(data = text2, aes(x = x, y = y, label = lab), color = "grey", size = 2, hjust = 0, vjust = 1)+
    geom_path(data = plot_iso_sub %>% filter(symp_delay == FALSE), aes(x = tot_iso_frac, y = asymp_iso_frac, color = d, group = d), linetype = "dashed")+
    geom_point(aes(x = tot_iso_frac, y = asymp_iso_frac, color = d, size = Ds), alpha = 0.8)+
    geom_label(data = text, aes(x = x,y=y,label = lab,color = d), fill = "white", label.size = 0, size = 2)+
    scale_color_distiller(palette = "Reds", direction = 0, name = "NPI Intensity", labels = percent_format(accuracy = 1))+
    scale_size_continuous(breaks = c(1/24,1,2,5), labels = c("Rapid","1 day", "2 days", "5 days"), name = c("Delay"))+
    scale_x_continuous(expand = c(0,0), name = "Total isolations\n(as % of all infections)", limits = c(0, max(line)+0.065), labels = percent_format(accuracy = 1))+
    scale_y_continuous(expand = c(0,0), name = "Isolations before reporting\n(as % of all infections)", labels = percent_format(accuracy = 1))+
    theme_classic()
  
  l_red = get_legend(p)
  
  p = p +
    theme(legend.position = "none")+
    annotation_custom(grob=l_red,ymin = max(line$x)/2.5, ymax=max(line$x)-0.020, xmin=0.05, xmax=0.2)
  
  if(input_example_flag == "TRUE"){
    exmpl = plot_iso %>% filter(symp_delay == FALSE, d == 0.35, tot.tests == 500)
    exmpl = rbind(exmpl,
                  plot_iso %>% filter(symp_delay == FALSE, d == 0.25, tot.tests == 2000))
    
    p = p + 
      geom_point(data = exmpl, aes(x = tot_iso_frac, y = asymp_iso_frac, shape = as.factor(asymp_iso_frac))) +
      scale_shape_manual(values = c(4,8)) + 
      guides(shape = FALSE)
  }
  
  ggsave(file.path(output_path), p, width = 5, height = 5)
  
}

# values for text
# plot_iso %>% filter(d == 0.35, tot.tests == 500, symp_delay == FALSE)
# plot_iso %>% filter(d == 0.35, tot.tests == 5000, symp_delay == FALSE)


#### difference in isolations between two models ####
if(input_fig_to_make == "diffiso"){
  plot_iso_tot_iso = cast(plot_iso%>% filter(Ds %in% c(1/24,seq(0.5,8,by = 0.5))), tot.tests + Ds~ symp_delay, value = "tot_iso_frac")
  colnames(plot_iso_tot_iso)[3:4] = c("nodelay", "delay")
  plot_iso_tot_iso$diff = plot_iso_tot_iso$delay - plot_iso_tot_iso$nodelay
  
  plot_iso_totiso_comb = merge(plot_iso_tot_iso, plot_iso[plot_iso$symp_delay == FALSE,c("tot.tests", "Ds", "d")])
  colnames(plot_iso_totiso_comb)[6] = "d_nodelay"
  plot_iso_totiso_comb = merge(plot_iso_totiso_comb, plot_iso[plot_iso$symp_delay == TRUE,c("tot.tests", "Ds", "d")])
  colnames(plot_iso_totiso_comb)[7] = "d_delay"
  plot_iso_totiso_comb$d_diff = with(plot_iso_totiso_comb, d_delay - d_nodelay)
  
  # values for text
  # plot_iso_totiso_comb %>% filter(d_diff == max(plot_iso_totiso_comb$d_diff))
  
  
  r = ggplot(data = plot_iso_totiso_comb)+
    geom_line(aes(x = diff, y = d_diff, color = as.factor(tot.tests)), size = 1)+
    geom_point(aes(x = diff, y = d_diff,  fill = as.factor(tot.tests), size = Ds), shape = 21)+
    guides(color = FALSE)+
    scale_color_brewer(palette = "YlGnBu", name = "Tests/day", labels = c("1%", "5%","20%","50%"))+
    scale_fill_brewer(palette = "YlGnBu", name = "Tests/day", labels = c("1%", "5%","20%","50%"))+
    scale_size_continuous(range = c(0.9,4), name = "Delay", breaks = c(1/24, 1, 2,5), labels = c("Rapid", "1 day", "2 days", "5 days"))+
    scale_x_continuous(name = "Difference in total isolations\n(as % of all infections)", labels = percent_format(accuracy = 1))+
    scale_y_continuous(labels = percent_format(accuracy = 1), name = "Difference in NPI intensity")+
    theme_bw()+
    theme(legend.direction = "vertical", 
          legend.position = "bottom")
  legend = get_legend(r)
  r = r + annotation_custom(grob=legend,xmin=min(plot_iso_totiso_comb$diff), 
                            xmax=min(plot_iso_totiso_comb$diff)+0.15, ymax = 0,  ymin = 0.05)+
    #annotation_custom(grob = l_size, ymin = -0.05, ymax =0, xmin = 250, xmax = 350)+
    theme(legend.position = "none")
  ggsave(file.path(output_path), r, width = 5, height = 5)
}
