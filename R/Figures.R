## Packages ##

library(ggplot2)
library(ggpubr)
library(gridGraphics)
library(simr)
library(rstatix)
library(devtools) 
library(nlcor)
library(lsmeans)
library(ggthemes)
library(RColorBrewer)
library(tidyr)
library(openxlsx)
library(extrafont)
library(cowplot)
library(forcats)
library(dplyr)
library(purrr)
library(stringr)
rm(list=ls())


#######################################################################################################
## Figure 4 & Figure S4 - VP Global effect & Correlation (Global effect ~ Verbal Coordination Index) ##
#######################################################################################################

## Left hemisphere ## 

dataset <- read.csv(paste('add path here','R_dataset_1_left.csv',  sep=''))

# to extract the number of patient by region
region_1 = subset(dataset, atlas_loca=="IFG BA44")
pat_in_region_1 = length(unique(region_1$patients))
region_2 = subset(dataset, atlas_loca=="STG BA22")
pat_in_region_2 = length(unique(region_2$patients))
region_3 = subset(dataset, atlas_loca=="STG BA41/42")
pat_in_region_3 = length(unique(region_3$patients))
region_4 = subset(dataset, atlas_loca=="IPL BA40")
pat_in_region_4 = length(unique(region_4$patients))

# to add this number in the region label
IFG_tick = paste("IFG BA44 \n(n=", pat_in_region_1, ")", sep='')
STG_A22_tick = paste("STG BA22 \n(n=", pat_in_region_2, ")", sep='')
STG_tick = paste("STG BA41/42 \n(n=", pat_in_region_3, ")", sep='')
IPL_tick = paste("IPL BA40 \n(n=", pat_in_region_4, ")", sep='')

# plot parameters
width_box = 0.65
size_tick_xaxis = 12
size_tick_yaxis = 12
size_title = 15
alpha_jitter = 0.9
size_jitter_points = 1
size_legend_text = 10
size_line = 2
roi_list = c("STG BA41/42", "STG BA22", "IPL BA40", "IFG BA44")

fig <- ggplot(transform(dataset, atlas_loca=factor(atlas_loca,levels=roi_list)), aes(x=atlas_loca, y=data)) + 
          geom_hline(yintercept=0, linetype='dashed', size=size_line, color="#F46D43") +
          geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=data_types), outlier.size = size_jitter_points) +
          geom_jitter(aes(group = data_types), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
          ylab("% power changes & rho values") + 
          scale_y_continuous(breaks=c(-1, 0, 1, 2), labels=c('-100', '0', '100', '200'), limits=c(-1, 2),
                             sec.axis = sec_axis(~ 1*., breaks = c(-1, 0, 1),labels=c('-1', '0', '1'))) +
          scale_x_discrete(breaks=c("STG BA41/42", "STG BA22", "IPL BA40", "IFG BA44"),
                           labels=c(STG_tick, STG_A22_tick, IPL_tick, IFG_tick), ) +
          xlab(element_blank()) + ylab(element_blank()) +
          labs(title=paste("HFa (70-125 Hz)", sep='')) + 
          guides(fill=guide_legend(title=element_blank(), byrow = TRUE), size = guide_legend(order = 3)) +
          theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                axis.text.x = element_text(face="plain",size=size_tick_xaxis),
                axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                legend.text = element_text(face="plain",size=size_legend_text),
                legend.position = "bottom", legend.spacing.x = unit(1, 'cm'))
fig
# left y-axis = power changes (%)
# right y-axis = r


## Right hemisphere ##

dataset <- read.csv(paste('add path here','R_dataset_1_right.csv',  sep=''))

# to extract the number of patient by region
region_1 = subset(dataset, atlas_loca=="IFG BA44")
pat_in_region_1 = length(unique(region_1$patients))
region_2 = subset(dataset, atlas_loca=="STG BA22")
pat_in_region_2 = length(unique(region_2$patients))
region_3 = subset(dataset, atlas_loca=="STG BA41/42")
pat_in_region_3 = length(unique(region_3$patients))
region_4 = subset(dataset, atlas_loca=="IPL BA40")
pat_in_region_4 = length(unique(region_4$patients))

# to add this number in the region label
IFG_tick = paste("IFG BA44 \n(n=", pat_in_region_1, ")", sep='')
STG_A22_tick = paste("STG BA22 \n(n=", pat_in_region_2, ")", sep='')
STG_tick = paste("STG BA41/42 \n(n=", pat_in_region_3, ")", sep='')
IPL_tick = paste("IPL BA40 \n(n=", pat_in_region_4, ")", sep='')

# plot parameters
width_box = 0.65
size_tick_xaxis = 12
size_tick_yaxis = 12
size_title = 15
alpha_jitter = 0.9
size_jitter_points = 1
size_legend_text = 10
size_line = 2
roi_list = c("STG BA41/42", "STG BA22", "IPL BA40", "IFG BA44")

fig <- ggplot(transform(dataset, atlas_loca=factor(atlas_loca,levels=roi_list)), aes(x=atlas_loca, y=data)) + 
          geom_hline(yintercept=0, linetype='dashed', size=size_line, color="#F46D43") +
          geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=data_types), outlier.size = size_jitter_points) +
          geom_jitter(aes(group = data_types), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
          ylab("% power changes & rho values") + 
          scale_y_continuous(breaks=c(-1, 0, 1, 2), labels=c('-100', '0', '100', '200'), limits=c(-1, 2),
                             sec.axis = sec_axis(~ 1*., breaks = c(-1, 0, 1),labels=c('-1', '0', '1'))) +
          scale_x_discrete(breaks=c("STG BA41/42", "STG BA22", "IPL BA40", "IFG BA44"),
                           labels=c(STG_tick, STG_A22_tick, IPL_tick, IFG_tick), ) +
          xlab(element_blank()) + ylab(element_blank()) +
          labs(title=paste("HFa (70-125 Hz)", sep='')) + 
          guides(fill=guide_legend(title=element_blank(), byrow = TRUE), size = guide_legend(order = 3)) +
          theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                axis.text.x = element_text(face="plain",size=size_tick_xaxis),
                axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                legend.text = element_text(face="plain",size=size_legend_text),
                legend.position = "bottom", legend.spacing.x = unit(1, 'cm'))
fig
# left y-axis = power changes (%)
# right y-axis = r


#########################################
## Figure 5 - Phase-Amplitude Coupling ##
#########################################


## Left hemisphere ##

dataset <- read.csv(paste('add path here','R_dataset_2_left.csv',  sep=''))

# choose the phase                 (Virtual Partner or Patient)
reg_choice = "Speech VP"           #or Speech Pat
legend_name = "VP Speech"          #or Patient

dataset = subset(dataset, (regressors == reg_choice |regressors == "Phase Diff"))

dataset_STG_BA4142 = subset(dataset, regions == "STG BA41/42" & ME.test != 0)
dataset_STG_BA22   = subset(dataset, regions == "STG BA22" & ME.test != 0)
dataset_IPL        = subset(dataset, regions == "IPL BA40" & ME.test != 0)
dataset_IFG        = subset(dataset, regions == "IFG BA44" & ME.test != 0)

dataset_STG_BA4142$regressors = factor(dataset_STG_BA4142$regressors, levels=c(reg_choice, 'Phase Diff'))
dataset_STG_BA22$regressors = factor(dataset_STG_BA22$regressors, levels=c(reg_choice, 'Phase Diff'))
dataset_IPL$regressors = factor(dataset_IPL$regressors, levels=c(reg_choice, 'Phase Diff'))
dataset_IFG$regressors = factor(dataset_IFG$regressors, levels=c(reg_choice, 'Phase Diff'))

# plot parameters
width_box = 0.65
size_tick_xaxis = 12
size_tick_yaxis = 12
size_title = 15
size_ylab = 13
alpha_jitter = 0.9
size_jitter_points = 1
size_legend_text = 10
roi_list = c("STG BA41/42", "STG BA22", "IPL BA40", "IFG BA44")

fig_STG_BA4142_pac <- ggplot(dataset_STG_BA4142, aes(x=regressors, y=PAC.values)) +
                        geom_hline(yintercept=0, linetype='dashed', size=0.8, color="#F46D43") +
                        geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=regressors), outlier.size = size_jitter_points) +
                        geom_jitter(aes(group = regressors), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                        scale_x_discrete(breaks=c(reg_choice, 'Phase Diff'), labels=c("Pat \nspeech", "Phase \ndifference")) +
                        scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3), labels=c("0","50","100","150","200", "250", "300")) +
                        xlab(element_blank()) + ylab("PAC (% compared to surrogates mean)") + labs(title="STG BA41/42") + 
                        scale_fill_manual(name = "Regressors", labels=c(legend_name, 'Phase difference'), 
                                          values = c("tomato", "dodgerblue")) +
                        guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                        theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                              axis.text.x = element_blank(),
                              axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                              axis.title.y = element_text(face="plain", size=size_ylab),
                              legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                              legend.position="bottom")

fig_STG_BA22_pac <- ggplot(dataset_STG_BA22, aes(x=regressors, y=PAC.values)) +
                      geom_hline(yintercept=0, linetype='dashed', size=0.8, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=regressors), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = regressors), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      scale_x_discrete(breaks=c(reg_choice, 'Phase Diff'), labels=c("Pat \nspeech", "Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5), labels=c("0","50","100","150","200", "250")) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="STG BA22") + 
                      scale_fill_manual(name = "Regressors", labels=c(legend_name, 'Phase difference'), 
                                        values = c("tomato", "dodgerblue")) +
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig_IPL_BA40_pac <- ggplot(dataset_IPL, aes(x=regressors, y=PAC.values)) +
                      geom_hline(yintercept=0, linetype='dashed', size=0.8, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=regressors), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = regressors), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      scale_x_discrete(breaks=c(reg_choice, 'Phase Diff'), labels=c("Pat \nspeech", "Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 0.9), labels=c("0","50","90"), limits = c(0,0.9)) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="IPL BA40") + 
                      scale_fill_manual(name = "Regressors", labels=c(legend_name, 'Phase difference'), 
                                        values = c("tomato", "dodgerblue")) +
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig_IFG_BA44_pac <- ggplot(dataset_IFG, aes(x=regressors, y=PAC.values)) +
                      geom_hline(yintercept=0, linetype='dashed', size=0.8, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), aes(fill=regressors), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = regressors), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      # scale_x_discrete(breaks=c(reg_choice, "Speech VP", 'Phase Diff'), labels=c("Patient \nspeech", "VP \nspeech", "Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 0.75), labels=c("0","50","75"), limits = c(0,0.75)) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="IFG BA44") + 
                      scale_fill_manual(name = "Regressors", labels=c(legend_name, 'Phase difference'), 
                                        values = c("tomato", "dodgerblue")) + #"cornflowerblue", 
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig <- ggarrange(fig_STG_BA4142_pac, fig_STG_BA22_pac, fig_IPL_BA40_pac, fig_IFG_BA44_pac, ncol = 4, nrow = 1, legend = "bottom", common.legend = TRUE) +
           theme(legend.title=element_text(size=10, face="plain"),
                 legend.text=element_text(size=10,  face="plain"))
fig


###############################################################################################
## Figure S6 - Phase-Amplitude Coupling according to the behavioural delay (left hemisphere) ##
###############################################################################################

dataset <- read.csv(paste('add path here','R_dataset_3.csv',  sep=''))

dataset_STG_BA4142 = subset(dataset, regions == "STG BA41/42" & ME.test != 0)
dataset_STG_BA22   = subset(dataset, regions == "STG BA22"    & ME.test != 0)
dataset_IPL        = subset(dataset, regions == "IPL BA40"    & ME.test != 0)
dataset_IFG        = subset(dataset, regions == "IFG BA44"    & ME.test != 0)

# plot parameters
width_box = 0.65
size_tick_xaxis = 12
size_tick_yaxis = 12
size_title = 15
size_ylab = 13
alpha_jitter = 0.9
size_jitter_points = 1
size_legend_text = 10


fig_STG_BA4142_pac <- ggplot(dataset_STG_BA4142, aes(x=regressors, y=PAC.values, fill = interaction_label)) +
                        geom_hline(yintercept=0, linetype='dashed', size=1.4, color="#F46D43") +
                        geom_boxplot(width=width_box, position=position_dodge(width=0.9), outlier.size = size_jitter_points) +
                        geom_jitter(aes(group = interaction_label), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                        scale_x_discrete(breaks=c('Phase Diff'), labels=c("Phase \ndifference")) +
                        scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2), labels=c("0","50","100","150","200")) +
                        xlab(element_blank()) + ylab("PAC (% compared to surrogates mean)") + labs(title="STG BA41/42") + 
                        scale_fill_manual(name = "Regressors & Delay", 
                                          values = c("dodgerblue", "dodgerblue4")) +
                        guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                        theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                              axis.text.x = element_blank(),
                              axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                              axis.title.y = element_text(face="plain", size=size_ylab),
                              legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                              legend.position="bottom")

fig_STG_BA22_pac <- ggplot(dataset_STG_BA22, aes(x=regressors, y=PAC.values, fill = interaction_label)) +
                      geom_hline(yintercept=0, linetype='dashed', size=1.4, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = interaction_label), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      scale_x_discrete(breaks=c('Phase Diff'), labels=c("Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 1, 1.5), labels=c("0","50","100","150")) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="STG BA22") + 
                      scale_fill_manual(name = "Regressors & Delay", 
                                        values = c("dodgerblue", "dodgerblue4")) +
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig_IPL_BA40_pac <- ggplot(dataset_IPL, aes(x=regressors, y=PAC.values, fill = interaction_label)) +
                      geom_hline(yintercept=0, linetype='dashed', size=1.4, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = interaction_label), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      scale_x_discrete(breaks=c('Phase Diff'), labels=c("Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 0.90), labels=c("0","50","90"), limits = c(0,0.90)) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="IPL BA40") + 
                      scale_fill_manual(name = "Regressors & Delay", 
                                        values = c("dodgerblue", "dodgerblue4")) +
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig_IFG_BA44_pac <- ggplot(dataset_IFG, aes(x=regressors, y=PAC.values, fill = interaction_label)) +
                      geom_hline(yintercept=0, linetype='dashed', size=1.4, color="#F46D43") +
                      geom_boxplot(width=width_box, position=position_dodge(width=0.9), outlier.size = size_jitter_points) +
                      geom_jitter(aes(group = interaction_label), color="grey20", alpha=alpha_jitter, position=position_dodge(width=0.9), size=size_jitter_points) + 
                      # scale_x_discrete(breaks=c("Speech Pat", "Speech VP", 'Phase Diff'), labels=c("Patient \nspeech", "VP \nspeech", "Phase \ndifference")) +
                      scale_y_continuous(breaks=c(0, 0.5, 0.90), labels=c("0","50","90"), limits = c(0,0.90)) +
                      xlab(element_blank()) + ylab(element_blank()) + labs(title="IFG BA44") + 
                      scale_fill_manual(name = "Regressors & Delay", 
                                        values = c("dodgerblue", "dodgerblue4")) + #"cornflowerblue", 
                      guides(fill=guide_legend(title=element_blank(), keywidth = 2, keyheight = 2), size = guide_legend(order = 20)) + 
                      theme(plot.title = element_text(size = size_title, hjust = 0.5,margin = margin(b = 15), face="plain"), 
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(face = "plain",size=size_tick_yaxis),
                            axis.title.y = element_text(face="plain", size=size_ylab),
                            legend.text = element_text(margin = margin(r = 20, unit = "pt"), face="plain",size=size_legend_text),
                            legend.position="bottom")

fig <- ggarrange(fig_STG_BA4142_pac, fig_STG_BA22_pac, fig_IPL_BA40_pac, fig_IFG_BA44_pac, ncol = 4, nrow = 1, legend = "bottom", common.legend = TRUE) +
          theme(legend.title=element_text(size=50, face="plain"), legend.text=element_text(size=50,face="plain"))
fig


