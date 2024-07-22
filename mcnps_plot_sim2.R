# plot for Simulation Study II
library(MASS)
library(cdcatR)
library(ggplot2)
library(ggpubr)
library(GDINA)
library(cdcatR)
library(cdmTools)

J = 300
K = 5
N = 1000 #MVN
ind <- 2 #MVN
q.fix <- 0.5
prop_B = 0.2
Nset_complete = 5
num.c = 2 #only consider two candidates
L = 4*K #stopping rule


N0.set <- c(30,100,500)
dependence_N0.set <- c(1,2)
names(dependence_N0.set) <- c("Independent","Moderately Dependent")
Option.set <- c(4,5,6)
J_Quality.set =  c(1,2)
names(J_Quality.set) <- c("Medium","High")


Matrix.con <- as.matrix(expand.grid(N0.set,dependence_N0.set,Option.set,J_Quality.set))
C.mcnps.con <- c()
p.list = vector(length = nrow(Matrix.con), mode = "list")
p2.list = vector(length = nrow(Matrix.con), mode = "list")
# sapply(1:nrow(Matrix.con),function(x){dir.create(paste0("Con",x))})

for (YY in 1:nrow(Matrix.con)) {
  setwd(paste0("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_formal_new_30&100&500"))
  cat("Condition:", YY,"\n")
  
  N0 <- Matrix.con[YY,1]
  ind_N0 <- Matrix.con[YY,2]
  H <- Matrix.con[YY,3]
  
  plot.df = read.csv(paste0("result_sim2_con", YY,".csv"))
  plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random","GNPS"))
  
  colors = c("MC-NPS" = "#D55E00", "JSD" = "#0072B2", "Random" = "#E69F00", "GNPS" = "#009E73")
  shapes = c("MC-NPS" = 15, "JSD" = 16, "Random" = 17, "GNPS" = 8)
  
  JSD.calibration = is.na(plot.df$PAR[which(plot.df$Method=="JSD")[1]])
  
  if(JSD.calibration){
    plot.df = subset(plot.df, Method != "JSD")
    colors = colors[-2]
    shapes = shapes[-2]
  }
  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method)) + #possible change size and colors...
    geom_point(aes(color = Method, shape = Method), size = 2) +
    labs(title = bquote(N[0] == .(N0)~"and"~.(names(dependence_N0.set)[ind_N0])~"Attributes"), x = "Cycle", y = "PAR") + 
    # labs(title = expression(paste(N_0," = ",N0,", and ", names(dependence_N0.set)[ind_N0], "structure"
                                   #, ", H = ", H,", and Item Quality: ", names(J_Quality.set)[Matrix.con[YY,4]])

    scale_y_continuous(breaks = seq(0,1,by=0.1), limits = c(0, 1)) + 
    scale_x_continuous(breaks = seq(1,L,by=1)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    theme_minimal() +  # Set a minimalistic theme with a clear background
    theme(text = element_text(size=20),
          legend.text = element_text(size = 25)) +
    theme(panel.grid = element_line(color = "gray", linetype = "dotted"))  # Add gridlines
  
  if (H > 4){
    p.list[[YY]] = p.list[[YY]] + geom_vline(xintercept = ceiling(K/2), linetype = "dotted", color = "red", size = 1.2)
  }
  
  
}

setwd("/Users/yuwang/Downloads")

#### 24 conditions: 30 and 100
# ag <- ggarrange(plotlist = p.list,ncol = 2, nrow = 2,common.legend = TRUE,legend = "bottom")
# 
# ggsave(paste0("result_sim2_4_Medium.png"),ag[[1]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_Medium.png"),ag[[2]],width = 16, height = 14)
# ggsave(paste0("result_sim2_6_Medium.png"),ag[[3]],width = 16, height = 14)
# ggsave(paste0("result_sim2_4_High.png"),ag[[4]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_High.png"),ag[[5]],width = 16, height = 14)
# ggsave(paste0("result_sim2_6_High.png"),ag[[6]],width = 16, height = 14)

### 36 conditions: 30, 100, and 500
ag <- ggarrange(plotlist = p.list,ncol = 3, nrow = 2,common.legend = TRUE,legend = "bottom")

ggsave(paste0("result_sim2_4_Medium.png"),ag[[1]],width = 24, height = 14)
ggsave(paste0("result_sim2_5_Medium.png"),ag[[2]],width = 24, height = 14)
ggsave(paste0("result_sim2_6_Medium.png"),ag[[3]],width = 24, height = 14)
ggsave(paste0("result_sim2_4_High.png"),ag[[4]],width = 24, height = 14)
ggsave(paste0("result_sim2_5_High.png"),ag[[5]],width = 24, height = 14)
ggsave(paste0("result_sim2_6_High.png"),ag[[6]],width = 24, height = 14)


# an_ag1 = annotate_figure(ag[[1]],top = text_grob("H = 4 & Item Quality: Medium", hjust = 1.53,size = 26))
# ggsave(paste0("result_sim2_4_Medium.png"),an_ag1,width = 16, height = 14)
# 
# an_ag2 = annotate_figure(ag[[2]],top = text_grob("H = 5 & Item Quality: Medium", hjust = 1.53,size = 26))
# ggsave(paste0("result_sim2_5_Medium.png"),an_ag2,width = 16, height = 14)
# 
# an_ag3 = annotate_figure(ag[[3]],top = text_grob("H = 6 & Item Quality: Medium", hjust = 1.53,size = 26))
# ggsave(paste0("result_sim2_6_Medium.png"),an_ag3,width = 16, height = 14)
# 
# an_ag4 = annotate_figure(ag[[4]],top = text_grob("H = 4 & Item Quality: High", hjust = 1.73,size = 26))
# ggsave(paste0("result_sim2_4_High.png"),an_ag4,width = 16, height = 14)
# 
# an_ag5 = annotate_figure(ag[[5]],top = text_grob("H = 5 & Item Quality: High", hjust = 1.73,size = 26))
# ggsave(paste0("result_sim2_5_High.png"),an_ag5,width = 16, height = 14)
# 
# an_ag6 = annotate_figure(ag[[6]],top = text_grob("H = 6 & Item Quality: High", hjust = 1.73,size = 26))
# ggsave(paste0("result_sim2_6_High.png"),an_ag6,width = 16, height = 14)


# # for Class K
# m = 0
# colors = c("MC-NPS" = "#D55E00", "JSD" = "#0072B2", "Random" = "#E69F00", "GNPS" = "#009E73")
# shapes = c("MC-NPS" = 15, "JSD" = 16, "Random" = 17, "GNPS" = 8)
# 
# for (YY in 1:nrow(Matrix.con)) {
#   
#   setwd(paste0("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_formal/Con",YY))
#   cat("Condition:", YY,"\n")
#   
#   N0 <- Matrix.con[YY,1]
#   ind_N0 <- Matrix.con[YY,2]
#   H <- Matrix.con[YY,3]
#   
#   n.classK <- read.csv(paste0("result_sim2_con", YY,"_classKN",".csv"))
#   
#   if(n.classK[m+1,2]==0){
#     next
#   } else{
#     plot.df = read.csv(paste0("result_sim2_con", YY,"_classK",m,".csv"))
#     plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random","GNPS"))
#     JSD.calibration = plot.df$ClassK_PAR[which(plot.df$Method=="JSD")[1]]==0
#     
#     if(JSD.calibration){
#       plot.df = subset(plot.df, Method != "JSD")
#       colors = colors[-2]
#       shapes = shapes[-2]
#     }
#     
#     p2.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = ClassK_PAR, group = Method)) + 
#       geom_line(aes(color = Method)) + #possible change size and colors...
#       geom_point(aes(color = Method, shape = Method), size = 2) +
#       labs(title = paste0("K",m, "\n",
#                           "N0 = ",N0,", Attribute Structure: ", names(dependence_N0.set)[ind_N0], ", H = ", H,
#                           ", and Item Quality: ", names(J_Quality.set)[Matrix.con[YY,4]]), x = "Cycle", y = "Merged Classwise PAR") + 
#       scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
#       scale_x_continuous(breaks = seq(1,L,by=1)) + 
#       annotate("text", x = 5, y = 1.02, label = paste0("Size of K",m,": ",n.classK[m+1,2]), size = 5, color = "Red") +
#       scale_color_manual(values = colors) +
#       scale_shape_manual(values = shapes)
#   }
#   
#   
# }
# 
# setwd("/Users/yuwang/Downloads")
# ag2 <- ggarrange(plotlist = p2.list,ncol = 2, nrow = 2,common.legend = TRUE,legend = "bottom",align = "v")
# ggsave(paste0("result_sim2_4_Medium_",m,".png"),ag2[[1]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_Medium_",m,".png"),ag2[[2]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_Medium_",m,".png"),ag2[[3]],width = 16, height = 14)
# ggsave(paste0("result_sim2_4_High_",m,".png"),ag2[[4]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_High_",m,".png"),ag2[[5]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_High_",m,".png"),ag2[[6]],width = 16, height = 14)

