# plot for Simulation Study I

library(MASS)
library(ggplot2)
library(ggpubr)

J = 300
K = 5
N = 1000 #MVN

q.fix <- 0.5
prop_B = 0.2
Nset_complete = 5
L = K  # in 
num.c = 2 #only consider two candidates
ind <- 2 #MVN

Option.set <- c(4,5,6)
J_Quality.set =  c(1,2)
names(J_Quality.set) <- c("Medium","High")

Matrix.con <- as.matrix(expand.grid(Option.set,J_Quality.set))
p.list = vector(length = length(Option.set), mode = "list")

setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim1_formal")
for (YY in 1:nrow(Matrix.con)) {
  cat("Condition:", YY)
  
  
  H <- Matrix.con[YY,1]
  corr.q = ifelse(rep(Matrix.con[YY,2]==1,J),runif(J, 0.77, 0.88),runif(J, 0.88, 1))
  
  colors = c("MC Optimality" = "#D55E00", "Binary Optimality" = "#0072B2", "MCNPS Without Optimality" = "#009E73")
  shapes = c("MC Optimality" = 16, "Binary Optimality" = 17, "MCNPS Without Optimality" = 15)
  # Shape 16: A filled circle; Shape 17: A filled diamond; Shape 15: A filled square.
  
  plot.df = read.csv(paste0("result_sim1_H", H,"Q",Matrix.con[YY,2],".csv"), header = T)
  plot.df$Method = factor(plot.df$Method, levels = c("Optimality_MC","Optimality_B","Without Optimality"))
  levels(plot.df$Method) <- c("MC Optimality", "Binary Optimality", "MCNPS Without Optimality")
  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method)) + #possible change size and colors...
    geom_point(aes(color = Method, shape = Method), size = 3) +
    labs(title = bquote(H[j] == .(H)), #, " & Item Quality: ",names(J_Quality.set)[Matrix.con[YY,2]]), 
                        x = "Cycle", y = "PAR") + 
    scale_y_continuous(breaks = seq(0,1,by=0.1), limits = c(0, 0.85)) + 
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    geom_vline(xintercept = ceiling(K/2), linetype = "dotted", color = "red", size = 1.2) +
    geom_vline(xintercept = K, linetype = "dotted", color = "#00CC00", size = 1.2) +
    theme_minimal() +  # Set a minimalistic theme with a clear background
    theme(text = element_text(size=20),
          legend.text = element_text(size = 25)) +
    theme(panel.grid = element_line(color = "gray", linetype = "dotted"))  # Add gridlines
  
}

setwd("/Users/yuwang/Downloads")
ag <- ggarrange(plotlist = p.list,ncol = 3, nrow = 1,common.legend = T,legend = "bottom")
ggsave(paste0("result_sim1_Medium.png"),ag[[1]],width = 20, height = 8)
ggsave(paste0("result_sim1_High.png"),ag[[2]],width = 20, height = 8)

# an_ag = annotate_figure(ag[[1]],top = text_grob("Item Quality: Medium", hjust = 2.715,size = 26))
# ggsave(paste0("result_sim1_Medium.png"),an_ag,width = 20, height = 8)

# an_ag2 = annotate_figure(ag[[2]],top = text_grob("Item Quality: High", hjust = 3.22,size = 26))
# top = text_grob("Item Quality: High", color = "black", face = "plain", size = 26))
# ggsave(paste0("result_sim1_High.png"),an_ag2,width = 20, height = 8)






