### formal simulations for the MC-NPS paper

### Simulation II: the comparisons between the MC-NPS method with the JSD, GNPS, and random selections
### Oct 6, 2023

## Delta Discrimination Index
## No Replication, more conditions
## Compared methods: MC-NPS (with Xu's optimality) and MC-NPS (with new optimality)
## Calibrated
## modified on Sep 26, 2023
setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim3_discussion_newJSD_true")

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


Option.set <- c(4,5,6)
J_Quality.set =  c(1,2)
names(J_Quality.set) <- c("Medium","High")


Matrix.con <- as.matrix(expand.grid(Option.set,J_Quality.set))
# write.csv(Matrix.con,"con.csv",row.names = FALSE)
C.mcnps.con <- c()
p.list = vector(length = nrow(Matrix.con), mode = "list")
# sapply(1:nrow(Matrix.con),function(x){dir.create(paste0("Con",x))})

for (YY in 1:nrow(Matrix.con)) {
  setwd(paste0("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim3_discussion_newJSD_true/Con",YY))
  cat("Condition:", YY,"\n")
  
  H <- Matrix.con[YY,1]
  
  corr.q = ifelse(rep(Matrix.con[YY,2]==1,J),runif(J, 0.77, 0.88),runif(J, 0.88, 1))
  
  all.mcq <- mcq.complete.generate(K,J,H,q.fix,prop_B,Nset_complete)
  q.matrix <- all.mcq$binaryQ
  mc.q <- all.mcq$mcQ
  write.csv(mc.q,paste0("mcq",YY,".csv"))
  key.all <- all.mcq$key
  save.m <- all.mcq$H_star
  
  LatentClass <- class.generate(K)
  eta.class <- epc.generate(mc.q,H,key.all) # "scored" option
  
  ############### Simulate Examinees
  true.att <- stu.generate(ind,K,N)
  mc.prob <- prob.generate(mc.q,H,true.att,corr.q)
  
  IB = c(1:J)
  
  # get all responses for all examinnes
  response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)
  
  process.J <- matrix(nrow=N,ncol=L) #administered item
  
  #============================================================#
  #                         MCNPS-OpMC                         #
  #============================================================#
  
  # Result 1: L(K)-list >> N*K matrix >> estimates of each examinee
  process.est.all_opmc <- vector(length = L,mode = "list")
  process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})
  
  if(H==4){
    
    for (i in 1:N) {
      cat("Condition:", YY,": MC-NPS_NOOP:", i,"\n")
      
      x.n.bank = response.bank[i,,drop=FALSE]
      
      obs.x <- c()
      administer.J <- sample(1:J,1)
      obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
      
      ini.J <- obs.x[,1]
      ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
      candidate.m = ini.est$`distance order`[1:num.c]
      
      process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
      
      rm.J = c(ini.J) # removed items
      remain.J = IB[-rm.J]
      
      for (j in 1:(L-1)) {
        next.J = MCNPS(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m,"Type II")
        
        ### update information
        rm.J = c(rm.J, next.J)
        remain.J = IB[-(rm.J)]
        
        cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
        candidate.m = cat.est$`distance order`[1:num.c]
        process.est.all_opmc[[(j+1)]][i,] = as.matrix(cat.est$est)
        
      }
      
    }
  } else {
    
    J_opmc = ifelse(H==4,K,ceiling(K/2))
    for (i in 1:N) {
      ## for each person, there exists:
      ## initial answer string
      ## remaining items
      ## remaining options
      cat("Condition:", YY,": MC-NPS_OPMC:", i,"\n")
      
      x.n.bank = response.bank[i,,drop=FALSE]
      
      obs.x <- c()
      process.est <- matrix(nrow = L, ncol = K)
      
      for (j in 1:J_opmc){
        # administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
        administer.J <- Q.optimal.mc(mc.q,obs.x,key.all,save.m,H)
        obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
        ini.J <- obs.x[,1]
        ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
        process.est[j, ] = as.matrix(ini.est$est)
        process.est.all_opmc[[j]][i,] = as.matrix(ini.est$est)
      }
      # obs.x <- matrix(c(1,response.bank[i,1],2,response.bank[i,2]),byrow = TRUE,nrow = K)
      
      # ini.J <- obs.x[,1]
      # ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
      candidate.m = ini.est$`distance order`[1:num.c]
      # process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
      
      rm.J = c(ini.J) # removed items
      remain.J = IB[-rm.J]
      
      for (j in 1:(L-J_opmc)) {
        # find the next item
        next.J = MCNPS(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m,"TypeII")
        
        ### update information
        rm.J = c(rm.J, next.J)
        remain.J = IB[-(rm.J)]
        
        cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
        candidate.m = cat.est$`distance order`[1:num.c]
        
        process.est[j+J_opmc, ] = as.matrix(cat.est$est)
        process.est.all_opmc[[(j+J_opmc)]][i,] = as.matrix(cat.est$est)
        
      }
      
    }
    
  }
  
  # true for JSD
  true.par2 = prob.generate(mc.q,H,LatentClass,corr.q)
  true.par = vector(mode = "list", length = J)
  true.par = lapply(true.par,function(x){x=matrix(nrow = H,ncol=2^K)})
  
  for(j in 1:J){
    for (m in 1:2^K){
      true.par[[j]][,m] = true.par2[[m]][j,] 
    }
  }
  
  #============================================================#
  #                            JSD                             #
  #============================================================#
  
  process.est.all_JSD <- vector(length = L,mode = "list")
  process.est.all_JSD <- lapply(process.est.all_JSD,function(x){x = matrix(nrow = N, ncol = K)})
  
  # process.J <- matrix(nrow=N,ncol=L) #administered item

  for (i in 1:N){
    
    cat("Condition:", YY,": JSD:", i,"\n")
    # process.est <- matrix(nrow = L, ncol = K)
    
    administer.J = sample(J,1)
    dat = response.bank[i,administer.J]
    remain.J = c(1:J)[-administer.J]
    prior = 1/(2^K)
    
    for (j in 2:(L+1)){  #+1 for estimation, not actually administered
      # jsd.selection = JSD(mc.q,save.m,O,p.eta = 0.82,dat,administer.J,key) 
      jsd.selection = JSD.func(true.par,prior,dat,administer.J,remain.J)
      process.est.all_JSD[[j-1]][i,] = LatentClass[jsd.selection$Est,]
      
      administer.J = c(administer.J,jsd.selection$Next)
      dat = response.bank[i,administer.J]
      remain.J = c(1:J)[-administer.J]
      prior = jsd.selection$UpPrior
      # process.est[j-1, ] = LatentClass[jsd.selection$Est,]
    }
    
  }
    

  
  #============================================================#
  #                           Random                           #
  #============================================================#
  
  process.est.all_random <- vector(length = L,mode = "list")
  process.est.all_random <- lapply(process.est.all_random,function(x){x = matrix(nrow = N, ncol = K)})
  
  for (i in 1:N){
    cat("Condition:", YY,": Random:", i,"\n")
    
    administer.J = sample(J,L,replace=FALSE)
    dat = response.bank[i,administer.J,drop=FALSE]
    
    for (j in 1:L){ 
      random.est = algo_mc.npc(sub.mcQ(mc.q,administer.J[1:j]),dat[,1:j,drop = FALSE],H)
      process.est.all_random[[j]][i,] = as.matrix(random.est$est)
      
    }
    
  }
  
  ###### ALL CAT ARE DONE
  
  ### PARs across all examinees
  
  process.par.all_opmc <- cbind(1:L,sapply(process.est.all_opmc,function(x){
    PAR(true.att,x)
  }),rep("MC-NPS",L))
  
  process.par.all_jsd<- cbind(1:L,sapply(process.est.all_JSD,function(x){
    PAR(true.att,x)
  }),rep("JSD",L))
  
  process.par.all_r <- cbind(1:L,sapply(process.est.all_random,function(x){
    PAR(true.att,x)
  }),rep("Random",L))
  
  plot.df = data.frame(rbind(process.par.all_opmc, process.par.all_jsd, process.par.all_r))
  colnames(plot.df) = c("Cycle","PAR","Method")
  plot.df$PAR = as.numeric(plot.df$PAR)
  plot.df$Cycle = as.numeric(plot.df$Cycle)
  plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random"))
  write.csv(plot.df,paste0("result_simTRUE_con", YY,".csv"),row.names=FALSE)
  
  colors = c("MC-NPS" = "#000000", "JSD" = "#E69F00", "Random" = "#56B4E9")
  shapes = c("MC-NPS" = 1, "JSD" = 2, "Random" = 3)
  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method)) + #possible change size and colors...
    geom_point(aes(color = Method, shape = Method), size = 2) +
    labs(title = bquote(H[j] == .(H)~"and Item Quality:"~.(names(J_Quality.set)[Matrix.con[YY,2]])), x = "Cycle", y = "PAR") + 
    scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
    scale_x_continuous(breaks = seq(1,L,by=1)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes)
  
  ggsave(paste0("result_simTRUE_con", YY,".png"),p.list[[YY]],width = 10, height = 8) # 8*8, if 3K
  
  
  # ### class size
  # n.class <- subset(PAR.class(true.att,process.est.all_opmc[[1]]),TRUE,ClassSize)
  # write.csv(n.class,paste0("result_sim2_con", YY,"_classN",".csv"))
  # 
  # ### PARs across each class
  # 
  # process.par.class_opmc <- sapply(process.est.all_opmc,function(x){
  #   PAR.class(true.att,x)$PAR
  # })
  # 
  # process.par.class_gnps <- sapply(process.est.all_GNPS,function(x){
  #   PAR.class(true.att,x)$PAR
  # })
  # 
  # process.par.class_jsd<- sapply(process.est.all_JSD,function(x){
  #   PAR.class(true.att,x)$PAR
  # })
  # 
  # process.par.class_r <- sapply(process.est.all_random,function(x){
  #   PAR.class(true.att,x)$PAR
  # })
  # 
  # for (m in 1:2^K){
  #   colors = c("MC-NPS" = "#000000", "JSD" = "#E69F00", "Random" = "#56B4E9", "GNPS" = "#009E73")
  #   shapes = c("MC-NPS" = 1, "JSD" = 2, "Random" = 3, "GNPS" = 7)
  #   
  #   if(n.class[m,]==0){
  #     next
  #   } else{
  #     plot.df = data.frame("Cycle" = rep(1:L,4),
  #                          "Class_PAR" = c(process.par.class_opmc[m,],
  #                                          process.par.class_gnps[m,], 
  #                                          process.par.class_jsd[m,], 
  #                                          process.par.class_r[m,]),
  #                          "Method" = rep(c("MC-NPS","GNPS","JSD","Random"),each = L))
  #     plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random","GNPS"))
  #     
  #     write.csv(plot.df,paste0("result_sim2_con", YY,"_class",m,".csv"),row.names=FALSE)
  #     
  #     if(length(fit_mcdina$item)==0){
  #       plot.df = subset(plot.df, Method != "JSD")
  #       colors = colors[-2]
  #       shapes = shapes[-2]
  #     }
  #     
  #     p = ggplot(data = plot.df, aes(x = Cycle, y = Class_PAR, group = Method)) + 
  #       geom_line(aes(color = Method)) + #possible change size and colors...
  #       geom_point(aes(color = Method, shape = Method), size = 2) +
  #       labs(title = paste0(rownames(n.class)[m], "\n",
  #                           "N0 = ",N0,", Attribute Structure: ", names(dependence_N0.set)[ind_N0], ", H = ", H,
  #                           ", and Item Quality: ", names(J_Quality.set)[Matrix.con[YY,4]]), x = "Cycle", y = "Classwise PAR") + 
  #       scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
  #       scale_x_continuous(breaks = seq(1,L,by=1)) + 
  #       annotate("text", x = 5, y = 1.02, label = paste0("Size of ",rownames(n.class)[m],": ",n.class[m,]), size = 5, color = "Red") +
  #       scale_color_manual(values = colors) +
  #       scale_shape_manual(values = shapes)
  #     
  #     ggsave(paste0("result_sim2_con", YY,"_class",m,".png"),p,width = 10, height = 8) # 8*8, if 3K
  #   }
  #   
  # }
  # 
  # ### PARs across each class that requires the same number of attributes
  # process.par.classK_opmc <- sapply(process.est.all_opmc,function(x){
  #   PAR.class.k(true.att,x)$PAR
  # })
  # 
  # process.par.classK_gnps <- sapply(process.est.all_GNPS,function(x){
  #   PAR.class.k(true.att,x)$PAR
  # })
  # 
  # process.par.classK_jsd<- sapply(process.est.all_JSD,function(x){
  #   PAR.class.k(true.att,x)$PAR
  # })
  # 
  # process.par.classK_r <- sapply(process.est.all_random,function(x){
  #   PAR.class.k(true.att,x)$PAR
  # })
  # 
  # ### class-K size
  # n.classK <- subset(PAR.class.k(true.att,process.est.all_opmc[[1]]),TRUE,ClassSize)
  # write.csv(n.classK,paste0("result_sim2_con", YY,"_classKN",".csv"))
  # 
  # 
  # for (m in 0:K){
  #   colors = c("MC-NPS" = "#000000", "JSD" = "#E69F00", "Random" = "#56B4E9", "GNPS" = "#009E73")
  #   shapes = c("MC-NPS" = 1, "JSD" = 2, "Random" = 3, "GNPS" = 7)
  #   
  #   if(n.classK[m+1,]==0){
  #     next
  #   } else{
  #     plot.df = data.frame("Cycle" = rep(1:L,4),
  #                          "ClassK_PAR" = c(process.par.classK_opmc[m+1,],
  #                                           process.par.classK_gnps[m+1,], 
  #                                           process.par.classK_jsd[m+1,], 
  #                                           process.par.classK_r[m+1,]),
  #                          "Method" = rep(c("MC-NPS","GNPS","JSD","Random"),each = L))
  #     plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random","GNPS"))
  #     
  #     write.csv(plot.df,paste0("result_sim2_con", YY,"_classK",m,".csv"),row.names=FALSE)
  #     
  #     if(length(fit_mcdina$item)==0){
  #       plot.df = subset(plot.df, Method != "JSD")
  #       colors = colors[-2]
  #       shapes = shapes[-2]
  #     }
  #     
  #     p = ggplot(data = plot.df, aes(x = Cycle, y = ClassK_PAR, group = Method)) + 
  #       geom_line(aes(color = Method)) + #possible change size and colors...
  #       geom_point(aes(color = Method, shape = Method), size = 2) +
  #       labs(title = paste0("K",m, "\n",
  #                           "N0 = ",N0,", Attribute Structure: ", names(dependence_N0.set)[ind_N0], ", H = ", H,
  #                           ", and Item Quality: ", names(J_Quality.set)[Matrix.con[YY,4]]), x = "Cycle", y = "Merged Classwise PAR") + 
  #       scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
  #       scale_x_continuous(breaks = seq(1,L,by=1)) + 
  #       annotate("text", x = 5, y = 1.02, label = paste0("Size of K",m,": ",n.classK[m+1,]), size = 5, color = "Red") +
  #       scale_color_manual(values = colors) +
  #       scale_shape_manual(values = shapes)
  #     
  #     ggsave(paste0("result_sim2_con", YY,"_classK",m,".png"),p,width = 10, height = 8) # 8*8, if 3K
  #   }
  # 
  # }
  
}

# setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_true")
# write.csv(C.mcnps.con,"numCal.csv",row.names = TRUE)
# 
# ag <- ggarrange(plotlist = p.list,ncol = 2, nrow = 2,common.legend = TRUE,legend = "bottom",align = "v")
# ggsave(paste0("result_sim2_4_Medium.png"),ag[[1]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_Medium.png"),ag[[2]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_Medium.png"),ag[[3]],width = 16, height = 14)
# ggsave(paste0("result_sim2_4_High.png"),ag[[4]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_High.png"),ag[[5]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_High.png"),ag[[6]],width = 16, height = 14)


