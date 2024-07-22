### PILOT: formal simulations for the MC-NPS paper

### Compare MC-NPS-OPMC, MC-NPS-OPB, MC-OPS-NPOP, MC-OPS-old, JSD; N0 = 100, N = 1000, a = 00000
setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_00000/fixed_op")

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


N0.set <- c(100)
dependence_N0.set <- c(2)
names(dependence_N0.set) <- c("moderately dependent")
Option.set <- c(4,5)
J_Quality.set =  c(1)
names(J_Quality.set) <- c("Medium")


Matrix.con <- as.matrix(expand.grid(N0.set,dependence_N0.set,Option.set,J_Quality.set))
C.mcnps.con <- c()
p.list = vector(length = nrow(Matrix.con), mode = "list")
# sapply(1:nrow(Matrix.con),function(x){dir.create(paste0("Con",x))})

for (YY in 1:nrow(Matrix.con)) {
  cat("Condition:", YY,"\n")
  
  N0 <- Matrix.con[YY,1]
  ind_N0 <- Matrix.con[YY,2]
  H <- Matrix.con[YY,3]
  
  
  all.mcq <- mcq.complete.generate(K,J,H,q.fix,prop_B,Nset_complete)
  q.matrix <- all.mcq$binaryQ
  mc.q <- all.mcq$mcQ
  write.csv(mc.q,paste0("mcq",YY,".csv"))
  key.all <- all.mcq$key
  save.m <- all.mcq$H_star
  
  LatentClass <- class.generate(K)
  eta.class <- epc.generate(mc.q,H,key.all) # "scored" option
  
  ############### Simulate Examinees
  true.att = matrix(c(0,0,0,0,0),nrow = N, ncol = K, byrow = TRUE)
  corr.q = rep(0.85,J)
  mc.prob <- prob.generate(mc.q,H,true.att,corr.q)
  
  IB = c(1:J)
  
  # get all responses for all examinnes
  response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)
  
  i.plot = vector(length = N, mode = "list")
  i.plot = lapply(i.plot,function(x){x=vector(length = 2, mode = "list")})
  
  
  process.est.all_opmc <- vector(length = L,mode = "list")
  process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})
  process.est.all_opmc2 <- vector(length = L,mode = "list")
  process.est.all_opmc2 <- lapply(process.est.all_opmc2,function(x){x = matrix(nrow = N, ncol = K)})
  
  process.J <- matrix(nrow=N,ncol=L) #administered item
  process.J2 <- matrix(nrow=N,ncol=L) #administered item
  
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
    process.est2 <- matrix(nrow = L, ncol = K)
    
    for (j in 1:J_opmc){
      # administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
      administer.J <- Q.optimal.mc(mc.q,obs.x,key.all,save.m,H)
      obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
      ini.J <- obs.x[,1]
      ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
      process.est[j, ] = as.matrix(ini.est$est)
      process.est2[j, ] = as.matrix(ini.est$est)
      process.est.all_opmc[[j]][i,] = as.matrix(ini.est$est)
      process.est.all_opmc2[[j]][i,] = as.matrix(ini.est$est)
    }
    
    obs.x2 = obs.x
    #============================================================#
    #                         MCNPS-OpMC                         #
    #============================================================#
    candidate.m = ini.est$`distance order`[1:num.c]
    # process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
    
    rm.J = c(ini.J) # removed items
    remain.J = IB[-rm.J]
    
    for (j in 1:(L-J_opmc)) {
      # find the next item
      next.J = MCNPS(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m)
      
      ### update information
      rm.J = c(rm.J, next.J)
      remain.J = IB[-(rm.J)]
      
      cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
      candidate.m = cat.est$`distance order`[1:num.c]
      
      process.est[j+J_opmc, ] = as.matrix(cat.est$est)
      process.est.all_opmc[[(j+J_opmc)]][i,] = as.matrix(cat.est$est)
      
    }
    process.J[i,] <- rm.J
    
    
    i.est.plot <- apply(process.est, 1, function(x){
      which(rowSums(abs(matrix(x,nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0)
    })
    i.est.plot <- cbind(1:L,i.est.plot,t(x.n.bank[,rm.J,drop=FALSE]))
    i.est.plot <- data.frame(i.est.plot)
    colnames(i.est.plot) <- c("Cycle","Membership","X")
    i.plot[[i]][[1]] = ggplot(i.est.plot, aes(x = Cycle, y = Membership, label = X)) +
      geom_point(size = 3, color = "blue") +
      labs(title = paste0("MC-NPS with OPMC: ","Examinee: ",i), x = "Cycle", y = "Membership") +
      scale_y_continuous(limits = c(1, 2^K), breaks = seq(1, 2^K, by = 1)) + 
      scale_x_continuous(breaks = seq(1, L, by = 1)) + 
      geom_hline(yintercept  = which(rowSums(abs(matrix(true.att[i,],nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0), color = "red")+
      geom_text(size = 5, nudge_x = 0.3, nudge_y = 0.5)
    # ggsave(paste0("MCNPS_", YY, "_",i,".png"),p1,width = 9, height = 9)
    
    write.csv(sub.mcQ(mc.q,rm.J),paste0("J_MCNPS_",YY,"_",i,".csv"),row.names = FALSE)
    
    #============================================================#
    #                    MCNPS-OpMC-Old Index                    #
    #============================================================#
    cat("Condition:", YY,": MC-NPS_Old:", i,"\n")
    

    
    candidate.m2 = ini.est$`distance order`[1:num.c]
    # process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
    
    rm.J2 = c(ini.J) # removed items
    remain.J2 = IB[-rm.J2]
    
    for (j in 1:(L-J_opmc)) {
      # find the next item
      next.J2 = MCNPS(LatentClass[candidate.m2[1],],LatentClass[candidate.m2[2],],eta.class,remain.J2,H,save.m,"Type II")
      
      ### update information
      rm.J2 = c(rm.J2, next.J2)
      remain.J2 = IB[-(rm.J2)]
      
      cat.est2 <- algo_mc.npc(sub.mcQ(mc.q,rm.J2),x.n.bank[,rm.J2,drop=FALSE],H)
      candidate.m2 = cat.est2$`distance order`[1:num.c]
      
      process.est2[j+J_opmc, ] = as.matrix(cat.est2$est)
      process.est.all_opmc2[[(j+J_opmc)]][i,] = as.matrix(cat.est2$est)
      
    }
    process.J2[i,] <- rm.J2
    
    
    i.est.plot <- apply(process.est2, 1, function(x){
      which(rowSums(abs(matrix(x,nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0)
    })
    i.est.plot <- cbind(1:L,i.est.plot,t(x.n.bank[,rm.J2,drop=FALSE]))
    i.est.plot <- data.frame(i.est.plot)
    colnames(i.est.plot) <- c("Cycle","Membership","X")
    i.plot[[i]][[2]] = ggplot(i.est.plot, aes(x = Cycle, y = Membership, label = X)) +
      geom_point(size = 3, color = "blue") +
      labs(title = paste0("MC-NPS with OPMC_Old: ","Examinee: ",i), x = "Cycle", y = "Membership") +
      scale_y_continuous(limits = c(1, 2^K), breaks = seq(1, 2^K, by = 1)) + 
      scale_x_continuous(breaks = seq(1, L, by = 1)) + 
      geom_hline(yintercept  = which(rowSums(abs(matrix(true.att[i,],nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0), color = "red")+
      geom_text(size = 5, nudge_x = 0.3, nudge_y = 0.5)
    # ggsave(paste0("MCNPS_", YY, "_",i,".png"),p1,width = 9, height = 9)
    
    write.csv(sub.mcQ(mc.q,rm.J2),paste0("J_MCNPS2_",YY,"_",i,".csv"),row.names = FALSE)
    
  }
  
  
  write.csv(process.J,paste0("J_MCNPS_",YY,".csv"))
  write.csv(process.J2,paste0("J_MCNPS2_",YY,".csv"))
  
  
  
  process.par.all_opmc <- cbind(1:L,sapply(process.est.all_opmc,function(x){
    PAR(true.att,x)
  }),rep("Optimality_MC",L))
  
  process.par.all_mc2 <- cbind(1:L,sapply(process.est.all_opmc2,function(x){
    PAR(true.att,x)
  }),rep("Old",L))
  
  # process.par.all_r <- cbind(1:L,sapply(process.est.all_random,function(x){
  #   PAR(true.att,x)
  # }),rep("Random",L))
  
  plot.df = data.frame(rbind(process.par.all_opmc, process.par.all_mc2))
  # process.par.all_r))
  colnames(plot.df) = c("Cycle","PAR","Method")
  plot.df$PAR = as.numeric(plot.df$PAR)
  plot.df$Cycle = as.numeric(plot.df$Cycle)
  plot.df$Method = factor(plot.df$Method, levels = c("Optimality_MC", "Old"))
  # "Random","GNPS"))
  
  write.csv(plot.df,paste0("result_sim2_con", YY,".csv"),row.names=FALSE)

  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method)) + #possible change size and colors...
    labs(title = paste0("N0 = ",N0,", Attribute Structure: ", names(dependence_N0.set)[ind_N0], ", H = ", H,
                        ", and Item Quality: ", names(J_Quality.set)[Matrix.con[YY,4]]), x = "Cycle", y = "PAR") + 
    scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
    scale_x_continuous(breaks = seq(1,L,by=1))
  
  ggsave(paste0("result_sim2_con", YY,".png"),p.list[[YY]],width = 14, height = 8) # 8*8, if 3K
  
  for (i in 1:N){
    ag <- ggarrange(plotlist = i.plot[[i]],ncol = 2, nrow = 2,common.legend = TRUE,legend = "bottom",align = "v")
    ggsave(paste0("Con_", YY, "_",i,".png"),ag,width = 16, height = 16) # 16*8, if 3K
  }
  
  
}

setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_00000/fixed_op")
ag <- ggarrange(plotlist = p.list,ncol = 2, nrow = 1,common.legend = TRUE,legend = "bottom",align = "v")
ggsave(paste0("result_sim2.png"),ag,width = 26, height = 8) # 16*8, if 3K

write.csv(C.mcnps.con,"numCal.csv",row.names = TRUE)

