### formal simulations for the MC-NPS paper

### Simulation I: the efficiency of the Q-optimality for MC items
### Compare PAR after the first K items: Q-optimality for MC items, Q-optimality, random1 + discrimination power, completely random
### Oct 5, 2023
library(MASS)
library(ggplot2)
library(ggpubr)

setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim1_formal")
# setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps")

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

for (YY in 1:nrow(Matrix.con)) {
  cat("Condition:", YY)
  
  
  H <- Matrix.con[YY,1]
  corr.q = ifelse(rep(Matrix.con[YY,2]==1,J),runif(J, 0.77, 0.88),runif(J, 0.88, 1))
  
  all.mcq <- mcq.complete.generate(K,J,H,q.fix,prop_B,Nset_complete)
  q.matrix <- all.mcq$binaryQ
  mc.q <- all.mcq$mcQ
  key.all <- all.mcq$key
  save.m <- all.mcq$H_star
  
  LatentClass <- class.generate(K)
  eta.class <- epc.generate(mc.q,H,key.all) # "scored" option
  
  ############### Simulate Examinees
  true.att <- stu.generate(ind,K,N)
  # true.att <- matrix(LatentClass[1,],nrow=N,ncol=K,byrow = TRUE)
  mc.prob <- prob.generate(mc.q,H,true.att,corr.q)

  # same for all examinees
  IB = c(1:J)
  
  # get all responses for all examinnes
  response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)
  
  # Condition I: Q-optimality for MC items
  # Result 1: L(K)-list >> N*K matrix >> estimates of each examinee
  process.est.all_opmc <- vector(length = K,mode = "list")
  process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})
  # 
  # process.J <- matrix(nrow=N,ncol=L) #administered item
  
  J_opmc = ifelse(H==4,K,ceiling(K/2))
  for (i in 1:N) {
    ## for each person, there exists:
    ## initial answer string
    ## remaining items
    ## remaining options
    cat("MC-OP:", i,"\n")
    
    x.n.bank = response.bank[i,,drop=FALSE]
    
    obs.x <- c()
    
    for (j in 1:J_opmc){
      # administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
      administer.J <- Q.optimal.mc(mc.q,obs.x,key.all,save.m,H) #Q.optimal.mc(mc.q,obs.x,key.all,save.m,H)
      obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
      ini.J <- obs.x[,1]
      ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
      process.est.all_opmc[[j]][i,] = as.matrix(ini.est$est)
    }
    # obs.x <- matrix(c(1,response.bank[i,1],2,response.bank[i,2]),byrow = TRUE,nrow = K)
    
    # ini.J <- obs.x[,1]
    # ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
    candidate.m = ini.est$`distance order`[1:num.c]
    # process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
    
    rm.J = c(ini.J) # removed items
    remain.J = IB[-rm.J]
    
    if(H==4){
      next
    } else {
      for (j in 1:(L-J_opmc)) {
        # find the next item
        next.J = MCNPS(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m)
        
        ### update information
        rm.J = c(rm.J, next.J)
        remain.J = IB[-(rm.J)]
        
        cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
        candidate.m = cat.est$`distance order`[1:num.c]
        
        process.est.all_opmc[[(j+J_opmc)]][i,] = as.matrix(cat.est$est)
        
      }
    }
    
  }
  
  
  # Condition II: Q-optimality for binary items
  # Result 2: L(K)-list >> N*K matrix >> estimates of each examinee
  process.est.all_opb <- vector(length = K,mode = "list")
  process.est.all_opb <- lapply(process.est.all_opb,function(x){x = matrix(nrow = N, ncol = K)})
  
  J_opb = K
  for (i in 1:N) {
    cat("B-OP:", i,"\n")
    
    x.n.bank = response.bank[i,,drop=FALSE]
    
    obs.x <- c()
    for (j in 1:J_opb){
      # administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
      administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
      obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
      ini.J <- obs.x[,1]
      ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
      process.est.all_opb[[j]][i,] = as.matrix(ini.est$est)
    }
    
  }
  
  
  # Condition III: NO Q-optimality for MC items
  # Result e: L(K)-list >> N*K matrix >> estimates of each examinee
  process.est.all_noop <- vector(length = K,mode = "list")
  process.est.all_noop <- lapply(process.est.all_noop,function(x){x = matrix(nrow = N, ncol = K)})
  # 
  # process.J <- matrix(nrow=N,ncol=L) #administered item
  
  for (i in 1:N) {
    cat("NO-OP:", i,"\n")
    
    x.n.bank = response.bank[i,,drop=FALSE]
    
    obs.x <- c()
    administer.J <- sample(1:J,1)
    obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
    
    ini.J <- obs.x[,1]
    ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
    candidate.m = ini.est$`distance order`[1:num.c]
    
    process.est.all_noop[[1]][i,] = as.matrix(ini.est$est)
    
    rm.J = c(ini.J) # removed items
    remain.J = IB[-rm.J]
    
    for (j in 1:(L-1)) {
      next.J = MCNPS(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m)
      
      ### update information
      rm.J = c(rm.J, next.J)
      remain.J = IB[-(rm.J)]
      
      cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
      candidate.m = cat.est$`distance order`[1:num.c]
      process.est.all_noop[[(j+1)]][i,] = as.matrix(cat.est$est)
      
    }
  
  }
  
  process.par.all_opmc <- cbind(1:L,sapply(process.est.all_opmc,function(x){
    PAR(true.att,x)
  }),rep("Optimality_MC",5))
  
  process.par.all_opb <- cbind(1:L,sapply(process.est.all_opb,function(x){
    PAR(true.att,x)
  }),rep("Optimality_B",5))
  
  process.par.all_noop <- cbind(1:L,sapply(process.est.all_noop,function(x){
    PAR(true.att,x)
  }),rep("Without Optimality",5))
  
  plot.df = data.frame(rbind(process.par.all_opmc, process.par.all_opb, process.par.all_noop))
  colnames(plot.df) = c("Cycle","PAR","Method")
  plot.df$PAR = as.numeric(plot.df$PAR)
  plot.df$Method = factor(plot.df$Method, levels = c("Optimality_MC", "Optimality_B", "Without Optimality"))
  write.csv(plot.df,paste0("result_sim1_H", H,"Q",Matrix.con[YY,2],".csv"),row.names=FALSE)

  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method),size = 1) + #possible change size and colors...
    labs(title = paste0("H = ", H, " & Item Quality: ",names(J_Quality.set)[Matrix.con[YY,2]]), x = "Cycle", y = "PAR") + ylim(0,1)
  ggsave(paste0("result_sim1_H", H,"_Q",Matrix.con[YY,2],".png"),p.list[[YY]],width = 8, height = 8)
  
}

ag <- ggarrange(plotlist = p.list,ncol = 3, nrow = 1,common.legend = TRUE,legend = "bottom",align = "v")
ggsave(paste0("result_sim1_Medium.png"),ag[[1]],width = 20, height = 8)
ggsave(paste0("result_sim1_High.png"),ag[[2]],width = 20, height = 8)




