### RUN Yigit's conditions

prob.generate.step = function(mcq, O, att, corr.q) {
  # decreased by 0.05 for each "step down"
  
  item.no <- unique(mcq[, 1])
  J = length(item.no)
  K = ncol(mcq) - 2
  
  N = nrow(att)
  
  LatentClass = class.generate(K)
  
  Ideal <- LatentClass %*% t(mcq[, 3:(2 + K)])
  Ideal.met <- 1 * (Ideal == (matrix(1, 2^K) %*% t(rowSums(mcq[, 3:(2 + K)]))))
  
  
  prob <- vector(length = 2^K, mode = "list")
  for (i in 1:2^K) {
    # 0.06; 0.82; 0.25: same as Jimmy
    #op <- mc.q[,1:2][which(Ideal.conj[i,]>0),,drop=FALSE]
    sub.r <- matrix(1 / O, nrow = J, ncol = O)
    
    i.eta <- Ideal.met[i, ]
    
    for (j in 1:J) {
      q = corr.q[j]
      
      work.q <- mcq[which(mcq[, 1] == item.no[j]), , drop = FALSE]
      work.eta <- i.eta[which(mcq[, 1] == item.no[j])]
      
      if (work.eta[1] == 1) {
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, work.q[work.eta[1], 2]] = q
      } else if (sum(work.eta) == 1){
        q = q-(which(work.eta == 1)-1)*0.05
        
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, work.q[which(work.eta == 1), 2]] = q
      } else if (sum(work.eta) > 1) { # not apply 0.05 to this condition
        possible.improper.Q = work.q[work.eta == 1, -(1:2)]
        
        ori.order <-
          order(rowSums(possible.improper.Q), decreasing = TRUE)
        all.include = all(apply(possible.improper.Q, 1, function(row)
          all(row <= possible.improper.Q[ori.order[1], ])))
        
        if (all.include) {
          work.eta[which(work.eta == 1)[-ori.order[1]]] = 0
        }
        
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, work.q[which(work.eta == 1), 2]] = (q + (1 - q) * (length(which(
          work.eta == 1
        )) - 1) / (O - 1)) / length(which(work.eta == 1))
      }
      
    }
    
    prob[[i]] <- sub.r
  }
  
  
  prob.each <- vector(length = N, mode = "list")
  for (i in 1:N){
    prob.each[[i]] = prob[[which(apply(LatentClass,1,function(y){
      all.equal(att[i,],y)
    })==TRUE)]]
  }
  
  return(prob.each)
  
}

mcq.generate.nest = function(O, Q) {
  mc.q <- c()
  key <- c()
  save.m <- c()
  
  J <- nrow(Q)
  K <- ncol(Q)
  q.matrix <- cbind(Q, rowSums(Q))
  
  for (i in 1:J) {
    num <- q.matrix[i, K + 1]
    
    sub.q <- matrix(0, nrow = num, ncol = K)
    od.option <- sample(1:O, O)
    id.key <- od.option[1]
    sub.q[1, ] <- q.matrix[i, 1:K]
    
    if (num > 1){
      for (m in 2:num){
        r.q <- which(sub.q[m-1,] == 1)
        x1 = sample(r.q,1)
        v <- sub.q[m-1,]
        v[x1] = 0
        sub.q[m, ] <- v
      }
        
    }
    
    item.q <- cbind(rep(i, num), od.option[1:num], sub.q)
    
    mc.q <- rbind(mc.q, item.q)
    key <- c(key, id.key)
    save.m <- c(save.m, num)
    
  }
  
  rownames(mc.q) <- NULL
  colnames(mc.q) <- c("Item", "Option", paste0("Att", seq(1:K)))
  
  return(list(
    "mcq" = mc.q,
    "key" = key,
    "H_star" = save.m
  ))
}

setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim4_Yigit")

library(MASS)
library(cdcatR)
library(ggplot2)
library(ggpubr)
library(GDINA)
library(cdcatR)
library(cdmTools)

J = 600
K = 5
N = 6000 #MVN
L = 4*K #stopping rule
H = 5
num.c = 2


# c(1,2,7,17,27,32)

LatentClass <- class.generate(K)
true.att <- t(cbind(replicate(1000,LatentClass[1,]),
                  replicate(1000,LatentClass[2,]),
                  replicate(1000,LatentClass[7,]),
                  replicate(1000,LatentClass[17,]),
                  replicate(1000,LatentClass[27,]),
                  replicate(1000,LatentClass[32,])))

discrimination.set <- c(0.8,0.6)
variance.set <- c(0.05,0.025)


# q.matrix = rbind(do.call("rbind",replicate(J/K/sum(rowSums(LatentClass)==1),LatentClass[rowSums(LatentClass)==1,],simplify = FALSE)),
#                  do.call("rbind",replicate(J/K/sum(rowSums(LatentClass)==2),LatentClass[rowSums(LatentClass)==2,],simplify = FALSE)),
#                  do.call("rbind",replicate(J/K/sum(rowSums(LatentClass)==3),LatentClass[rowSums(LatentClass)==3,],simplify = FALSE)),
#                  do.call("rbind",replicate(J/K/sum(rowSums(LatentClass)==4),LatentClass[rowSums(LatentClass)==4,],simplify = FALSE)),
#                  do.call("rbind",replicate(J/K/sum(rowSums(LatentClass)==5),LatentClass[rowSums(LatentClass)==5,],simplify = FALSE))
# )

q.matrix = t(cbind(replicate(J/5,LatentClass[2,]),
                   replicate(J/5,LatentClass[7,]),
                   replicate(J/5,LatentClass[17,]),
                   replicate(J/5,LatentClass[27,]),
                   replicate(J/5,LatentClass[32,])))
mc.q = rbind(t(replicate(J/5,LatentClass[2,])),
             do.call("rbind",replicate(J/5,LatentClass[c(7,2),],simplify = FALSE)),
             do.call("rbind",replicate(J/5,LatentClass[c(17,7,2),],simplify = FALSE)),
             do.call("rbind",replicate(J/5,LatentClass[c(27,17,7,2),],simplify = FALSE)),
             do.call("rbind",replicate(J/5,LatentClass[c(32,27,17,7,2),],simplify = FALSE)))

mc.q = cbind(c(1:(J/5),rep((J/5+1):(2*J/5),each = 2),rep((2*J/5+1):(3*J/5),each = 3),rep((3*J/5+1):(4*J/5),each = 4),rep((4*J/5+1):J,each = 5)),
             c(rep(1,J/5),rep(c(1,2),J/5),rep(c(1,2,3),J/5),rep(c(1:4),J/5),rep(c(1:5),J/5)),
             mc.q)
key = rep(1,J)
save.m = c(rep(1,J/5),rep(2,J/5),rep(3,J/5),rep(4,J/5),rep(5,J/5))
write.csv(mc.q,"ex_mcQ.csv",row.names = FALSE)


Matrix.con <- as.matrix(expand.grid(discrimination.set,variance.set))
p.list = vector(length = nrow(Matrix.con), mode = "list")
LS = rbind(LatentClass[1,],LatentClass[2,],
          LatentClass[7,],LatentClass[17,],
          LatentClass[27,],LatentClass[32,])

for (YY in 1:nrow(Matrix.con)) {
  # setwd(paste0("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_formal_new/Con",YY))
  cat("Condition:", YY,"\n")
 
  dis.p = Matrix.con[YY,1]
  var.p = Matrix.con[YY,2]

 
  eta.class <- epc.generate(mc.q,H,key.all,LS) # "scored" option
  
  corr.q = runif(J, dis.p-var.p, dis.p+var.p)
  
  ############### Simulate Examinees
  mc.prob <- prob.generate.step(mc.q,H,true.att,corr.q)
  
  IB = c(1:J)
  
  # get all responses for all examinnes
  response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)
  
  process.J <- matrix(nrow=N,ncol=L) #administered item
  
  #============================================================#
  #                  MCNPS-NOOp (bc the Q-matrix)              #
  #============================================================#
  
  # Result 1: L(K)-list >> N*K matrix >> estimates of each examinee
  process.est.all_opmc <- vector(length = L,mode = "list")
  process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})
  
    
  for (i in 1:N) {
    ## for each person, there exists:
    ## initial answer string
    ## remaining items
    ## remaining options
    cat("Condition:", YY,": MC-NPS_NOOP:", i,"\n")
    
    x.n.bank = response.bank[i,,drop=FALSE]
    
    ini.J = sample(J,1)
    ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H,LS)
    candidate.m = ini.est$`distance order`[1:num.c]
    process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
    
    rm.J = c(ini.J) # removed items
    remain.J = IB[-rm.J]
    
    for (j in 2:L) {
      # find the next item
      next.J = MCNPS(LS[candidate.m[1],],LS[candidate.m[2],],eta.class,remain.J,H,save.m,"TypeII",LS)
      
      ### update information
      rm.J = c(rm.J, next.J)
      remain.J = IB[-(rm.J)]
      
      cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H,LS)
      candidate.m = cat.est$`distance order`[1:num.c]
      
      process.est.all_opmc[[j]][i,] = as.matrix(cat.est$est)
      
    }
    
  }
  
  
  
  # true for JSD
  true.par2 = prob.generate(mc.q,H,LatentClass,corr.q)
  true.par = vector(mode = "list", length = J)
  true.par = lapply(true.par,function(x){x=matrix(nrow = H,ncol=6)})
  
  for(j in 1:J){
    for (i in 1:6){
      m = c(1,2,7,17,27,32)[i]
      true.par[[j]][,i] = true.par2[[m]][j,] 
    }
  }
  
  
  #============================================================#
  #                    JSD with true parameters                #
  #============================================================#
  
  process.est.all_JSD <- vector(length = L,mode = "list")
  process.est.all_JSD <- lapply(process.est.all_JSD,function(x){x = matrix(nrow = N, ncol = K)})
  
    
   
  item.par = true.par
    
  for (i in 1:N){
    
    cat("Condition:", YY,": JSD:", i,"\n")
    # process.est <- matrix(nrow = L, ncol = K)
    
    administer.J = sample(J,1)
    dat = response.bank[i,administer.J]
    remain.J = c(1:J)[-administer.J]
    prior = rep(1/6,6)
    
    save.list = vector(mode = "list",length = 19)
    
    for (j in 2:(L+1)){  #+1 for estimation, not actually administered
      # jsd.selection = JSD(mc.q,save.m,O,p.eta = 0.82,dat,administer.J,key) 
      jsd.selection = JSD.func(item.par,prior,dat,administer.J,remain.J)
      save.list[[j-1]] = jsd.selection
      process.est.all_JSD[[j-1]][i,] = LatentClass[c(1,2,7,17,27,32)[jsd.selection$Est],]
      
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
      random.est = algo_mc.npc(sub.mcQ(mc.q,administer.J[1:j]),dat[,1:j,drop = FALSE],H,LS)
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
  
  plot.df = data.frame(rbind(process.par.all_opmc,process.par.all_jsd, process.par.all_r))
  colnames(plot.df) = c("Cycle","PAR","Method")
  plot.df$PAR = as.numeric(plot.df$PAR)
  plot.df$Cycle = as.numeric(plot.df$Cycle)
  plot.df$Method = factor(plot.df$Method, levels = c("MC-NPS","JSD","Random"))
  write.csv(plot.df,paste0("result_sim4_con", YY,".csv"),row.names=FALSE)
  
  colors = c("MC-NPS" = "#000000", "JSD" = "#E69F00", "Random" = "#56B4E9")
  shapes = c("MC-NPS" = 1, "JSD" = 2, "Random" = 3)
  
  
  p.list[[YY]] = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
    geom_line(aes(color = Method)) + #possible change size and colors...
    geom_point(aes(color = Method, shape = Method), size = 2) +
    labs(title = paste0("Dis = ",Matrix.con[YY,1],"Var = ",Matrix.con[YY,2]), x = "Cycle", y = "PAR") + 
    scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
    scale_x_continuous(breaks = seq(1,L,by=1)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes)
  
  ggsave(paste0("result_sim4_con", YY,".png"),p.list[[YY]],width = 10, height = 8) # 8*8, if 3K
  
  
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
  #                                          process.par.classK_gnps[m+1,], 
  #                                          process.par.classK_jsd[m+1,], 
  #                                          process.par.classK_r[m+1,]),
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

# setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_formal_new")
# write.csv(C.mcnps.con,"numCal.csv",row.names = TRUE)
# 
# ag <- ggarrange(plotlist = p.list,ncol = 2, nrow = 2,common.legend = TRUE,legend = "bottom",align = "v")
# ggsave(paste0("result_sim2_4_Medium.png"),ag[[1]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_Medium.png"),ag[[2]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_Medium.png"),ag[[3]],width = 16, height = 14)
# ggsave(paste0("result_sim2_4_High.png"),ag[[4]],width = 16, height = 14)
# ggsave(paste0("result_sim2_5_High.png"),ag[[5]],width = 16, height = 14) 
# ggsave(paste0("result_sim2_6_High.png"),ag[[6]],width = 16, height = 14)


