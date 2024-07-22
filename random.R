# the first item is selected randomly, and then the discrimination_random power is used
# so no Q-optimality rule was involved
## modified on Oct 4, 2023

library(gtools)
library(plyr)
library(MASS)
library(cdcatR)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(data.table)
setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/discrimination_random")

N <- 2000
K <- 4
J <- c(300)

H = 5
q.fix <- 0.5   ### how to set q in MC

L = 3*K  ### stopping rule: fixed length
num.c = 2 #only consider two candidate items

# q.matrix <- generate.Q(K,J,q.fix)
# all.mcq <- mcq.generate(H,q.matrix)
# mc.q <- all.mcq$mcq
# key.all <- all.mcq$key
# save.m <- all.mcq$`coded options`

all.mcq <- mcq.complete.generate(K,J,H,q.fix)
q.matrix <- all.mcq$binaryQ
mc.q <- all.mcq$mcQ
key.all <- all.mcq$key
save.m <- all.mcq$`coded options`

write.csv(mc.q,"simQ.csv",row.names = FALSE)

LatentClass <- class.generate(K)
eta.class <- epc.generate(mc.q,H,key.all) # "scored" option

J_star = K/2

for (m in 1:nrow(LatentClass)){
  cat("LatentClass:", m)
  # sapply(1:16,function(x){dir.create(paste0("m",x))})
  setwd(paste0("/Users/yuwang/Desktop/Research/MC-NPS/nps/discrimination_random/m",m))
  
  true.att <- matrix(LatentClass[m,],nrow=N,ncol=K,byrow = TRUE)
  name <- paste0(LatentClass[m,],collapse = "_")
  corr.q = runif(J,0.99-1/H, 0.9) # gs parameters
  # corr.q = rep(0.9,J)
  
  mc.prob <- prob.generate(mc.q,H,true.att,corr.q)
  
  ### No need for iterations for now
  ## only one test
  
  # initial status: K/2 or (K+1)/2 items
  # same for all examinees
  IB = c(1:J)
  
  # get all responses for all examinnes
  etaID.response = etaID.generate(mc.q,H,key.all)[1:2,m]
  
  response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)
  write.csv(response.bank,paste0("response_",m,"-", name,".csv"))
  # response.bank[,1] <- etaID.response[1]
  # response.bank[,2] <- etaID.response[2]
  
  process.est.person <- vector(length = N, mode = "list")
  process.est.all <- vector(length = L,mode = "list")
  
  process.est.all <- lapply(process.est.all,function(x){x = matrix(nrow = N, ncol = K)})
  
  process.J <- matrix(nrow=N,ncol=L) #administered item
  
  for (i in 1:N) {
    ## for each person, there exists:
    ## initial answer string
    ## remaining items
    ## remaining options
    cat("MC-NPS-Examinee:", i,"\n")
    
    x.n.bank = response.bank[i,,drop=FALSE]
    
    obs.x <- c()
    administer.J <- sample(1:J,1)
    obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
    
    ini.J <- obs.x[,1]
    ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
    candidate.m = ini.est$`distance order`[1:num.c]
    
    process.est <- as.matrix(ini.est$est)
    process.est.all[[1]][i,] = process.est
    
    rm.J = c(ini.J) # removed items
    remain.J = IB[-rm.J]
    
    for (j in 1:(L-1)) {
      # find the next item
      dist.set = Delta(LatentClass[candidate.m[1],],LatentClass[candidate.m[2],],eta.class,remain.J,H,save.m)
      dist.max = max(dist.set)
      candidate.J = which(dist.set == dist.max)
      next.J = remain.J[ifelse(length(candidate.J)==1,candidate.J,sample(candidate.J))]
      
      ### update information
      rm.J = c(rm.J, next.J)
      remain.J = IB[-(rm.J)]
      
      cat.est <- algo_mc.npc(sub.mcQ(mc.q,rm.J),x.n.bank[,rm.J,drop=FALSE],H)
      candidate.m = cat.est$`distance order`[1:num.c]
      process.est <- rbind(process.est,as.matrix(cat.est$est))
      process.est.all[[(j+1)]][i,] = as.matrix(cat.est$est)
      
    }
    
    process.est.person[[i]] <- process.est
    
    i.est.plot <- apply(process.est, 1, function(x){
      which(rowSums(abs(matrix(x,nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0)
    })
    i.est.plot <- cbind(1:L,i.est.plot)
    i.est.plot <- data.frame(i.est.plot)
    colnames(i.est.plot) <- c("Cycle","Membership")
    
    p1 = ggplot(i.est.plot, aes(x = Cycle, y = Membership)) +
      geom_point(size = 3, color = "blue") +
      labs(title = paste0("(Optimal) Estimation Path:
                          ", "Group: ", m," Examinee: ",i), x = "Cycle", y = "Membership") +
      scale_y_continuous(limits = c(1, 16), breaks = seq(1, 16, by = 1)) + 
      scale_x_continuous(breaks = seq(1, L, by = 1)) + 
      geom_hline(yintercept  = m, color = "red")
    ggsave(paste0("OpEP_", m,"_p",i,".png"),p1,width = 9, height = 9)
    
    process.J[i,] <- rm.J
    
  }
  
  setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/discrimination_random")
  
  p.all <- vector(length = length(process.est.all), mode = "list")
  str.class <- apply(LatentClass, 1, paste, collapse = ",")
  for (j in 1:length(process.est.all)) {
    est_all <- apply(process.est.all[[j]], 1, paste, collapse = ",")
    n.class <- sapply(str.class,function(x){
      length(which(est_all==x))
    })
    plot.df <- data.frame("class" = str.class,
                          "prop" = n.class/N)
    p.all[[j]] = ggplot(data = plot.df, aes(x = class, y = prop)) + 
      geom_point() +
      geom_hline(yintercept  = 1/2^K, color = "red") +
      # labs(title = paste0("(Optimal) Estimation Path:", "Group: ", m," Item: ",(j+1)), x = "Cycle", y = "Membership") +
      labs(title = paste0("(Optimal) Estimation Path:", "Group: ", m," Item: ",(j)), x = "Cycle", y = "Membership") +
      scale_y_continuous(limits = c(0, 1))
  }
  ag <- ggarrange(plotlist = p.all,ncol = 3, nrow = 4)
  ggsave(paste0("OpEP_", m,"All",".png"),ag,width = 20, height = 20)
  
  
  write.csv(process.J,paste0("Optimal_J_administered_",m,"-", name,".csv"))
  
}




