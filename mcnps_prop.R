### 
setwd("/Users/yuwang/Downloads")

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

H <- 5

corr.q = runif(J, 0.77, 0.88)

all.mcq <- mcq.complete.generate(K,J,H,q.fix,prop_B,Nset_complete)
q.matrix <- all.mcq$binaryQ
mc.q <- all.mcq$mcQ
write.csv(mc.q,paste0("mcq_m.csv"))
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
process.Hj_star <- matrix(nrow=N,ncol=L) #administered item

#============================================================#
#                         MCNPS-OpMC                         #
#============================================================#

# Result 1: L(K)-list >> N*K matrix >> estimates of each examinee
process.est.all_opmc <- vector(length = L,mode = "list")
process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})

if(H==4){
  
  for (i in 1:N) {
    cat("MC-NPS_OPMC:", i,"\n")
    
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
    cat("MC-NPS_OPMC:", i,"\n")
    
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
    
    process.J[i,] = rm.J
    process.Hj_star[i,] = save.m[rm.J]
    
    
  }
  
}


write.csv(process.J,"J.csv",row.names = FALSE)
write.csv(process.Hj_star,"Hj_star.csv",row.names = FALSE)

write.csv(do.call(rbind,apply(process.Hj_star,2,table))/1000,"Hj_star_Prop.csv",row.names = FALSE)
