### PILOT: formal simulations for the MC-NPS paper

### Compare MC-NPS-OPMC, MC-NPS-OPB, MC-OPS-NPOP, MC-OPS-old, JSD; N0 = 100, N = 1000, a = 00000
setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/sim2_K45_fixedOP")

library(MASS)
library(cdcatR)
library(ggplot2)
library(ggpubr)
library(GDINA)
library(cdcatR)
library(cdmTools)

J = 300
K = 5
N = 2000 #MVN
ind <- 2 #MVN
q.fix <- 0.5
prop_B = 0.2
Nset_complete = 5
num.c = 2 #only consider two candidates
L = 4*K #stopping rule


H = 5
# J_Quality.set =  c(1)
# names(J_Quality.set) <- c("Medium")
corr.q = runif(J, 0.77, 0.88)

all.mcq <- mcq.complete.generate(K,J,H,q.fix,prop_B,Nset_complete)
q.matrix <- all.mcq$binaryQ
mc.q <- all.mcq$mcQ
write.csv(mc.q,"mcq.csv")
key.all <- all.mcq$key
save.m <- all.mcq$H_star

LatentClass <- class.generate(K)
eta.class <- epc.generate(mc.q,H,key.all) # "scored" option

############### Simulate Examinees
true.att <- matrix(LatentClass[30,], N, K, byrow = T)
mc.prob <- prob.generate(mc.q,H,true.att,corr.q)

IB = c(1:J)

C.mcnps.con <- c()
C.mcnps = 1
conv = FALSE
ind_N0 = 2
N0 = 100

while (!conv && C.mcnps <= 10000){
  true.att_N0 <- stu.generate(ind_N0,K,N0)
  mc.prob_N0 <- prob.generate(mc.q,H,true.att_N0,corr.q)
  
  N0.response.ls = lapply(mc.prob_N0,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  N0.response <- matrix(unlist(N0.response.ls),nrow=N0,ncol=J,byrow=TRUE)
  
  fit_mcdina = MCmodel_Yu(N0.response, mc.q, model = "MCDINA", key = key.all, H,conv.crit = .001,maxitr=2000,conv.type="pr",SE=FALSE)
  conv = !(fit_mcdina$person[1,1]==9)
  C.mcnps = C.mcnps + 1
  cat("Generating the calibration sample",C.mcnps,"\n")
}
C.mcnps.con <- c(C.mcnps.con, C.mcnps)



# get all responses for all examinnes
response.bank.ls = lapply(mc.prob,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
response.bank <- matrix(unlist(response.bank.ls),nrow=N,ncol=J,byrow=TRUE)

#============================================================#
#                         MCNPS-OpMC                         #
#============================================================#

# Result 1: L(K)-list >> N*K matrix >> estimates of each examinee
process.est.all_opmc <- vector(length = L-2,mode = "list")
process.est.all_opmc <- lapply(process.est.all_opmc,function(x){x = matrix(nrow = N, ncol = K)})
process.J <- matrix(nrow=N,ncol=L) #administered item

process.est.all_JSD <- vector(length = L-2,mode = "list")
process.est.all_JSD <- lapply(process.est.all_JSD,function(x){x = matrix(nrow = N, ncol = K)})
process.J2 <- matrix(nrow=N,ncol=L) #administered item

# process.J <- matrix(nrow=N,ncol=L) #administered item

J_opmc = ifelse(H==4,K,ceiling(K/2))
for (i in 1:N) {
  ## for each person, there exists:
  ## initial answer string
  ## remaining items
  ## remaining options
  cat("MC-NPS_OPMC:", i,"\n")
  
  x.n.bank = response.bank[i,,drop=FALSE]
  
  obs.x <- c()
  process.est <- matrix(nrow = L-2, ncol = K)
  
  for (j in 1:J_opmc){
    # administer.J <- Q.optimal.binary(q.matrix,obs.x,key.all)
    administer.J <- Q.optimal.mc(mc.q,obs.x,key.all,save.m,H)
    obs.x <- rbind(obs.x,c(administer.J,x.n.bank[,administer.J]))
    ini.J <- obs.x[,1]
    ini.est <- algo_mc.npc(sub.mcQ(mc.q,ini.J),x.n.bank[,ini.J,drop=FALSE],H)
  }
  
  process.est[1, ] = as.matrix(ini.est$est)
  process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
  
  # process.est.all_opmc[[1]][i,] = as.matrix(ini.est$est)
  
  candidate.m = ini.est$`distance order`[1:num.c]
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
    
    process.est[j+1, ] = as.matrix(cat.est$est)
    process.est.all_opmc[[(j+1)]][i,] = as.matrix(cat.est$est)
    
  }
  process.J[i,] <- rm.J

  
  i.est.plot <- apply(process.est, 1, function(x){
    which(rowSums(abs(matrix(x,nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0)
  })
  i.est.plot <- cbind(1:(L-2),i.est.plot,t(x.n.bank[,rm.J[-c(1:2)],drop=FALSE]))
  i.est.plot <- data.frame(i.est.plot)
  colnames(i.est.plot) <- c("Cycle","Membership","X")
  p1 = ggplot(i.est.plot, aes(x = Cycle, y = Membership, label = X)) +
    geom_point(size = 3, color = "blue") +
    labs(title = paste0("MC-NPS with OPMC: ","Examinee: ",i), x = "Cycle", y = "Membership") +
    scale_y_continuous(limits = c(1, 2^K), breaks = seq(1, 2^K, by = 1)) + 
    scale_x_continuous(breaks = seq(1, L-2, by = 1)) + 
    geom_hline(yintercept  = which(rowSums(abs(matrix(true.att[i,],nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0), color = "red")+
    geom_text(size = 5, nudge_x = 0.3, nudge_y = 0.5)
  ggsave(paste0("MCNPS_",i,".png"),p1,width = 9, height = 9)
  
  write.csv(sub.mcQ(mc.q,rm.J),paste0("J_MCNPS_",i,".csv"),row.names = FALSE)
  
  #============================================================#
  #                            JSD                             #
  #============================================================#
  
  cat("JSD:", i,"\n")
  process.est2 = matrix(nrow = L-2, ncol = K)
  administer.J = ini.J
  dat = response.bank[i,administer.J]
  
  for (j in (1+J_opmc):(L+1)){  #+1 for estimation, not actually administered
    # jsd.selection = JSD(mc.q,save.m,O,p.eta = 0.82,dat,administer.J,key) 
    jsd.selection = JSD(mc.q,save.m,H,p.eta = fit_mcdina$item,dat,administer.J,key.all)
    process.est.all_JSD[[j-J_opmc]][i,] = LatentClass[jsd.selection$Est,]
    
    administer.J = c(administer.J,jsd.selection$Next)
    dat = response.bank[i,administer.J]
    process.est2[j-J_opmc, ] = LatentClass[jsd.selection$Est,]
  }
  
  i.est.plot <- apply(process.est2, 1, function(x){
    which(rowSums(abs(matrix(x,nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0)
  })
  i.est.plot <- cbind(1:(L-2),i.est.plot, dat[-c(1,2,L+1)])
  i.est.plot <- data.frame(i.est.plot)
  colnames(i.est.plot) <- c("Cycle","Membership","X")
  p1 = ggplot(i.est.plot, aes(x = Cycle, y = Membership,label = X)) +
    geom_point(size = 3, color = "blue") +
    labs(title = paste0("JSD: ","Examinee: ",i), x = "Cycle", y = "Membership") +
    scale_y_continuous(limits = c(1, 2^K), breaks = seq(1, 2^K, by = 1)) + 
    scale_x_continuous(breaks = seq(1, L-2, by = 1)) + 
    geom_hline(yintercept  = which(rowSums(abs(matrix(true.att[i,],nrow = 2^K, ncol = K,byrow = TRUE)-LatentClass))==0), color = "red") +
    geom_text(size = 5, nudge_x = 0.3, nudge_y = 0.5)
  ggsave(paste0("JSD_", i,".png"),p1,width = 9, height = 9)
  
  write.csv(sub.mcQ(mc.q,administer.J[-(L+1)]),paste0("J_JSD_",i,".csv"),row.names = FALSE)
  
  process.J2[i,] <- administer.J[-(L+1)]
  
  
  
}
write.csv(process.J,"J_MCNPS.csv")
write.csv(process.J2,paste0("J_JSD.csv"))


#=============#
# Calibration Sample for GNPS and JSD
#=============#
# Since the MC-NPS may not be able to converge, I may need to generate several samples to calibrate the item parameters for the JSD method

process.par.all_opmc <- cbind(1:(L-J_opmc+1),sapply(process.est.all_opmc,function(x){
  PAR(true.att,x)
}),rep("Optimality_MC",L-J_opmc+1))

process.par.all_jsd<- cbind(1:(L-J_opmc+1),sapply(process.est.all_JSD,function(x){
  PAR(true.att,x)
}),rep("JSD",L-J_opmc+1))


plot.df = data.frame(rbind(process.par.all_opmc, process.par.all_jsd))
colnames(plot.df) = c("Cycle","PAR","Method")
plot.df$PAR = as.numeric(plot.df$PAR)
plot.df$Cycle = as.numeric(plot.df$Cycle)
plot.df$Method = factor(plot.df$Method, levels = c("Optimality_MC", "JSD"))

write.csv(plot.df,paste0("result_sim2_K45.csv"),row.names=FALSE)


p = ggplot(data = plot.df, aes(x = Cycle, y = PAR, group = Method)) + 
  geom_line(aes(color = Method)) + #possible change size and colors...
  labs(title = "10111", x = "Cycle", y = "PAR") + 
  scale_y_continuous(breaks = seq(0,1,by=0.1)) + 
  scale_x_continuous(breaks = seq(1,L,by=1))

ggsave("result_sim2_K45.png",p,width = 14, height = 8) # 8*8, if 3K






