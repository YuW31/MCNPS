J = K = 5
H = 4
N = 1000


Q1 = matrix(
  c(1,1,0,0,0,
    0,1,0,0,0,
    0,1,1,0,0,
    0,1,0,0,0,
    0,0,1,1,0,
    0,0,1,0,0,
    0,0,0,1,1,
    0,0,0,1,0,
    1,0,0,0,1,
    0,0,0,0,1
  ),
  nrow = 2*J,ncol=K,byrow=T
)
Q1 = cbind(rep(1:J,each = 2),
           rep(1:2,J),
           Q1)
head(Q1)

Q2_B = diag(K)
Q2 = mcq.generate(H,Q2_B)$mcq

LatentClass = class.generate(K)

par.m <- matrix(nrow = 2^K, ncol = 2)
for (m in 1:2^K){
  cat("Class",m,"\n")
  
  true.att = matrix(LatentClass[m,], N, K, byrow = T)
  
  # true.att = stu.generate(1,5,1000)
  
  corr.q = runif(J, 0.77, 0.88)
  mc.prob1 <- prob.generate(Q1,H,true.att,corr.q)
  mc.prob2 <- prob.generate(Q2,H,true.att,corr.q)
  
  obs1.ls = lapply(mc.prob1,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  obs1 <- matrix(unlist(obs1.ls),nrow=N,ncol=J,byrow=TRUE)
  
  obs2.ls = lapply(mc.prob2,function(y){apply(y,1,function(x){sample(1:H,1,prob=x)})})
  obs2 <- matrix(unlist(obs2.ls),nrow=N,ncol=J,byrow=TRUE)
  
  est1 = algo_mc.npc(Q1,obs1,H)$est
  est2 = algo_mc.npc(Q2,obs2,H)$est
  
  par1 = PAR(est1,true.att)
  par2 = PAR(est2,true.att)
  par.m[m,] = c(par1,par2) 
  
}

par.m
