### investigate the response trajectory

### Initial Findings (without running the investigation, just looking results by eyes):
# 1. For 0000, it seems more items can reduce the number of mastered skills;
#    The ideal pattern: 1100 to 1000 to 0000 (0000-302)


setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/discrimination2")

mc.q <- as.matrix(read.csv("simQ.csv", header = TRUE))
head(mc.q)
K = 4
L = 3*K
H = 5

LatentClass <- class.generate(K)

save.m = c()
binary.q = matrix(nrow = 300,ncol = K)
key = c()
for (j in 1:300){
  qj = mc.q[which(mc.q[,1]==j),,drop=FALSE]
  save.m = c(save.m,nrow(qj))
  key = c(key,qj[1,2])
  binary.q[j,] = qj[1,-c(1,2)]
}

eta.class <- epc.generate(mc.q,H,key.all) # "scored" option

i0.set <- c(9,10,20)
length(i0.set)

# What we want to take a look:
# 1. Whether 0000 is considered as the second possible profile? Specifically, I don't understand how the latent class jumps to 0010 from 1010?
# 2. How did the distance between each latent class change?
# 3. Why the estimate is stuck? whether more items can solve this issue? [79/176] If so, should we use a combined stopping rule? Fixed length + 3 stable estimates
# 4. Based on the current estimation, which q-vector is the best to .... Does this q-vector exist in the Q-matrix? Is this selected?
# 5. Case 126: abnormal
# 6. [166] If the estimated profile is 0000 before 3L, the following items cannot secure the estimated profile
# 183: 
# 300: longer test, mastering more attributes?
# 303: In which case, the examinee will be estimated to master more and more skill

administer.J <- as.matrix(read.csv("Optimal_J_administered_1-0_0_0_0.csv", header = TRUE))
head(administer.J)

setwd("/Users/yuwang/Desktop/Research/MC-NPS/nps/discrimination2/m1")
response.bank <- as.matrix(read.csv("response_1-0_0_0_0.csv", header = TRUE))
response.bank <- response.bank[,-1]
head(response.bank)

### begin iinvestigation
i = 3

### loop begins here
paste0("This is the ",i0.set[i], "th examinee belonging to the latent class 0000.")

J.i <- administer.J[i0.set[i],-1]

dist.i <- matrix(nrow = 2^K,ncol = L)
for (j in 1:L){
  cat("\n","\n","Item",j,"\n")
  Q.j <- sub.mcQ(mc.q,J.i[1:j])
  print(Q.j[which(Q.j[,1]==J.i[j]),])
  
  est.j <- algo_mc.npc(Q.j,response.bank[i0.set[i],J.i[1:j],drop=FALSE],H)
  
  print(response.bank[i0.set[i],J.i[1:j],drop=FALSE])
  print(est.j$est)
  consider.m = est.j$`distance order`[1:2]
  cat("\n","Next cycle will consider:",consider.m,"\n")
  print(est.j$distance)

  dist.i[,j] <- est.j$distance
}


which(apply(binary.q,1,function(x){sum(abs(x-LatentClass[3,]))})==0)

v1 = LatentClass[10,]
v2 = LatentClass[15,]
etas = eta.class
j2 = (1:300)[-c(1,8,189)]
O = H
Hj_star = save.m
ws = ""


dis.P = Delta(v1,v2,etas,j2,O,Hj_star,ws = "")
dist.max = max(dis.P)
dist.max
candidate.J = which(dis.P == dist.max)
j2[candidate.J]

mc.q[which(mc.q[,1]==32),]


i = i + 1

### verify my explanations
