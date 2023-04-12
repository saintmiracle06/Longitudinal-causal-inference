library("simcausal")

D <- DAG.empty() + 
  node("b", distr="rcat.b1", probs = c(0.25,0.75)) +
  node("l", t=0, distr="rcat.b1", probs = c(0.20,0.80)) +
  node("h", t=0, distr="rcat.b1", probs = c(0.81,0.19)) + 
  node("z", t=0, distr="rbern", prob=plogis(1 - 0.3*b + 0.5*h[t] -0.3*l[t])) +
  
  node("l", t=1:3, distr="rbern", prob=plogis(1- 0.3*b + 1.5*l[t-1])) +
  node("h", t=1:3, distr="rbern", prob=plogis(1.5 - 0.5*b + 0.5*h[t-1])) +
  node("z", t=1:3, distr="rbern", prob=plogis(0.5 - 0.3*b + 0.5*h[t]-0.4*l[t])) +
  node("d", t=0:3, distr="rbern", prob=plogis(-1+0.5*z[t])) +
  node("Yhat", t=0:3, distr="rnorm", mean=b + 0.9*l[t] + 0.5*h[t] + 0.2*d[t], EFU=TRUE) +
  node("Y", t=0:3, distr="rnorm", mean=b + 0.9*l[t] + 0.5*h[t] + 0.2*z[t], EFU=TRUE)

D <- set.DAG(D)

dat.long <- sim(D,n=100000)
colnames(dat.long)
(ATE1<-mean(dat.long$Yhat_0))
##ATE1:4.019579
(ATE2<-mean(dat.long$Yhat_1))
##ATE2:3.052384
(ATE3<-mean(dat.long$Yhat_2))
##ATE3:2.951065
(ATE4<-mean(dat.long$Yhat_3))
##ATE4:2.925953



