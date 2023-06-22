library(simcausal)
library(lmtp)

id <- 10000000
n_obs <- 3
D <- DAG.empty() + 
  node("b1", distr="rcat.b1", probs = c(0.25,0.75)) +
  node(c("b2","b3"), distr = "mvtnorm::rmvnorm",
       asis.params = list(mean = "c(0,1)")) +
  node("h", t=0, distr="rnorm", mean=3.5, sd = 1) +
  node("l", t=0, distr="rnorm", mean=1.2, sd = 1) + 
  node("z", t=0, distr="rbeta", shape1=1, shape2 =5) +
  node("Y", t=0, distr="rnorm", mean=b1 + 0.2*b2*b3) + 
  
  node("l", t=1:n_obs, distr="rnorm", mean=plogis(1 - 0.3*b1 + 0.4*b2*b3 + 0.2*b1*b2*b3+ 1.5*l[t-1])) +
  node("h", t=1:n_obs, distr="rnorm", mean=plogis(1.5 - 0.5*b1 + 0.3*b1*b3 + 0.3*b1*b2*b3 + 0.5*h[t-1])) +
  node("z", t=1:n_obs, distr="rbeta", shape1=plogis(0.5 - 0.3*b1 + 0.5*h[t]-0.4*l[t] + 0.3*h[t]*l[t]-0.2*b1*b2*b3), 
       shape2=plogis(0.5 - 0.3*b2 + 0.5*h[t]-0.4*l[t] + 0.3*h[t]*l[t]-0.2*b1*b2*b3)) +
  node("Y", t=1:n_obs, distr="rnorm", mean=b1 + 0.3*b2*b3+ 0.9*l[t]-0.65*l[t-1] + 0.5*h[t] +0.4*h[t-1] + 0.2*h[t-1]*l[t-1] 
       + 0.3*z[t]*h[t] + 0.2*z[t-1]*l[t-1], EFU=TRUE)
D <- set.DAG(D)
dat <- sim(D,n=id)

policy <- function(data, trt) {
  a <- data[trt]
  (a - 0.1) * (a - 0.9 > 0) + a * (a - 0.9 < 0)
}


A<-c("z_0","z_1","z_2","z_3")
d<-policy(dat,A)
z_matrix<-dat[A]
Y3<-mean(mean(d$z_0)/mean(z_matrix$z_0)*mean(d$z_1)/mean(z_matrix$z_1)*mean(d$z_2)/mean(z_matrix$z_2)*mean(d$z_2)/mean(z_matrix$z_2)*dat$Y_3)


