library(simcausal)
library(lmtp)

id<-1000
#number of visits (K+1)
n_obs<-2
set.seed(1256)
D <- DAG.empty() + 
  node("b", distr="rcat.b1", probs = c(0.25,0.75)) +
  node("l", t=0, distr="rcat.b1", probs = c(0.20,0.80)) +
  node("h", t=0, distr="rnorm", mean=1.2) + 
  node("z", t=0, distr="rbern", prob=plogis(1 - 0.3*b + 0.5*h[t] -0.3*l[t])) +
  
  node("l", t=1:n_obs, distr="rbern", prob=plogis(1 - 0.3*b + 1.5*l[t-1])) +
  node("h", t=1:n_obs, distr="rbern", prob=plogis(1.5 - 0.5*b + 0.5*h[t-1])) +
  node("z", t=1:n_obs, distr="rbern", prob=plogis(0.5 - 0.3*b + 0.5*h[t]-0.4*l[t])) +
  node("Y", t=0:n_obs, distr="rnorm", mean=b + 0.9*l[t] + 0.5*h[t] + 0.2*z[t], EFU=TRUE)

D <- set.DAG(D)

dat <- sim(D,n=id)
dat_long<-reshape(data = dat,sep = "_",,direction="long",idvar="ID",varying = -c(1:2))

## Estimates
A<-c("z_0","z_1","z_2")
B<-c("b")
L <- list(c("l_0","h_0"),c("l_1","h_1"), c("l_2","h_2"))
Y<-c("Y_1")
SL.library <- c("SL.glmnet", "SL.glm", "SL.mean","SL.earth")

## Sequentially doubly robust estimator for the effects of traditional causal effects 
lmtp_sdr(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
         learners_outcome = SL.library,
         learners_trt = "SL.glmnet",
         outcome_type = "continuous", mtp = TRUE, folds = 2)


## Inverse probability of treatment weighting estimator 
lmtp_ipw(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
         learners_outcome = SL.library,
         learners_trt = "SL.glmnet",
         outcome_type = "continuous", mtp = FALSE, folds = 2)

## TMLE estimator 
lmtp_tmle(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
          learners_outcome = SL.library,
          learners_trt = "SL.glmnet",
          outcome_type = "continuous", mtp = FALSE, folds = 2)
