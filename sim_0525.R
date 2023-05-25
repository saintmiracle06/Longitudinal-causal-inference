## Setting simulation: 4 estimators
## sample sizes: 200, 1000, 5000
## DGP 3 time points, L continous, trt=cont.(complex DGP depends on previous), 
## 2-3 time-varying covariates
## 500*simulation
## learner trt and outcome ML models (glmnet, xgboot,rf, mars, gams), hyperparamters tuning
## fold 5-10

library(simcausal)
library(lmtp)

simfun<-function(id = c(200,1000,5000), n_obs =3){
  D <- DAG.empty() + 
    node("b1", distr="rcat.b1", probs = c(0.25,0.75)) +
    node(c("b2","b3"), distr = "mvtnorm::rmvnorm",
         asis.params = list(mean = "c(0,1)")) +
    node("h", t=0, distr="rnorm", mean=3.5, sd = 1) +
    node("l", t=0, distr="rnorm", mean=1.2, sd = 1) + 
    node("z", t=0, distr="rbeta", shape1=1, shape2 =5) +
    
    node("l", t=1:n_obs, distr="rnorm", mean=plogis(1 - 0.3*b1 + 0.4*b2*b3 + 0.2*b1*b2*b3+ 1.5*l[t-1])) +
    node("h", t=1:n_obs, distr="rnorm", mean=plogis(1.5 - 0.5*b1 + 0.3*b1*b3 + 0.3*b1*b2*b3 + 0.5*h[t-1])) +
    node("z", t=1:n_obs, distr="rbeta", shape1=plogis(0.5 - 0.3*b1 + 0.5*h[t]-0.4*l[t] + 0.3*h[t]*l[t]-0.2*b1*b2*b3), 
         shape2=plogis(0.5 - 0.3*b2 + 0.5*h[t]-0.4*l[t] + 0.3*h[t]*l[t]-0.2*b1*b2*b3)) +
    node("Y", t=0:n_obs, distr="rnorm", mean=b1 + 0.3*b2*b3+ 0.9*l[t] + 0.5*h[t] + 0.2*z[t] + 0.3*z[t]*h[t], EFU=TRUE)
  D <- set.DAG(D)
  dat <- sim(D,n=id)
  return(dat)
}

## 500 Simulation sets
for (i in 1:500){
dat<-simfun(id=200)
A<-c("z_0","z_1","z_2","z_3")
B<-c("b1","b2","b3")
L <- list(c("l_0","h_0"),c("l_1","h_1"), c("l_2","h_2"),c("l_3","h_3"))
Y<-c("Y_3")
#SL.library <- c("SL.glmnet", "SL.xgboost","SL.randomForest","SL.glm","SL.polymars", "SL.gam","SL.mean","SL.earth")
SL.library <- c("SL.glmnet","SL.xgboost","SL.randomForest","SL.glm","SL.gam","SL.mean","SL.earth")

## Estimates
## Sequentially doubly robust estimator for the effects of traditional causal effects 
lmtp_sdr(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
         learners_outcome = SL.library,
         learners_trt = SL.library,
         outcome_type = "continuous", mtp = TRUE, folds = 5)

## add causal effect where treatment is set to 1 at all time points for all observations
lmtp_sdr(dat, A, Y, baseline=B, time_vary = L, shift = static_binary_on, 
         learners_outcome = SL.library,
         learners_trt = SL.library,
         outcome_type = "continuous", mtp = TRUE, folds = 5)

## Inverse probability of treatment weighting estimator 
lmtp_ipw(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
         learners_outcome = SL.library,
         learners_trt = SL.library,
         outcome_type = "continuous", mtp = FALSE, folds = 5)
## Standard errors aren't provided for the IPW estimator.

## TMLE estimator 
lmtp_tmle(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
          learners_outcome = SL.library,
          learners_trt = SL.library,
          outcome_type = "continuous", mtp = FALSE, folds = 5)
}


