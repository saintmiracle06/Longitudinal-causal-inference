## Setting simulation: 4 estimators
## sample sizes: 200, 1000, 5000
## DGP 3 time points, L continous, trt=cont.(complex DGP depends on previous), 
## 2-3 time-varying covariates
## 500*simulation
## learner trt and outcome ML models (glmnet, xgboot,rf, mars, gams), hyperparamters tuning
## fold 5-10

## true value 10 million

library(simcausal)
library(lmtp)

simfun<-function(id = c(200,1000,5000), n_obs =3)
  {
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
  return(dat)
}

##DGP turn formula into functions
## Bias of the estimators

## 500 Simulation sets
est1<-rep(NA,500)
est2<-rep(NA,500)
est3<-rep(NA,500)
est4<-rep(NA,500)

sd1<-rep(NA,500)
sd2<-rep(NA,500)
sd3<-rep(NA,500)
sd4<-rep(NA,500)

low1<-rep(NA,500)
low2<-rep(NA,500)
low3<-rep(NA,500)
low4<-rep(NA,500)

high1<-rep(NA,500)
high2<-rep(NA,500)
high3<-rep(NA,500)
high4<-rep(NA,500)

##[beta distribution (0,1)]
## check policy
## page 14, longitudinal treatment effect ##
policy <- function(data, trt) {
  a <- data[[trt]]
  (a - 0.1) * (a - 0.9 > 0) + a * (a - 0.9 < 0)
}

SL.library <- c("SL.glmnet","SL.xgboost","SL.randomForest","SL.glm","SL.mean","SL.earth")

for (i in 1:500)
{
dat<-simfun(id=200)
A<-c("z_0","z_1","z_2","z_3")
B<-c("b1","b2","b3")
L <- list(c("l_0","h_0","Y_0"),c("l_1","h_1","Y_1"), c("l_2","h_2","Y_2"),c("l_3","h_3","Y_3"))
Y<-c("Y_3")
#SL.library <- c("SL.glmnet", "SL.xgboost","SL.randomForest","SL.glm","SL.polymars", "SL.gam","SL.mean","SL.earth")
#SL.library <- c("SL.glmnet","SL.xgboost","SL.randomForest","SL.glm","SL.mean","SL.earth")



## Estimates
## Sequentially doubly robust estimator for the effects of traditional causal effects 
#lmtp_sdr(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
#         learners_outcome = SL.library,
#         learners_trt = SL.library,
#         outcome_type = "continuous", mtp = TRUE, folds = 5)

## add causal effect where treatment is set to 1 at all time points for all observations
a<-lmtp_sdr(dat, A, Y, baseline=B, time_vary = L, shift = policy, 
         learners_outcome = SL.library,
         learners_trt = SL.library,
         outcome_type = "continuous", mtp = TRUE, folds = 5)
est1[i]<-a$theta
sd1[i]<-a$standard_error
low1[i]<-a$low
high1[i]<-a$high

## Inverse probability of treatment weighting estimator 
b<-lmtp_ipw(dat, A, Y, baseline=B, time_vary = L, shift = policy, 
         learners_outcome = SL.library,
         learners_trt = SL.library,
         outcome_type = "continuous", mtp = FALSE, folds = 5)
## Standard errors aren't provided for the IPW estimator.
est2[i]<-b$theta

## TMLE estimator 
c<-lmtp_tmle(dat, A, Y, baseline=B, time_vary = L, shift = policy, 
          learners_outcome = SL.library,
          learners_trt = SL.library,
          outcome_type = "continuous", mtp = FALSE, folds = 5)
est3[i]<-c$theta
sd3[i]<-c$standard_error
low3[i]<-c$low
high3[i]<-c$high

## G-computation estimator
d<-lmtp_sub(dat, A, Y, baseline=B, time_vary = L, shift = NULL, 
         learners = SL.library,
         outcome_type = "continuous", folds = 5)
est4[i]<-d$theta
print(i)

}

write.csv(est1, "est_sdr_200.csv", row.names=FALSE)
write.csv(est2, "est_ipw_200.csv", row.names=FALSE)
write.csv(est3, "est_tmle_200.csv", row.names=FALSE)
write.csv(est4, "est_sub_200.csv", row.names=FALSE)

write.csv(sd1, "sd_sdr_200.csv", row.names=FALSE)
write.csv(sd2, "sd_ipw_200.csv", row.names=FALSE)
write.csv(sd3, "sd_tmle_200.csv", row.names=FALSE)
write.csv(sd4, "sd_sub_200.csv", row.names=FALSE)

write.csv(low1, "low_sdr_200.csv", row.names=FALSE)
write.csv(low2, "low_ipw_200.csv", row.names=FALSE)
write.csv(low3, "low_tmle_200.csv", row.names=FALSE)
write.csv(low4, "low_sub_200.csv", row.names=FALSE)

write.csv(high1, "high_sdr_200.csv", row.names=FALSE)
write.csv(high2, "high_ipw_200.csv", row.names=FALSE)
write.csv(high3, "high_tmle_200.csv", row.names=FALSE)
write.csv(high4, "high_sub_200.csv", row.names=FALSE)

