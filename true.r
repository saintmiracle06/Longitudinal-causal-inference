library(simcausal)
library(lmtp)
library(mc2d)
## PERT DISTRIBUTION
## https://search.r-project.org/CRAN/refmans/mc2d/html/pert.html
## https://en.wikipedia.org/wiki/PERT_distribution

id <- 100000
n_obs <- 3
D <- DAG.empty() +
  node("b1", distr="rcat.b1", probs = c(0.25,0.75)) +
  node(c("b2","b3"), distr = "mvtnorm::rmvnorm",
       asis.params = list(mean = "c(0,1)")) +
  node("h", t=0, distr="rnorm", mean=3.5, sd = 1) +
  node("l", t=0, distr="rnorm", mean=1.2, sd = 1) +
  node("z", t=0, distr="rpert", min=0, mode=plogis(0.3 + 0.2*b1 - 0.3*b2),max=1,shape=plogis(0.5 + 0.3*b3)) +
  node("Y", t=0, distr="rnorm", mean=b1 + 0.2 * b2 * b3) +
  
  node("l", t=1:n_obs, distr="rnorm", mean=plogis(0.5 - 0.3*b1 + 0.4*b2*b3 + 0.2*b1*b2*b3+ 0.6*l[t-1])) +
  node("h", t=1:n_obs, distr="rnorm", mean=plogis(1.5 - 0.5*b1 + 0.5*b1*b3 + 0.2*b1*b2*b3 + 1.5*h[t-1])) +
  
  node("z", t=1:n_obs, distr="rpert", min=0, mode=plogis(0.5 - 0.3*b1 + 0.2*h[t]-0.2*l[t] + 0.3*h[t]*l[t]- 0.2*b1*b2 + 0.2*b3),
       max=1,shape = plogis(0.5 + 0.3*b2 + 0.5*h[t]+0.4*l[t] + 0.3*h[t]*l[t] + 0.3*b1*b2 + 0.7*b2*b3)) +
  
  node("Y", t=1:n_obs, distr="rnorm", mean=b1 + 0.3*b2*b3+ 0.9*l[t]-0.65*l[t-1] + 0.5*h[t] +0.4*h[t-1] + 0.2*h[t-1]*l[t-1]
       + 0.3*z[t]*h[t] + 0.2*z[t-1]*l[t-1], EFU=TRUE)

D <- set.DAG(D)
dat <- sim(D,n=id)


g0 <- function(z0, b1, b2, b3) dpert(z0, min=0, mode=plogis(0.3 + 0.2*b1 - 0.3*b2),max=1,shape=plogis(0.5 + 0.3*b3))
g1 <- function(z1, b1, b2, b3, h1, l1) {
  dpert(z1,
        min=0, mode=plogis(0.5 - 0.3*b1 + 0.2*h1-0.2*l1 + 0.3*h1*l1- 0.2*b1*b2 + 0.2*b3),
        max=1,shape = plogis(0.5 + 0.3*b2 + 0.5*h1+0.4*l1 + 0.3*h1*l1 + 0.3*b1*b2 + 0.7*b2*b3))
}

gd0 <- function(z0,b1, b2, b3) g0(z0 + 0.1,b1, b2, b3) * as.numeric(z0 < 0.8) + g0(z0,b1, b2, b3) * as.numeric(z0 - 0.1 >= 0.8)
gd1 <- function(z1, b1, b2, b3, h1, l1) {
    g1(z1 - 0.1, b1, b2, b3, h1, l1) * as.numeric(z1 < 0.8) +
        g1(z1, b1, b2, b3, h1, l1) * as.numeric(z1 - 0.1 >= 0.8)
}

weights <- with(dat,
                cbind(gd0(z_0,b1, b2, b3) / g0(z_0,b1, b2, b3),
                gd1(z_1, b1, b2, b3, h_1, l_1) / g1(z_1, b1, b2, b3, h_1, l_1),
                gd1(z_2, b1, b2, b3, h_2, l_2) / g1(z_2, b1, b2, b3, h_2, l_2),
                gd1(z_3, b1, b2, b3, h_3, l_3) / g1(z_3, b1, b2, b3, h_3, l_3))
)

library(matrixStats)
we <- rowProds(weights)
true <- mean(dat$Y_3* we)

