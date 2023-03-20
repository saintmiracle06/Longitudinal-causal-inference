## import packages
library(reshape2)
library(dplyr)
library(truncnorm)
library(data.table)

id<-10000000
#number of visits (K+1)
n_obs<-3

# generate the time-dependent variables h and l
# generate the time-dependent treatment z, which depends on d
# outcome y
h<-matrix(nrow=id,ncol=n_obs)
l<-matrix(nrow=id,ncol=n_obs)

z<-matrix(nrow=id,ncol=n_obs)
y<-matrix(nrow=id,ncol=n_obs)

expit<-function(x){exp(x)/(1+exp(x))}

# Baseline
#random error u
u<-rnorm(id,0,0.1)
## Baseline variable l_0
l_0 <- rbinom(id,1,expit(1.5+u))
#Time varying covariates l,h, get the time 1 value
l[,1]<-rbinom(id,1,expit(1+u))
h[,1]<-rbinom(id,1,expit(0.5+u))
## treatment
z[,1] <-rbinom(id,1,expit(0.8+u))
## Baseline observed y
y[,1]<-rnorm(id,1+u+l_0+z[,1]+h[,1]+l[,1],1)

# Time-varying covariates matrix
for(k in 2:n_obs){
  h[,k]=rbinom(id,1,expit(-1+0.5*h[,k-1]+0.1*(k-1)+u))
  l[,k]=rbinom(id,1,expit(-1+0.4*l[,k-1]+0.2*(k-1)+u))
  z[,k]=rbinom(id,1,expit(-1+0.5*h[,k]-0.2*l[,k]+z[,k-1]))
  y[,k]=rnorm(id,0.9*h[,k]-0.8*h[,k-1]-0.2*l[,k]+0.2*l[,k-1]+1.2*z[,k]-z[,k-1]+u+l_0,1)
}

## Wide to long transformation
h.dat=as.data.frame(h)
names(h.dat)=paste0("h.",0:2)

l.dat=as.data.frame(l)
names(l.dat)=paste0("l.",0:2)

z.dat=as.data.frame(z)
names(z.dat)=paste0("z.",0:2)

y.dat=as.data.frame(y)
names(y.dat)=paste0("y.",0:2)

dat<-data.frame(id=1:id,h.dat,z.dat,l.dat,y.dat,u,l_0)
head(dat)
## add intervention of d
d0 <-rbinom(dat$id,1,0.5*dat$z.0)
## d applied to z
dat_1<-dat
dat_1$d0<-d0
dat_1$y<-rnorm(id,1+u+l_0+dat_1$d0+dat_1$h.0+dat_1$l.0,1)
## Treatment effect is d in time 1
ATE1<-mean(dat_1$y)
## ATE1: 3.514948
## add intervention of d1
d1 <-rbinom(dat$id,1,0.5*dat_1$z.1)
dat_2<-dat_1
dat_2$d1<-d1
dat_2$y1<-rnorm(id,0.9*dat_2$h.1-0.8*dat_2$h.0-0.2*dat_2$l.1+0.2*dat_2$l.0+1.2*dat_2$d1-dat_2$d0+u+l_0,1)
## Treatment effect is d in time 1
ATE2<-mean(dat_2$y1)
## ATE2: 0.6406066
##Time 3
d2 <-rbinom(dat$id,1,0.5*dat_1$z.2)
dat_2<-dat_1
dat_2$d1<-d1
dat_2$y1<-rnorm(id,0.9*dat_2$h.1-0.8*dat_2$h.0-0.2*dat_2$l.1+0.2*dat_2$l.0+1.2*dat_2$d1-dat_2$d0+u+l_0,1)
## Treatment effect is d in time 1
ATE2<-mean(dat_2$y1)

