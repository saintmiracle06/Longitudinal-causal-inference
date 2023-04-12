## import packages
library(reshape2)
library(dplyr)
library(truncnorm)
library(data.table)
library(simcausal)

id<-100000
#number of visits (K+1)
n_obs<-3

# generate the time-dependent variables h and l
# generate the time-dependent treatment z, which depends on d
# outcome y
h<-matrix(nrow=id,ncol=n_obs)
l<-matrix(nrow=id,ncol=n_obs)

z<-matrix(nrow=id,ncol=n_obs)
d<-matrix(nrow=id,ncol=n_obs)

y<-matrix(nrow=id,ncol=n_obs)
y_inv<-matrix(nrow=id,ncol=n_obs)

expit<-function(x){exp(x)/(1+exp(x))}

# Baseline
#random error u
set.seed(123456)
u1<-rnorm(id,0,0.1)
u2<-rnorm(id,0,0.1)
u3<-rnorm(id,0,0.1)
u4<-rnorm(id,0,0.1)
u5<-rnorm(id,0,0.1)

## Baseline variable l_0
l_0 <- rbinom(id,1,expit(1.5+u1))
#Time varying covariates l,h, get the time 1 value
l[,1]<-rbinom(id,1,expit(1+u2))
h[,1]<-rbinom(id,1,expit(0.5+u3))
## treatment
z[,1] <-rbinom(id,1,expit(0.8+u4))
d[,1] <-rbinom(id,1,0.5*z[,1])

## Baseline observed y
y[,1]<-rnorm(id,1+u5+l_0+z[,1]+h[,1]+l[,1],1)
y_inv[,1]<-rnorm(id,1+u5+l_0+d[,1]+h[,1]+l[,1],1)


# Time-varying covariates matrix
for(k in 2:n_obs){
  h[,k]=rbinom(id,1,expit(-1+0.5*h[,k-1]+0.1*(k-1)+rnorm(1,0,0.1)))
  l[,k]=rbinom(id,1,expit(-1+0.4*l[,k-1]+0.2*(k-1)+rnorm(1,0,0.1)))
  z[,k]=rbinom(id,1,expit(-1+0.5*h[,k]-0.2*l[,k]+z[,k-1]))
  y[,k]=rnorm(id,0.9*h[,k]-0.8*h[,k-1]-0.2*l[,k]+0.2*l[,k-1]+1.2*z[,k]-z[,k-1]+rnorm(1,0,0.1)+l_0,1)
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

y.dat=as.data.frame(y)
names(y.dat)=paste0("y.",0:2)

dat<-data.frame(id=1:id,h.dat,z.dat,l.dat,y.dat,l_0)
head(dat)
## add intervention of d
## d applied to z
dat$d0<-d[,1]
dat$y_inv<-y_inv[,1]
## Treatment effect is d in time 1
ATE1<-mean(dat$y_inv)
## ATE1: 3.516358
## add intervention of d in time point 2
z_inv1<-rbinom(id,1,expit(-1+0.5*h[,2]-0.2*l[,2]+d[,1]))
d1<-rbinom(dat$id,1,0.5*z_inv1)
dat$d1<-d1
dat$y_inv1<-rnorm(id,0.9*dat$h.1-0.8*dat$h.0-0.2*dat$l.1+0.2*dat$l.0+1.2*dat$d1-dat$d0+rnorm(1,0,0.1)+l_0,1)
## Treatment effect is d in time 1
ATE2<-mean(dat$y_inv1)
## ATE2: 0.5313002
##Time 3
z_inv2<-rbinom(id,1,expit(-1+0.5*h[,2]-0.2*l[,2]+d1))
d2 <-rbinom(dat$id,1,0.5*z_inv2)
dat$d2<-d2
dat$y_inv2<-rnorm(id,0.9*dat$h.1-0.8*dat$h.0-0.2*dat$l.1+0.2*dat$l.0+1.2*dat$d2-dat$d1+rnorm(1,0,0.1)+l_0,1)
## Treatment effect is d in time 1
ATE3<-mean(dat$y_inv2)
## ATE3: 0.9488658

