## checking sl3 library
## function loop through time
## glue for calling the string
## test code in simulations: discrete data, superlearner with interactions(saturated model)
## past info
## generate more variables L1
## true value : 10 million
## generate saturated model

## how to check truth 10 million
## time 1
## don't mask y
## H, Z, D, d(z)~I(\epsilon)
## y~ bernolli(q(z, H))
## y(d)~bernolli(q(d(z),H))

## then do DGP of 500 samples and check the distribution

## time 2
## H1=f(U)
## z1~Ber(P)


## library
library(reshape2)
library(dplyr)
library(SuperLearner)
library(truncnorm)
library(data.table)


## simulation study:
## generate longitudinal data
#sample size
## bias 
## variance
## screen regression function
glmnet3 <- function(X, y, family = c("gaussian", "binomial"), id = NULL) {
  if (!is.null(id)) {
    # need to match index with fold number for cv.glmnet
    folds <- origami::make_folds(
      nrow(X), fold_fun = origami::folds_vfold,
      cluster_ids = id, V = 10
    )
    foldid <- vector("numeric", nrow(X))
    for (i in 1:nrow(X)) {
      for (v in 1:10) {
        if (i %in% folds[[v]]$validation_set) {
          foldid[i] <- v
          break
        }
      }
    }
  } else {
    foldid <- NULL
  }
  ans <- list(covars = names(X))
  if (ncol(X) == 1) {
    x <- as.matrix(X)
    ans$fit <- glm(y ~ ., data = cbind(y = y, X), family = match.arg(family))
  } else {
    f <- as.formula(paste0("~ .^", ncol(X)))
    x <- model.matrix(f, X)[, -1]
    ans$fit <- glmnet::cv.glmnet(x, y, family = match.arg(family), foldid = foldid)
  }
  # x <- try(model.matrix(f, X)[, -1])
  # if (inherits(x, "try-error")) browser()
  # ans$fit <- glmnet::cv.glmnet(x, y, family = match.arg(family), foldid = foldid)
  ans
}
predict.glmnet3 <- function(object, newx) {
  if (inherits(object$fit, "glm")) {
    return(predict(object$fit, newx, type = "response"))
  }
  X <- newx[, object$covars, drop = TRUE]
  f <- as.formula(paste0("~ .^", ncol(X)))
  X <- model.matrix(f, X)[, -1]
  as.vector(predict(object$fit, X, type = "response")[, 1])
}

ATE1<-rep(NA,200)
ATE2<-rep(NA,200)

## MC simulation
for (i in 1:200){
id<-3000
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
## Baselin variable l_0
l_0 <- rbinom(id,1,expit(1.5+u))
#h[,1]<-rnorm(id,0+u,1)
l[,1]<-rbinom(id,1,expit(1+u))
h[,1]<-rbinom(id,1,expit(0.5+u))

## function D for z  
z[,1] <-rbinom(id,1,expit(0.8+u))
y[,1]<-rnorm(id,1+u+l_0+z[,1]+h[,1]+l[,1],1)

# Time-varying variables
for(k in 2:n_obs){
  #h[,k]=rnorm(id,0.8*h[,k-1]-z[,k-1]+0.1*(k-1)+u,1)
  h[,k]=rbinom(id,1,expit(-1+0.5*h[,k-1]+0.1*(k-1)+u))
  l[,k]=rbinom(id,1,expit(-1+0.4*l[,k-1]+0.2*(k-1)+u))
  z[,k]=rbinom(id,1,expit(-1+0.5*h[,k]+z[,k-1]))
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

## generate lagged variables based on carry-over the first observations
## zlag1.dat=as.data.frame(cbind(rep(0,id),z.dat[,1:2]))
zlag1.dat=as.data.frame(cbind(z.dat[,1],z.dat[,1:2]))
zlag2.dat=as.data.frame(cbind(z.dat[,1:2],z.dat[,1]))
names(zlag1.dat)=paste0("zlag1.",0:2)
names(zlag2.dat)=paste0("zlag2.",0:2)

hlag1.dat=as.data.frame(cbind(h.dat[,1],h.dat[,1:2]))
hlag2.dat=as.data.frame(cbind(h.dat[,1:2],h.dat[,1]))
names(hlag1.dat)=paste0("hlag1.",0:2)
names(hlag2.dat)=paste0("hlag2.",0:2)

llag1.dat=as.data.frame(cbind(rep(0,id),l.dat[,1:2]))
llag2.dat=as.data.frame(cbind(rep(0,id),rep(0,id),l.dat[,1]))
names(llag1.dat)=paste0("llag1.",0:2)
names(llag2.dat)=paste0("llag2.",0:2)

dat<-data.frame(id=1:id,h.dat,z.dat,l.dat, zlag1.dat,zlag2.dat,hlag1.dat,hlag2.dat,llag1.dat,llag2.dat,y.dat,u,l_0)
#dat<-data.frame(id=1:id,h.dat,z.dat,y.dat,u)
#dat<-data.frame(id=1:id,h.dat,z.dat,zlag1.dat,zlag2.dat,hlag1.dat,hlag2.dat,y.dat,u)

#dat.long<-reshape(data = dat,varying=c(paste0("z.",0:2),paste0("h.",0:2),
#                                       paste0("y.",0:2)),direction="long",idvar="id")
dat.long<-reshape(data = dat,varying=c(paste0("h.",0:2),paste0("z.",0:2),paste0("l.",0:2), paste0("zlag1.",0:2),
                                       paste0("zlag2.",0:2),paste0("hlag1.",0:2),paste0("hlag2.",0:2),
                                       paste0("llag1.",0:2),paste0("llag2.",0:2),
                                       paste0("y.",0:2)),direction="long",idvar="id")

dat.long<-dat.long[order(dat.long$id,dat.long$time),]
dat.long$r<-rbinom(id,1,0.85)
dat.long$y_0<-ifelse(dat.long$r==1, dat.long$y, NA)

## Regression ##
## add time in the regression and carried over (Add interaction term)
## ATE estimators
##colnames(dat.long)

my_vars = c("time","z","h","l","zlag1", "zlag2", "hlag1", "hlag2", "llag1", "llag2","u","l_0")
y = dat.long[dat.long$r ==1,]$y_0
X = subset(dat.long[dat.long$r ==1,], select = my_vars)
## fit specified model
model1<-glmnet3(X, y, family ="gaussian", id = dat.long[dat.long$r ==1,]$id)
## predict on specified dataset where d is specified by user
## d function z
d <-rbinom(id,1,0.5*dat.long$z)
## d applied to z
dat.long_1<-dat.long
dat.long_1$z <-d
## Treatment effect is d 
dat.long_1$y_1 <- predict.glmnet3(model1, dat.long_1)
ATE1[i]<-mean(dat.long_1[dat.long_1$time==0,]$y_1)
## -0.0440883##
# reg_0 <- lm(y_0 ~ z + h + l + zlag1 + zlag2 + hlag1 + hlag2 + llag1 + llag2 + time +
#                   z*h + z*l + z*time + z*hlag1 + z*hlag2 + z*llag1 + z*llag2 +
#                   h*time + h*zlag1 + h*zlag2 + h*l + h*llag1 + h*llag2 +
#                   zlag1*time + hlag1*time + llag1*time + llag1*time, data = dat.long)
# summary(reg_0)
# dat.long_1$y_reg <- predict(reg_0, newdata = dat.long_1) ## warning message?#

# ATE1[i]<- mean(dat.long[which(dat.long$z==1),]$y_1)-mean(dat.long[which(dat.long$z==0),]$y_1)
# ATE1<- mean(dat.long_1[which(dat.long_1$z==1),]$y_1)-mean(dat.long_1[which(dat.long_1$z==0),]$y_1)
#ATE1 = 1.000109
##make the z observed z
dat.long_1$z<-dat.long$z
dat1<-dat.long_1[which(dat.long_1$time>0),]

## limited to z_(s-k), h_(s-k)
my_vars2 = c("time","zlag1", "zlag2", "hlag1", "hlag2", "llag1", "llag2","u","l_0")
y = dat1$y_1
X = subset(dat1, select = my_vars2)
model2<-glmnet3(X, y, family ="gaussian", id = dat1$id)

## just including the time-varying variables before this time step (make clear)
# reg_1 <- lm(y_1 ~ l + zlag1 + hlag1 + llag1 + llag2 + time + 
#               z*h + z*l + z*time + z*hlag1 + z*hlag2 + z*llag1 + z*llag2 +
#               h*time + h*zlag1 + h*zlag2 + h*l + h*llag1 + h*llag2 + 
#               zlag1*time + hlag1*time + llag1*time + llag1*time, data = dat1)
# summary(reg_1)
## d function z
## d applied to z
d <-rbinom(id,1,0.5*dat1$z)
## d applied to z
dat2<-dat1
dat2$z <-d
dat2$y_2 <- predict.glmnet3(model2, dat2)
# dat1$y_2 <- predict(reg_1, newdata = dat1)
#ATE2[i]<- mean(dat1[which(dat1$z==1),]$y_2)-mean(dat1[which(dat1$z==0),]$y_2)
ATE2[i]<- mean(dat2[dat2$time==1,]$y_2)
}
hist(ATE1, xlim=c(2.5, 4))
abline(v=3.514948, col='red', lwd=3, lty='dashed')

hist(ATE2, xlim=c(0, 2))
abline(v=0.6406066, col='red', lwd=3, lty='dashed')

## automatically generate lag function k is number of lags
k=2
for (i in 1:k){
  ## Use regression by can use any other functions
  ## k regression contain all data in the history
  reg<-lm(DT[[paste0("y",sep="_",i-1)]] ~ DT[[paste0("z",sep="_",i)]] + DT[[paste0("h",sep="_",i)]])
  summary(reg)
  prediction <- predict(reg, newdata = DT)
  #df[,paste0("MktCap",j,sep="")]<-df$Shares*df[,paste0("Price",j,sep="")]
  DT[,paste0("y",sep="_",i)] <- prediction
}

