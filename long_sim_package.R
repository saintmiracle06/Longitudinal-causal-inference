## library
library(reshape2)
library(dplyr)
library(SuperLearner)
library(truncnorm)
library(data.table)
library(simcausal)


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

ATE1<-rep(NA,500)
ATE2<-rep(NA,500)
ATE3<-rep(NA,500)

for (i in 1:500){
  id<-1000
  #number of visits (K+1)
  n_obs<-2
  
  D <- DAG.empty() + 
    node("b", distr="rcat.b1", probs = c(0.25,0.75)) +
    node("l", t=0, distr="rcat.b1", probs = c(0.20,0.80)) +
    node("h", t=0, distr="rcat.b1", probs = c(0.81,0.19)) + 
    node("z", t=0, distr="rbern", prob=plogis(1 - 0.3*b + 0.5*h[t] -0.3*l[t])) +
    
    node("l", t=1:n_obs, distr="rbern", prob=plogis(1 - 0.3*b + 1.5*l[t-1])) +
    node("h", t=1:n_obs, distr="rbern", prob=plogis(1.5 - 0.5*b + 0.5*h[t-1])) +
    node("z", t=1:n_obs, distr="rbern", prob=plogis(0.5 - 0.3*b + 0.5*h[t]-0.4*l[t])) +
    node("d", t=0:n_obs, distr="rbern", prob=plogis(-1+0.5*z[t])) +
    node("Y", t=0:n_obs, distr="rnorm", mean=b + 0.9*l[t] + 0.5*h[t] + 0.2*z[t], EFU=TRUE)
  
   D <- set.DAG(D)
  
   dat <- sim(D,n=id)
  
   dat_long<-reshape(data = dat,sep = "_",,direction="long",idvar="ID",varying = -c(1:2))
   
   dat_long$r<-rbinom(id,1,0.85)
   dat_long$y_0<-ifelse(dat_long$r==1, dat_long$Y, NA)
  
   data_lag <- dat_long %>%                            # Add lagged column
     group_by(ID) %>%
     dplyr::mutate(h_lag1 = dplyr::lag(h, n = 1, default =first(h)),
                   h_lag2 = dplyr::lag(h_lag1, n = 1, default =first(h_lag1)),
                   l_lag1 = dplyr::lag(l, n = 1, default =first(l)),
                   l_lag2 = dplyr::lag(l_lag1, n = 1, default =first(l_lag1)),
                   z_lag1 = dplyr::lag(z, n = 1, default =first(z)),
                   z_lag2 = dplyr::lag(z_lag1, n = 1, default =first(z_lag1)))%>%
     as.data.frame()
   
   my_vars = c("time","z","h","l","z_lag1", "z_lag2", "h_lag1", "h_lag2", "l_lag1", "l_lag2","b")
   y = data_lag[data_lag$r ==1,]$y_0
   X = subset(data_lag[data_lag$r ==1,], select = my_vars)
   ## fit specified model
   model1<-glmnet3(X, y, family ="gaussian", id = data_lag[data_lag$r ==1,]$ID)
   ## predict on specified dataset where d is specified by user
   data_lag_1<-data_lag
   data_lag_1$z <-data_lag$d
   ## Treatment effect is d 
   data_lag_1$y_1 <- predict.glmnet3(model1, data_lag_1)
   ATE1[i]<-mean(data_lag_1[data_lag_1$time==0,]$y_1)
   
   ## First treatment
   dat1<-data_lag_1[which(data_lag_1$time!=0),]
   my_vars2 = c("time","z","z_lag1", "z_lag2", "h_lag1", "h_lag2", "l_lag1", "l_lag2","b")
   y = dat1$y_1
   X = subset(dat1, select = my_vars2)
   model2<-glmnet3(X, y, family ="gaussian", id = dat1$ID)
   dat2<-dat1
   dat2$z <-dat1$d
   dat2$y_2 <- predict.glmnet3(model2, dat2)
   ATE2[i]<-mean(dat2$y_2)
   
   ## Second treatment
   dat3<-dat2[which(dat2$time!=1),]
   my_vars3 = c("time","z","z_lag2","h_lag2","l_lag2","b")
   y = dat3$y_1
   X = subset(dat3, select = my_vars3)
   model3<-glmnet3(X, y, family ="gaussian", id = dat3$ID)
   dat4<-dat3
   dat4$z <-dat3$d
   dat4$y_3 <- predict.glmnet3(model3, dat4)
   ATE3[i]<-mean(dat4$y_3)
   print(i)
}

png(filename="ATE1.png")
hist(ATE1, xlim=c(3.8, 4.1))
abline(v=4.019579, col='red', lwd=3, lty='dashed')
dev.off()

png(filename="ATE2.png")
hist(ATE2, xlim=c(2.9, 3.2))
abline(v=3.052384, col='red', lwd=3, lty='dashed')
dev.off()

png(filename="ATE3.png")
hist(ATE3, xlim=c(2.9, 3.2))
abline(v=2.951065, col='red', lwd=3, lty='dashed')
dev.off()