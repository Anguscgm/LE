# library(randnet)
library(PRROC)

library(igraph)



library(pROC)
library(latentnet)
library(network)
library(pcaPP)
library(doParallel)
library(RSpectra)
library(rmatio)




net.gen.from.P <- function(P,mode="undirected"){
  n <- nrow(P)
  if(mode=="undirected"){
    upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1


    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)

  }else{
    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    diag(A) <- 0
  }
  return(A)
}



ego.sample <- function(P,rho,deg){
  P <- P*deg/mean(rowSums(P))
  A <- net.gen.from.P(P)
  n <- nrow(P)
  sample.n <- ceiling(n*rho)
  sample.index <- sample(1:n,size=sample.n,replace=FALSE)
  out.index <- (1:n)[-sample.index]
  out.n <- length(out.index)
  new.index <- c(sample.index,out.index)
  A <- A[new.index,new.index]
  P <- P[new.index,new.index]
  A.obs <- A
  A.obs[-(1:sample.n),-(1:sample.n)] <- NA
  A.out <- A[-(1:sample.n),-(1:sample.n)]
  P.out <- P[-(1:sample.n),-(1:sample.n)]
  result <- list(sample.n=sample.n,out.n=out.n,A=A,P=P,A.obs=A.obs,A.out=A.out,P.out=P.out)
  return(result)
}


real.sample <- function(A,rho){
  # P <- P*deg/mean(rowSums(P))
  # A <- net.gen.from.P(P)
  n <- nrow(A)
  sample.n <- ceiling(n*rho)
  sample.index <- sample(1:n,size=sample.n,replace=FALSE)
  out.index <- (1:n)[-sample.index]
  out.n <- length(out.index)
  new.index <- c(sample.index,out.index)
  A <- A[new.index,new.index]
  # P <- P[new.index,new.index]
  A.obs <- A
  A.obs[-(1:sample.n),-(1:sample.n)] <- NA
  A.out <- A[-(1:sample.n),-(1:sample.n)]
  # P.out <- P[-(1:sample.n),-(1:sample.n)]
  # result <- list(sample.n=sample.n,out.n=out.n,A=A,P=P,A.obs=A.obs,A.out=A.out,P.out=P.out)
  result <- list(sample.n=sample.n,out.n=out.n,A=A,A.obs=A.obs,A.out=A.out)
  return(result)
}

nonuniformbyA.real.sample <- function(A,rho,n.ratio = c(.33,.34,.33),
                                 sample.ratio= c(1.5,1,.5)){
  if (length(n.ratio)!=length(sample.ratio)){
    stop("n.ratio and sample.ratio must have the same legnth.")
  }
  # P <- P*deg/mean(rowSums(P))
  # A <- net.gen.from.P(P)
  A.index <- order(rowSums(A),decreasing = T)
  # P <- P[A.index,A.index]
  A <- A[A.index,A.index]
  n <- round(nrow(A)*n.ratio)
  # Add/ subtract from middle layer if n does not match after rounding.
  n[2] = n[2] + nrow(A) - sum(n)
  sample.n <- pmin(ceiling(n*rho*sample.ratio),n)
  sample.index <- c()
  temp.ind <- 0
  for (i in 1:length(n.ratio)){
    sample.index <- c(sample.index,
                      temp.ind + sample(1:n[i],size=sample.n[i],replace=FALSE))
    temp.ind <- temp.ind + n[i]
  }
  n <- sum(n)
  sample.n <- sum(sample.n)
  out.index <- (1:n)[-sample.index]
  out.n <- length(out.index)
  new.index <- c(sample.index,out.index)
  A <- A[new.index,new.index]
  # P <- P[new.index,new.index]
  A.obs <- A
  A.obs[-(1:sample.n),-(1:sample.n)] <- NA
  A.out <- A[-(1:sample.n),-(1:sample.n)]
  # P.out <- P[-(1:sample.n),-(1:sample.n)]
  result <- list(sample.n=sample.n,out.n=out.n,A=A,A.obs=A.obs,A.out=A.out)
  return(result)
}

Global.CC <- function(A){
  diag(A) <- 0
  A2 <- A%*%A
  A3 <- A2%*%A
  return(sum(diag(A3))/sum(A2))
}


ERGM.link.pred <- function(A){
  na.index <- which(rowSums(is.na(A))>0)
  deg.norm <- mean(rowSums(A[-na.index,]))
  deg.NA <- mean(rowSums(A[na.index,],na.rm=T))
  mid.deg <- ceiling(0.5*deg.NA+0.5*deg.norm)
  data.net <- as.network(A,directed=FALSE)

  n <- nrow(A)

  while(TRUE){
    test1 <- ergm(data.net~edges+triangle,estimate="MPLE" )
    if(sum(is.na(test1$coef))==0){
      break
    }
  }
  summary(test1)

  sim.data <- simulate(test1,nsim=100)
  pred.net <- matrix(0,n,n)
  for(T in 1:100){
    pred.net <- pred.net + as.matrix(sim.data[[T]])
  }
  pred.net <- pred.net/100
  return(pred.net)
}




ERGM.gwesp.pred <- function(A){
  na.index <- which(rowSums(is.na(A))>0)
  deg.norm <- mean(rowSums(A[-na.index,]))
  deg.NA <- mean(rowSums(A[na.index,],na.rm=T))
  mid.deg <- ceiling(0.5*deg.NA+0.5*deg.norm)
  data.net <- as.network(A,directed=FALSE)

  n <- nrow(A)

  while(TRUE){
                                        #test1 <- ergm(data.net~edges + gwesp(0.5,fixed=T),control=control.ergm(MCMLE.maxit = 1))
      test1 <- ergm(data.net~edges + gwesp(0.5,fixed=T),estimate="MPLE")
    if(sum(is.na(test1$coef))==0){
      break
    }
  }
  summary(test1)

  sim.data <- simulate(test1,nsim=100)
  pred.net <- matrix(0,n,n)
  for(T in 1:100){
    pred.net <- pred.net + as.matrix(sim.data[[T]])
  }
  pred.net <- pred.net/100
  return(pred.net)
}



library(irlba)
library(MASS)

Block.egocentric.fit <- function(A,train.index,K){
  n <- nrow(A)
  holdout.index <- setdiff(1:n,train.index)
  A.new <- A
  n.train <- length(train.index)
  n.holdout <- n - n.train
  A1 <- A.new[train.index,]

  if(K==1){
    A0 <- A.new[train.index,train.index]
    pb <- sum(A0)/n.train^2
    if(pb < 1e-6) pb <- 1e-6
    if(pb > 1- 1e-6) pb <- 1-1e-6
    Phat <- matrix(pb,n,n)
    return(Phat)
  }

  A1.svd <- irlba(A1,nv=K,nu=K)
  V <- A1.svd$v
  km <- kmeans(V,centers=K,nstart=30,iter.max=30)

  degrees <- colSums(A1)
  no.edge <- sum(degrees==0)


  B <- matrix(0,K,K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      N.1i <- intersect(train.index,which(km$cluster==i))
      N.2i <- intersect(holdout.index,which(km$cluster==i))
      N.1j <- intersect(train.index,which(km$cluster==j))
      N.2j <- intersect(holdout.index,which(km$cluster==j))
      B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(length(N.1i)*length(N.1j)+length(N.1j)*length(N.2i)+length(N.1i)*length(N.2j)+1)
      #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
    }
  }
  B <- B+t(B)
  Theta <- matrix(0,n,K)
  for(i in 1:K){
    N.1i <- intersect(train.index,which(km$cluster==i))
    N.2i <- intersect(holdout.index,which(km$cluster==i))
    B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(length(N.1i)*(length(N.1i)-1)/2+length(N.1i)*length(N.2i)+1)
    Theta[which(km$cluster==i),i] <- 1

  }

  Phat <- Theta %*%B %*% t(Theta)
  Phat[Phat > (1-1e-5)] <- 1-1e-5
  Phat[Phat < 1e-5] <- 1e-5

  return(Phat)

}



Block.egocentric.validation.auc <- function(A,train.index,test.index,K){
  Phat <- Block.egocentric.fit(A,train.index,K=K)
  predictors <- Phat[test.index,-train.index]
  labels <- A[test.index,-train.index]
  # test.index <- which(upper.tri(labels))
  # predictors <- predictors[test.index]
  # labels <- labels[test.index]
  predictors <- as.numeric(predictors)
  labels <- as.numeric(labels)
  if(length(unique(labels))>1){
      auc.value <- as.numeric(auc(response=labels,predictor=predictors))
  }else{
      auc.value <- 0.5
  }
  return(auc.value)

}

Block.egocentric.CV <- function(A,train.index,Kmax=5,B=3,rho=0.3){
  n <- nrow(A)
  n.train <- length(train.index)
  Kmax = min(Kmax,n.train-ceiling(rho*n.train)-1)
  cv.mat <- matrix(NA,B,Kmax-1)
  for(bb in 1:B){
    validation.index <- sample(train.index,size=ceiling(rho*n.train))
    sub.train.index <- setdiff(train.index,validation.index)

    for(k in 2:Kmax){
      cv.mat[bb,k-1] <- Block.egocentric.validation.auc(A,sub.train.index,validation.index,k)
    }
  }
  return(which.max(colMeans(cv.mat))+1)

}


library(VGAM)
library(amen)

ego.perf <- function(A.impute,A.true,A.impute.out,A.true.out,PR=TRUE){
    #print(length(unique(as.numeric(A.true.out))))
    if((sum(A.impute.out)==0)||(length(unique(as.numeric(A.true.out)))<2)){
        tmp.auc=0.5
        tmp.fpr=NA
        tmp.tpr=NA
        tmp.pr.auc=0
        tmp.precision=NA
        tmp.recall=NA
        # tmp.norm=NA
        # tmp.rel.err=NA
        # tmp.rel.err.full=NA

      #   dhat <-  mean(A.impute)
      #   derr <- abs(dhat-mean(A.true))/mean(A.true)
      # eig.P <- eigs_sym(A.true,k=1)
      # eig.Phat <- eigs_sym(A.impute,k=1)
      # ec.P <- sum(max(abs(eig.P$vectors[,1])) - abs(eig.P$vectors[,1]))
      # ec.Phat <- sum(max(abs(eig.Phat$vectors[,1])) - abs(eig.Phat$vectors[,1]))
      # cerr <- abs(ec.Phat-ec.P)/ec.P
      # cc.P <- Global.CC(A.true)
      # cc.Phat <- Global.CC(A.impute)
      # cc.err <- abs(cc.Phat-cc.P)/cc.P

    }else{
        
          tmp.ROC <- PRROC::pr.curve(scores.class0=as.numeric(A.impute.out),weights.class0=as.numeric(A.true.out),curve=TRUE)
          tmp.precision <- tmp.ROC$curve[,1]#1-tmp.ROC$specificities
          tmp.recall <- tmp.ROC$curve[,2]#tmp.ROC$sensitivities #### This is wrong, we need to swap.
          tmp.pr.auc <- as.numeric(tmp.ROC$auc.integral)
          
        
          tmp.ROC <- pROC::roc(as.numeric(A.true.out)~as.numeric(A.impute.out))
          tmp.fpr <- 1-tmp.ROC$specificities
          tmp.tpr <- tmp.ROC$sensitivities
          tmp.auc <- as.numeric(tmp.ROC$auc)
          # tmp.norm <- norm(A.true.out-A.impute.out,"F")
          # tmp.rel.err <- norm(A.true.out-A.impute.out,"F")/norm(A.true.out,"F")
          # tmp.rel.err.full <- norm(A.true-A.impute,"F")/norm(A.true,"F")
        
      #   dhat <- mean(A.impute)
      #   derr <- abs(dhat-mean(A.true))/mean(A.true)
      # eig.P <- eigs_sym(A.true,k=1)
      # eig.Phat <- eigs_sym(A.impute,k=1)
      # ec.P <- sum(max(abs(eig.P$vectors[,1])) - abs(eig.P$vectors[,1]))
      # ec.Phat <- sum(max(abs(eig.Phat$vectors[,1])) - abs(eig.Phat$vectors[,1]))
      # cerr <- abs(ec.Phat-ec.P)/ec.P
      # cc.P <- Global.CC(A.true)
      # cc.Phat <- Global.CC(A.impute)
      # cc.err <- abs(cc.Phat-cc.P)/cc.P
    }
    # return(list(fpr=tmp.fpr,tpr=tmp.tpr,auc=tmp.auc,derr=derr,cerr=cerr,precision=tmp.precision,recall=tmp.recall,pr.auc=tmp.pr.auc,cc=cc.err,norm=tmp.norm,rel.err=tmp.rel.err))
  return(list(fpr=tmp.fpr,tpr=tmp.tpr,auc=tmp.auc,precision=tmp.precision,
              recall=tmp.recall,pr.auc=tmp.pr.auc))
}







ego.fit.DCBM <- function(A,train.index,K){
  n <- nrow(A)
  holdout.index <- setdiff(1:n,train.index)
  A.new <- A[c(train.index,holdout.index),c(train.index,holdout.index)]
  n.holdout <- length(holdout.index)
  n.train <- n-n.holdout
  A1 <- A.new[1:n.train,]
  A1.svd <- irlba(A1,nu=K,nv=K)
  #V <- A1.svd$v[,1:K]
  V <- A1.svd$v
  if(K==1) {V.norms <- abs(V)}else{
    V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
  }
  iso.index <- which(V.norms==0)
  Psi <- V.norms
  Psi <- Psi / max(V.norms)
  Psi.outer <- outer(Psi,Psi)
  inv.V.norms <- 1/V.norms # Only up to a scale anyway.
  inv.V.norms[iso.index] <- 1
  V.normalized <- diag(inv.V.norms)%*%V
  
  if(K==1){
    N.1i <- 1:n.train
    N.2i <- (n.train+1):n
    # pb <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)
    pb <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer)[N.1i])/2+1)
    
    P.hat.holdout <-  diag(Psi[(n.train+1):n])%*%matrix(1,ncol=(n-n.train),nrow=(n-n.train))%*%diag(Psi[(n.train+1):n])*pb
    P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
    P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
    A.2 <- A.new[(n.train+1):n,(n.train+1):n]
    sum.index <- lower.tri(A.2)
    loglike <- -sum(A.2[sum.index]*log(P.hat.holdout[sum.index])) - sum((1-A.2[sum.index])*log(1-P.hat.holdout[sum.index]))
    
    l2 <- sum((A.2[sum.index]-P.hat.holdout[sum.index])^2)
    return(list(loglike=loglike,l2=l2))
  }
  
  
  km <- kmeans(V.normalized,centers=K,nstart=30,iter.max=30)
  
  degrees <- colSums(A1)
  no.edge <- sum(degrees==0)
  
  B <- matrix(0,K,K)
  
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      N.1i <- intersect(1:n.train,which(km$cluster==i))
      N.2i <- intersect((n.train+1):n,which(km$cluster==i))
      N.1j <- intersect(1:n.train,which(km$cluster==j))
      N.2j <- intersect((n.train+1):n,which(km$cluster==j))
      B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j]) + sum(A.new[N.1j,N.2i])+1)/(sum(Psi.outer[N.1i,N.1j]) + sum(Psi.outer[N.1i,N.2j]) + sum(Psi.outer[N.1j,N.2i])+1)
      #B[i,j] <- (sum(A.new[N.1i,N.1j]) + sum(A.new[N.1i,N.2j])+1)/(length(N.1i)*length(N.1j)+length(N.1i)*length(N.2j)+1)
    }
  }
  B <- B+t(B)
  Theta <- matrix(0,n,K)
  for(i in 1:K){
    N.1i <- intersect(1:n.train,which(km$cluster==i))
    N.2i <- intersect((n.train+1):n,which(km$cluster==i))
    # B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer))+1)
    B[i,i] <- (sum(A.new[N.1i,N.1i])/2 + sum(A.new[N.1i,N.2i])+1)/(sum(Psi.outer[N.1i,N.1i])/2 + sum(Psi.outer[N.1i,N.2i]) - sum(diag(Psi.outer)[N.1i])/2+1)
    Theta[which(km$cluster==i),i] <- 1
    
  }
  tmp.imt.mat <- Theta*Psi
  P.hat.holdout <-  tmp.imt.mat%*%B%*%t(tmp.imt.mat)
  P.hat.holdout[P.hat.holdout<1e-6] <- 1e-6
  P.hat.holdout[P.hat.holdout>(1-1e-6)] <- 1-1e-6
  return(P.hat.holdout)
  
}




DCBM.egocentric.validation.auc <- function(A,train.index,test.index,K){
  Phat <- ego.fit.DCBM(A,train.index,K=K)
  # predictors <- Phat[test.index,test.index]
  # labels <- A[test.index,test.index]
  # test.index <- which(upper.tri(labels))
  # predictors <- predictors[test.index]
  # labels <- labels[test.index]
  predictors <- as.numeric(Phat[test.index,-train.index])
  labels <- as.numeric(A[test.index,-train.index])
  if(length(unique(labels))>1){
    auc.value <- as.numeric(auc(response=labels,predictor=predictors))
  }else{
    auc.value <- 0.5
  }
  return(auc.value)
  
}

DCBM.egocentric.CV <- function(A,train.index,Kmax=5,B=3,rho=0.3){
  n <- nrow(A)
  n.train <- length(train.index)
  Kmax = min(Kmax,n.train-ceiling(rho*n.train)-1)
  cv.mat <- matrix(NA,B,Kmax-1)
  for(bb in 1:B){
    validation.index <- sample(train.index,size=ceiling(rho*n.train))
    sub.train.index <- setdiff(train.index,validation.index)
    
    for(k in 2:Kmax){
      cv.mat[bb,k-1] <- DCBM.egocentric.validation.auc(A,sub.train.index,validation.index,k)
    }
  }
  return(which.max(colMeans(cv.mat))+1)
  
}


library(MASS)

#### general low-rank estimation
#### A is the adjacency matrix
#### train.index: the index of the observed nodes
#### r is the rank
SE.egocentric <- function(A,train.index,r=3){
  n <- nrow(A)
  n.train <- length(train.index)
  test.index <- setdiff(1:n,train.index)
  # A11 <- A[train.index,train.index]
  # A12 <- A[train.index,test.index]
  # A21 <- A[test.index,train.index]
  A1 <- A[train.index,]
  SVD1 <- irlba(A1,nv=r,nu=r)
  if(r>1){
      Ahat1 <- SVD1$u[,1:r]%*%(t(SVD1$v[,1:r])*SVD1$d[1:r])
  }else{
      Ahat1 <- matrix(SVD1$u[,1:r],ncol=1)%*%(t(matrix(SVD1$v[,1:r],ncol=1))*SVD1$d[1:r])
  }
  Ahat11 <- Ahat1[,train.index]
  Ahat12 <- Ahat1[,test.index]
  Ahat21 <- t(Ahat12)
  
  
  A11.est.inv <- ginv(Ahat11)
  A11.est.inv <- (A11.est.inv+t(A11.est.inv))/2
  A22.pred <- Ahat21%*%A11.est.inv%*%Ahat12
  A.new <- A
  A.new[train.index,train.index] <- Ahat11
  A.new[train.index,test.index] <- Ahat12
  A.new[test.index,train.index] <- Ahat21
  A.new[test.index,test.index] <- A22.pred
  return(A.new)

}


##### internal function: used in cross-validation
SE.egocentric.validation.auc <- function(A,train.index,test.index,r){
  Phat <- SE.egocentric(A,train.index,r=r)
  predictors <- as.numeric(Phat[test.index,-train.index])
  labels <- as.numeric(A[test.index,-train.index])
  
  #return(norm(labels-predictors,"F")^2)
  # test.index <- which(upper.tri(labels))
  # predictors <- predictors[test.index]
  # labels <- labels[test.index]
   if(length(unique(as.numeric(labels)))>1){
     auc.value <- as.numeric(auc(response=labels,predictor=predictors))
   }else{
     auc.value <- 0.5
   }
   if(auc.value<0.5) auc.value <- 0.5
   loss <- -auc.value
   
   return(loss)
  
}


###### cross-validation function to tune the rank
###### rmax: you search rank from r=1 to r=rmax
###### rho is the hold-out proportion for validation
###### B is the replication
SE.egocentric.CV <- function(A,train.index,rmax=5,B=3,rho=0.3){
  n <- nrow(A)
  n.train <- length(train.index)
  rmax = min(rmax,n.train-ceiling(rho*n.train)-1)
  cv.mat <- matrix(NA,B,rmax)
  for(bb in 1:B){
    validation.index <- sample(train.index,size=ceiling(rho*n.train))
    sub.train.index <- setdiff(train.index,validation.index)
    
    for(r in 1:rmax){
      cv.mat[bb,r] <- SE.egocentric.validation.auc(A,sub.train.index,validation.index,r)
    }
  }
  return(which.min(colMeans(cv.mat)))
  
}

################################################################################
# LE

#### general low-rank estimation
#### A is the adjacency matrix
#### train.index: the index of the observed nodes
#### r is the rank
LE.egocentric <- function(A,train.index,r=3){
  n <- nrow(A)
  n.train <- length(train.index)
  test.index <- setdiff(1:n,train.index)
  A11 <- A[train.index,train.index]
  A12 <- A[train.index,test.index]
  A21 <- A[test.index,train.index]
  A1 <- A[train.index,]
  SVD1 <- irlba(A11,nv=r,nu=r)
  if(r>1){
    Ahat11 <- SVD1$u[,1:r]%*%(t(SVD1$v[,1:r])*SVD1$d[1:r])
  }else{
    Ahat11 <- matrix(SVD1$u[,1:r],ncol=1)%*%(t(matrix(SVD1$v[,1:r],ncol=1))*SVD1$d[1:r])
  }
  # Ahat11 <- Ahat1[,train.index]
  # Ahat12 <- Ahat1[,test.index]
  # Ahat21 <- t(Ahat12)
  
  
  A11.est.inv <- ginv(Ahat11)
  A22.pred <- A21%*%A11.est.inv%*%A12
  A.new <- A
  A.new[train.index,train.index] <- A11
  A.new[train.index,test.index] <- A12
  A.new[test.index,train.index] <- A21
  A.new[test.index,test.index] <- A22.pred
  return(A.new)
  
}


##### internal function: used in cross-validation
LE.egocentric.validation.auc <- function(A,train.index,test.index,r){
  Phat <- LE.egocentric(A,train.index,r=r)
  predictors <- as.numeric(Phat[test.index,-train.index])
  labels <- as.numeric(A[test.index,-train.index])
  
  #return(norm(labels-predictors,"F")^2)
  # test.index <- which(upper.tri(labels))
  # predictors <- predictors[test.index]
  # labels <- labels[test.index]
  if(length(unique(as.numeric(labels)))>1){
    auc.value <- as.numeric(auc(response=labels,predictor=predictors))
  }else{
    auc.value <- 0.5
  }
  if(auc.value<0.5) auc.value <- 0.5
  loss <- -auc.value
  
  return(loss)
  
}


###### cross-validation function to tune the rank
###### rmax: you search rank from r=1 to r=rmax
###### rho is the hold-out proportion for validation
###### B is the replication
LE.egocentric.CV <- function(A,train.index,rmax=5,B=3,rho=0.3){
  n <- nrow(A)
  n.train <- length(train.index)
  rmax = min(rmax,n.train-ceiling(rho*n.train)-1)
  cv.mat <- matrix(NA,B,rmax)
  for(bb in 1:B){
    validation.index <- sample(train.index,size=ceiling(rho*n.train))
    sub.train.index <- setdiff(train.index,validation.index)
    
    for(r in 1:rmax){
      cv.mat[bb,r] <- LE.egocentric.validation.auc(A,sub.train.index,validation.index,r)
    }
  }
  return(which.min(colMeans(cv.mat)))
  
}
################################################################################
# LE+

#### general low-rank estimation
#### A is the adjacency matrix
#### train.index: the index of the observed nodes
#### r is the rank
LEplus.egocentric <- function(A,train.index,r=3){
  n <- nrow(A)
  n.train <- length(train.index)
  test.index <- setdiff(1:n,train.index)
  A11 <- A[train.index,train.index]
  A12 <- A[train.index,test.index]
  A21 <- A[test.index,train.index]
  A1 <- A[train.index,]
  # SE
  SVD_SE <- irlba(A1,nv=r,nu=r)
  if(r>1){
    Ahat1 <- SVD_SE$u[,1:r]%*%(t(SVD_SE$v[,1:r])*SVD_SE$d[1:r])
  }else{
    Ahat1 <- matrix(SVD_SE$u[,1:r],ncol=1)%*%(t(matrix(SVD_SE$v[,1:r],ncol=1))*SVD_SE$d[1:r])
  }
  Ahat11 <- Ahat1[,train.index]
  Ahat12 <- Ahat1[,test.index]
  Ahat21 <- t(Ahat12)
  
  # LE
  SVD_LE <- irlba(A11,nv=r,nu=r)
  if(r>1){
    Atilde11 <- SVD_LE$u[,1:r]%*%(t(SVD_LE$v[,1:r])*SVD_LE$d[1:r])
  }else{
    Atilde11 <- matrix(SVD_LE$u[,1:r],ncol=1)%*%(t(matrix(SVD_LE$v[,1:r],ncol=1))*SVD_LE$d[1:r])
  }
  
  A11.est.inv <- ginv(Ahat11)
  # A22.pred <- (Ahat21%*%ginv(Ahat11)%*%Ahat12+A21%*%ginv(Atilde11)%*%A12)/2
  A.new <- A
  A.new[train.index,train.index] <- (Ahat11+A11)/2
  A.new[train.index,test.index] <- (Ahat12+A12)/2
  A.new[test.index,train.index] <- (Ahat21+A21)/2
  A.new[test.index,test.index] <-
    (Ahat21%*%ginv(Ahat11)%*%Ahat12+A21%*%ginv(Atilde11)%*%A12)/2
  return(A.new)
  
}


##### internal function: used in cross-validation
LEplus.egocentric.validation.auc <- function(A,train.index,test.index,r){
  Phat <- LEplus.egocentric(A,train.index,r=r)
  predictors <- as.numeric(Phat[test.index,-train.index])
  labels <- as.numeric(A[test.index,-train.index])
  
  #return(norm(labels-predictors,"F")^2)
  # test.index <- which(upper.tri(labels))
  # predictors <- predictors[test.index]
  # labels <- labels[test.index]
  if(length(unique(as.numeric(labels)))>1){
    auc.value <- as.numeric(auc(response=labels,predictor=predictors))
  }else{
    auc.value <- 0.5
  }
  if(auc.value<0.5) auc.value <- 0.5
  loss <- -auc.value
  
  return(loss)
  
}


###### cross-validation function to tune the rank
###### rmax: you search rank from r=1 to r=rmax
###### rho is the hold-out proportion for validation
###### B is the replication
LEplus.egocentric.CV <- function(A,train.index,rmax=5,B=3,rho=0.3){
  n <- nrow(A)
  n.train <- length(train.index)
  rmax = min(rmax,n.train-ceiling(rho*n.train)-1)
  cv.mat <- matrix(NA,B,rmax)
  for(bb in 1:B){
    validation.index <- sample(train.index,size=ceiling(rho*n.train))
    sub.train.index <- setdiff(train.index,validation.index)
    
    for(r in 1:rmax){
      cv.mat[bb,r] <- LEplus.egocentric.validation.auc(A,sub.train.index,validation.index,r)
    }
  }
  return(which.min(colMeans(cv.mat)))
  
}