# Please make sure to set your directory, for example,
# setwd("/sample/Real_Data_Examples")
rho.seq <- c(0.5,0.3,0.2,0.1) # Sequence of rhos
M <- 100 # Number of iterations
source("RmethodsHeader-Real.R")

library(doParallel)
registerDoParallel(cores=20)
# registerDoParallel(cores=2)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least one argument must be supplied.n", call.=FALSE)
}
if (length(args)>=2){
  M <- args[2]
}
A <- read.mat(paste0("data/",args[1],"-A.mat"))[[1]]

rho.result <- list()
method = "LE"
for(j in 1:length(rho.seq)){
  rho <- rho.seq[j]
  # tmp.result <- list()
  rho.result <- list()
  rho.result <- foreach(i = 1:M)%dopar%{
    source("RModelsNoSS-Implementation.R")
    dt <- real.sample(A,rho=rho)
    
    A.obs <- dt$A.obs
    A.impute <- dt$A
    A.true <- dt$A
    A.true.out <- A.true[-(1:dt$sample.n),-(1:dt$sample.n)]
    
    n <- nrow(A)
    holdout.nodes <- (1:n)[-(1:dt$sample.n)]
    
    
    train.index <- 1:dt$sample.n
    Rmat <- A.obs[1:dt$sample.n,]
    SR <- norm(Rmat,"F")^2/norm(Rmat,"2")^2
    k.seq <- rep(NA,20)
    for(k in 1:20){
      k.seq[k] <- suppressMessages(
        LE.egocentric.CV(A.obs,train.index,rho = 0.1,rmax=min(ceiling(SR+10),ceiling(n^(1/2))))
      )
    }
    kblock <- round(mean(k.seq))
    if(kblock<2){
      kblock <- 2
    }
    start.time = Sys.time()
    block.fit <- suppressMessages(LE.egocentric(A.obs,train.index,kblock))
    runtime = Sys.time()-start.time
    A.impute.out <- block.fit[-(1:dt$sample.n),-(1:dt$sample.n)]
    A.impute.out <- block.fit[-(1:dt$sample.n),-(1:dt$sample.n)]
    A.impute <- dt$A
    A.impute[-(1:dt$sample.n),-(1:dt$sample.n)] <- block.fit[-(1:dt$sample.n),-(1:dt$sample.n)]
    print("Evaluating cross-validated LE ...")
    SBM.perf <- ego.perf(A.impute,A.true,A.impute.out,A.true.out)
    SBM.perf[["runtime"]] <- runtime
    # tmp.result[[j]] <- SBM.perf
    tmp <- SBM.perf
  }
  save(rho.result,file=paste0(method,"_Results/",args[1],"-",method,"-rho=",
                              rho,".Rda"))
}