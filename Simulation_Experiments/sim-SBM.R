# Please make sure to set your directory, for example,
# setwd("/sample/Simulation_Experiments")
rho.seq <- c(0.7,0.9,0.5,0.2,0.1,0.3) # Sequence of rhos
deg.seq <- c(20,50) # Sequence of degrees
M <- 100 # Number of iterations
source("RModelsNoSS-Implementation.R")

library(doParallel)
registerDoParallel(cores=20)
# registerDoParallel(cores=2)

method = "SBM"

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least one argument must be supplied.n", call.=FALSE)
}
if (length(args)>=2){
  M <- args[2]
}

P <- read.mat(paste0("data/",args[1],"-P.mat"))[[1]]


for(j in 1:length(rho.seq)){
  rho <- rho.seq[j]
  for (deg in deg.seq){
    rho.result <- list()
    rho.result <- foreach(i = 1:M)%dopar%{
    # for (i in 1:M){
      source("RModelsNoSS-Implementation.R")
      # tmp.result <- list()
  
      dt <- ego.sample(P,rho=rho,deg=deg)
  
      A <- dt$A.obs
      A.impute <- dt$A
      A.true <- dt$A
      A.true.out <- A.true[-(1:dt$sample.n),-(1:dt$sample.n)]
  
      n <- nrow(A)
      holdout.nodes <- (1:n)[-(1:dt$sample.n)]
  
  
      train.index <- 1:dt$sample.n
      Rmat <- A[1:dt$sample.n,]
      SR <- norm(Rmat,"F")^2/norm(Rmat,"2")^2
      k.seq <- rep(NA,20)
      for(k in 1:20){
        k.seq[k] <- suppressMessages(
          Block.egocentric.CV(A,train.index,Kmax=SR+10)
        )
      }
      kblock <- floor(mean(k.seq))
      if(kblock<2){
        kblock <- 2
      }
      start.time = Sys.time()
      block.fit <- Block.egocentric.fit(A,train.index,K=kblock)
      runtime = Sys.time()-start.time
      A.impute.out <- block.fit[-(1:dt$sample.n),-(1:dt$sample.n)]
      A.impute <- dt$A
      A.impute[-(1:dt$sample.n),-(1:dt$sample.n)] <- block.fit[-(1:dt$sample.n),-(1:dt$sample.n)]
      print("Evaluating cross-validated SBM ...")
      SBM.perf <- ego.perf(A.impute,A.true,dt$P,A.impute.out,A.true.out,dt$P.out)
      SBM.perf[["runtime"]] <- runtime
      # tmp.result[[j]] <- SBM.perf
    
      # tmp<- tmp.result
      tmp <- SBM.perf
      # rho.result[[i]] <- tmp.result
    }
    
    save(rho.result,file=paste0(method,"_Results/",args[1],"-",method,"-rho=",
                                rho,"-deg=",deg,"-time.Rda"))
  }
}