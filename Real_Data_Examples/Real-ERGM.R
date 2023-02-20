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
P <- read.mat(paste0("data/",args[1],"-A.mat"))[[1]]

method = "ERGM"
for(j in 1:length(rho.seq)){
  rho <- rho.seq[j]
  rho.result <- list()
  rho.result <- foreach(i = 1:M)%dopar%{
    # for (i in 1:M){
    source("RModelsNoSS-Implementation.R")
    dt <- real.sample(P,rho=rho)

    A <- dt$A.obs
    start.time = Sys.time()
    bp <- ERGM.gwesp.pred(A)
    runtime = Sys.time()-start.time
    A.impute.out <- bp[-(1:dt$sample.n),-(1:dt$sample.n)]
    A.impute <- dt$A
    A.impute[-(1:dt$sample.n),-(1:dt$sample.n)] <- bp[-(1:dt$sample.n),-(1:dt$sample.n)]
    A.true <- dt$A
    A.true.out <- A.true[-(1:dt$sample.n),-(1:dt$sample.n)]
    print("Evaluating ERGM...")
    perf <- ego.perf(A.impute,A.true,A.impute.out,A.true.out)
    # tmp.result[[j]] <- perf
    perf[["runtime"]] <- runtime
    
    # tmp <- tmp.result
    tmp <- perf

  }
  # tmp<- tmp.result
  save(rho.result,file=paste0(method,"_Results/",args[1],"-",method,"-rho=",
                              rho,".Rda"))
}