# Please make sure to set your directory, for example,
# setwd("/sample/Simulation_Experiments")
rho.seq <- c(0.7,0.9,0.5,0.2,0.1,0.3) # Sequence of rhos
deg.seq <- c(20,50) # Sequence of degrees
M <- 100 # Number of iterations
source("RmethodsHeader-sim.R")

library(doParallel)
registerDoParallel(cores=20)
# registerDoParallel(cores=2)

method = "ERGM"

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
      start.time = Sys.time()
      bp <- ERGM.gwesp.pred(A)
      runtime = Sys.time()-start.time
      A.impute.out <- bp[-(1:dt$sample.n),-(1:dt$sample.n)]
      A.impute <- dt$A
      A.impute[-(1:dt$sample.n),-(1:dt$sample.n)] <- bp[-(1:dt$sample.n),-(1:dt$sample.n)]
      A.true <- dt$A
      A.true.out <- A.true[-(1:dt$sample.n),-(1:dt$sample.n)]
      print("Evaluating ERGM...")
      perf <- ego.perf(A.impute,A.true,dt$P,A.impute.out,A.true.out,dt$P.out)
      # tmp.result[[j]] <- perf
      perf[["runtime"]] <- runtime
      
      # tmp <- tmp.result
      tmp <- perf
      # rho.result[[i]]<- tmp.result
    }
    if (!dir.exists(paste0(method,"_Results/"))){
      dir.create(paste0(method,"_Results/"))
    }
    save(rho.result,file=paste0(method,"_Results/",args[1],"-",method,"-rho=",
                                rho,"-deg=",deg,"-time.Rda"))
  }
}