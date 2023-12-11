setwd("D:/experiment/Conference Paper/AAAI/AAAI2024/code")
rm(list = ls())

d_index <- 3

dpath           <- file.path("D:/experiment/AAAI/AAAI2024/dataset")

Dataset         <- c("w8a", "ijcnn1_all","a9a_all","phishing","cod-rna","SUSY50000")

savepath1       <- paste0("D:/experiment/Conference Paper/AAAI/AAAI2024/Result/",paste0("NOGD-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1]                                             

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

# -4 -3 -2 -1 0 1 2 3 4
sigma    <- 2^(0)
B        <- 600
r        <- 0.2*B

reptimes <- 10
runtime  <- c(rep(0, reptimes))
errorrate<- c(rep(0, reptimes))
eta      <- 1000/sqrt(length_tr)

for(re in 1:reptimes)
{
  order    <- sample(1:length_tr,length_tr,replace = F)   #dis
  svpara   <- array(0,1)    # store the parameter for each support vector
  svmat    <- matrix(0,nrow = feature_tr,ncol=1)
  Gram     <- matrix(1,nrow = B,ncol=B)  
  
  t1       <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  svpara[1]  <- eta*ylabel[order[1]]
  k          <- 1
  
  for (t in 2:length_tr)
  {
    err <- 0
    if(k<B)
    {
      diff_S_i <- svmat - trdata[order[t], ]
      tem     <- colSums(diff_S_i*diff_S_i)
      kt      <- exp(tem/(-2*(sigma^2)))
      fx      <- crossprod(svpara[1:k],kt)[1,1]
      hatyi   <- 1     
      if(fx < 0 )  
        hatyi <- -1
      if(hatyi != ylabel[order[t]])
      {
        err   <- 1
      }
      if(fx*ylabel[order[t]]<1)
      {
        svmat       <- cbind(svmat,trdata[order[t],])
        svpara[k+1] <- eta*ylabel[order[t]]
        k           <- k+1
        Gram[k,1:(k-1)] <- kt
        Gram[1:(k-1),k] <- kt
      }
    }
    if(k==B)
    {
      T <- svd(Gram)
#      U <- T$u
      D <- T$d  ##3 vector
      V <- T$v

      D <- D^{-0.5}
      
      tem <- D*V
      tem <- solve(tem)
      w_0 <- svpara%*%tem
      w_0 <- t(w_0)
      w_t <- w_0[1:r]
      
#      Diag <- diag(D[1:r])
      V    <- t(V[,1:r])
      D    <- D[1:r]
    }
    if(k>=B)
    {
      k   <- k+1
      
      diff_St <- svmat - trdata[order[t], ]
      tem     <- colSums(diff_St*diff_St)
      kt      <- exp(tem/(-2*(sigma^2)))
      tem     <- V%*%kt
      ph_t    <- D*tem
      
      fx      <- crossprod(w_t,ph_t)[1,1]
      hatyt   <- 1     
      if(fx < 0 )  
        hatyt <- -1
      if(hatyt != ylabel[order[t]])
      {
        err   <- 1
      }
      if(ylabel[order[t]]*fx<1)
        w_t   <- w_t+eta*ylabel[order[t]]*ph_t
    }
    error <- error + err   #record the err of selected superarm
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
}

save_result <- list(
  note     = c(" the next term are:alg_name--dataname--sam_num--sigma--sv_num--gamma--eta--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_error"),
  alg_name = c("NOGD"),
  dataname = paste0(Dataset[d_index], ".train"),
  sam_num  = length_tr,
  sv_num   = B,
  run_time = as.character(runtime),
  err_num  = as.character(errorrate), 
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = F) 

sprintf("the candidate kernel parameter is %f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("the average running time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of average error rate is %.5f in dataset", sd(errorrate))
