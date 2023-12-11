setwd("D:/experiment/Conference Paper/AAAI/AAAI2024/code")
rm(list = ls())

d_index <- 2

dpath           <- file.path("D:/experiment/AAAI/AAAI2024/dataset")  

Dataset         <- c("w8a", "ijcnn1_all","a9a_all","phishing","cod-rna","SUSY50000")

savepath1      <- paste0("D:/experiment/Conference Paper/AAAI/AAAI2024/Result/",
                         paste0("Projectron-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

##############################################################################

# -4 -3 -2 -1 0 1 2 3 4
sigma     <- 2^(4)
#B         <- 600
#U         <- sqrt(B)/3

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  k          <- 0
  norm_ft    <- 0
  error      <- 0
  
  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  
  Inver_K    <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram       <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  d_ast      <- array(0,1)                          # The optimal parameter d
  delta      <- 0                                   # The difference of f''-f'
  kt         <- array(0,1)
  
  sv_index   <- array(0,1)
  svpara     <- array(0,1)
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  sv_index[1]<- order[1]
  svpara[1]  <- ylabel[order[1]]
  k          <- 1
  Inver_K[1,1] <- 1
  Gram[1,1]    <- 1 
  
  ### from the second instance
  for (i in 2:length_tr)
  {
    err     <- 0
    diff    <- svmat[,1:k]- trdata[order[i], ]
    if(k>1)
    {
      tem   <- colSums(diff*diff)
    }else{
      tem   <- sum((svmat[,1:k]- trdata[order[i], ])*(svmat[,1:k]- trdata[order[i], ]))
    }
    kt    <- exp(tem/(-2*(sigma)^2))
    sum   <- crossprod(svpara[1:k],kt)
    fx <- sum[1,1]
    hatyi <- 1
    if(fx < 0)
      hatyi  <- -1
    if(hatyi != ylabel[order[i]])
    {
      error <- error + 1
      #### compute delta
      d_ast <-Inver_K%*%kt
      inter <- 1-crossprod(d_ast,kt)[1,1]
      if(abs(1-crossprod(d_ast,kt)[1,1])<1e-10)
        inter <- 0
      delta <- sqrt(inter)
      eta <- 0.9
      if(delta <= eta)
      {
        svpara <- svpara + ylabel[order[i]]*as.vector(d_ast)
      }else{
        k           <- k+1
        svmat       <- cbind(svmat,trdata[order[i],])
        sv_index[k] <- order[i]
        svpara[k]   <- ylabel[order[i]]

        #update the inverse kernel matrix 
        tem_d       <- d_ast
        tem_d[k]    <- -1
        incre       <- tem_d %*% t(tem_d)/(1-crossprod(d_ast,kt)[1,1])
        incre[1:(k-1),1:(k-1)] <- incre[1:(k-1),1:(k-1)]+Inver_K
        Inver_K     <- incre
        Gram        <- cbind(Gram,kt)
        Gram        <- rbind(Gram,c(kt,1))      
      }
    }
  }
  t2 <- proc.time()
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  All_bud[re]   <- k
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("Projectron-"),
  dataname = paste0(Dataset[d_index], ".train"),
  eta  = eta,
  ker_para = sigma,
  sv_num   = sum(All_bud)/re,
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of support vectors is %d", round(sum(All_bud)/re))
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of error is %.5f in dataset", sd(errorrate))

