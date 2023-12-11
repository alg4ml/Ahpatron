setwd("D:/experiment/Conference Paper/AAAI/AAAI2024/code")
rm(list = ls())


d_index <- 6

dpath           <- file.path("D:/experiment/AAAI/AAAI2024/dataset") 

Dataset         <- c("w8a", "ijcnn1_all","a9a_all","phishing","cod-rna","SUSY50000")

savepath1       <- paste0("D:/experiment/Conference Paper/AAAI/AAAI2024/Result/",
                          paste0("Ahpatron-",Dataset[d_index],".txt"))

traindatapath   <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata          <- traindatamatrix[ ,-1]
ylabel          <- traindatamatrix[ ,1]                                        

length_tr       <- nrow(trdata)                                               
feature_tr      <- ncol(trdata)              

##############################################################################

sigma     <- 2^(0)
B         <- 600
U         <- sqrt(B)/2
eta       <- 0.0005
varepsilon <- 0.5

# eta = 0.01 w8a

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
N         <- c(rep(0, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  k          <- 0
  
  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  svpara     <- array(0,1)
  inver_Gram_1   <- matrix(0,nrow = B/2,ncol=B/2)
  Gram       <- matrix(1,nrow = B,ncol=B) 
  Gram_1     <- matrix(1,nrow = B/2,ncol=B/2)          
  Gram_12    <- matrix(0,nrow = B/2,ncol=B/2)  
  Diag       <- diag(eta,nrow=B/2,ncol=B/2)
  
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  lambda     <- 1*U/sqrt(4*B)
#  varepsilon <- lambda/0.5
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  svpara[1]  <- lambda*ylabel[order[1]]
  k          <- 1
  norm_f     <- lambda
  bar_t      <- 0
  
  
  ### from the second instance
  for(t in 2:length_tr)
  {
    diff  <- svmat- trdata[order[t], ]
    kt    <- exp(colSums(diff*diff)/(-2*(sigma)^2))
    fx    <- crossprod(svpara[1:k],kt)[1,1]
    hatyi <- 1
    if(fx < 0)
      hatyi <- -1
    if(hatyi != ylabel[order[t]])
    {
      error <- error + 1
    }
    if(fx*ylabel[order[t]]<1-varepsilon)
    {
      if(fx*ylabel[order[t]]>0)
        N[re] <- N[re] +1
      if(k<B)
      {
        k           <- k+1
        svmat       <- cbind(svmat,trdata[order[t],])
        svpara[k]   <- lambda*ylabel[order[t]]
        norm_f      <- sqrt(norm_f^2+lambda^2+2*lambda*ylabel[order[t]]*fx)
        if(norm_f>U)
        {
          svpara    <- svpara*U/norm_f
          norm_f    <- U
        }
        Gram[k,1:(k-1)] <- kt
        Gram[1:(k-1),k] <- kt
      }else{
        k  <- B/2
        
        abs_svpara <- abs(svpara)
        med        <- median(abs_svpara)
        
        index        <- which(abs_svpara>=med,arr.ind = TRUE)
        large_index  <- index[1:(B/2)]
        small_index  <- setdiff(c(1:B),large_index)
        
        large_svpara <- svpara[large_index]
        small_svpara <- svpara[small_index]
        
        svmat        <- svmat[,large_index]
        
        Gram_1       <- Gram[large_index,large_index]

        inver_Gram_1 <- solve(Gram_1+Diag)
        Gram_12      <- Gram[large_index,small_index]
        
        tem     <- Gram_12%*%small_svpara
        theta   <- inver_Gram_1%*%tem
        
        svpara  <- large_svpara+theta   
        
        ############# projection f_1+f_2 ##########
        tem     <- Gram_1%*%svpara
        
        Norm    <- sqrt(crossprod(svpara,tem)[1,1])

        svpara  <- svpara*norm_f/Norm
        Norm    <- norm_f
                  
        kt_large       <- kt[large_index] 
        tem_fx         <- crossprod(svpara,kt_large)[1,1]
        
        svmat          <- cbind(svmat,trdata[order[t],])
        svpara[k+1]    <- lambda*ylabel[order[t]]
        k              <- k+1
        
        norm_f     <- sqrt(Norm^2+lambda^2+2*lambda*ylabel[order[t]]*tem_fx)
        if(norm_f>U)
        {
          svpara    <- svpara*U/norm_f
          norm_f    <- U
        }
        Gram[1:(B/2),1:(B/2)]  <-Gram_1
        Gram[k,1:(k-1)] <- kt_large
        Gram[1:(k-1),k] <- kt_large
      }
    }
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("Ahpatron"),
  dataname = paste0(Dataset[d_index], ".train"),
  varepsilon   = varepsilon,
  ker_para = sigma,
  sv_num   = B,
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
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of error is %.5f in dataset", sd(errorrate))
sprintf("average updating time is %.5f in dataset", mean(N))
