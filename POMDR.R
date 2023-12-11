setwd("D:/experiment/Conference Paper/AAAI/AAAI2024/code")
rm(list = ls())

d_index <- 3

dpath           <- file.path("D:/experiment/AAAI/AAAI2024/dataset")  

Dataset         <- c("w8a", "ijcnn1_all","a9a_all","phishing","cod-rna","SUSY50000")

savepath1      <- paste0("D:/experiment/Conference Paper/AAAI/AAAI2024/Result/",
                         paste0("POMDR-M-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

##############################################################################

sigma     <- 1
B         <- 600
B0        <- round(10*log(length_tr))
U         <- 25
coe       <- 0.1
M         <- 15

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))
sumtdelta <- c(rep(0, reptimes))
bar_t     <- c(rep(length_tr, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  #  order      <- c(1:length_tr)
  k          <- 0
  error      <- 0
  alpha      <- 1/length_tr^(1/2)
  alpha_t    <- 0
  sum_delta  <- 0
  sum_t_delta<- 0
  lambda_t   <- U/sqrt(3)
  Norm       <- 0
  Norm2      <- 0
  f_t_1      <- 0

  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  
  Inver_K    <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram       <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram_B     <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  beta_ast   <- array(0,1)                          # The optimal parameter d
  delta      <- 0                                   # The difference of f''-f'
  kt         <- array(0,1)
  indx       <- array(0,1)
  
  #sv_index   <- array(0,1)
  svpara     <- array(0,1)
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  #sv_index[1]<- order[1]
  svpara[1]  <- lambda_t*ylabel[order[1]]
  k          <- 1
  Inver_K[1,1] <- 1
  Gram[1,1]    <- 1 
  Gram_B[1,1]  <- 1
  sum_t_delta  <- sum_t_delta+1
  Norm  <- lambda_t
  i     <- 1
  
  ### from the second instance
  while(bar_t[re]>=length_tr && i< length_tr)
  {
    i=i+1
    if(i <= M)
    {
      diff  <- t(trdata[order[1:(i-1)], ])- trdata[order[i], ]
      if(i==2)
      {
        tem <- crossprod(diff[1,],diff[1,])[1,1]
      }
      else
      {
        tem <- colSums(diff*diff)
      }
      kt_1  <- exp(tem/(-2*(sigma)^2))
      k_t_M <- 1/(i-1)*crossprod(ylabel[order[1:(i-1)]],kt_1)[1,1]
    }else{
      diff  <- t(trdata[order[(i-M):(i-1)], ])- trdata[order[i], ]
      tem   <- colSums(diff*diff)
      kt_1  <- exp(tem/(-2*(sigma)^2))
      k_t_M <- 1/M*crossprod(ylabel[order[(i-M):(i-1)]],kt_1)[1,1]
    }

    lambda_t   <- U/sqrt(3+sum_t_delta)

    diff  <- svmat- trdata[order[i], ]
    tem   <- colSums(diff*diff)
    kt    <- exp(tem/(-2*(sigma)^2))
    f_t_1 <- crossprod(svpara[1:k],kt)[1,1]
    
    fx    <- f_t_1 + lambda_t*k_t_M
    hatyi <- 1
    if(fx < 0)
      hatyi  <- -1
    if(hatyi != ylabel[order[i]])
    {
      error <- error + 1
    }
    if(ylabel[order[i]]*fx<1)
    {
      #### compute alpha_t

      delta_t    <- max(1-2*ylabel[order[i]]*k_t_M,0)
      sumtdelta[re] <- sumtdelta[re] + delta_t
      
      beta_ast   <- Inver_K%*%kt
      alpha_t    <- 1-crossprod(beta_ast,kt)[1,1]
      if(alpha_t<0)
        alpha_t  <- 0
      sq_alpha_t <- sqrt(alpha_t)
      if(sq_alpha_t <= alpha)
#      if(sq_alpha_t <= 0.3)
      {
        tem   <- Gram%*%beta_ast
        tem1  <- crossprod(tem,beta_ast)[1,1]
        tem2 <- 0
        if(i<=M)
        {
          if(i==2)
          {
            tem2 <- tem2 + ylabel[order[1]]*beta_ast
          }else{
            for(r in 1:(i-1))
            {
              diff  <- svmat- trdata[order[r], ]
              kt_    <- exp(colSums(diff*diff)/(-2*(sigma)^2))
              tem2  <- tem2+ylabel[order[r]]*crossprod(kt_,beta_ast)[1,1]
            }
          }
          tem2 <- tem2/(i-1)
        }else{
          for(r in 1:M)
          {
            diff  <- svmat- trdata[order[i-r], ]
            kt_   <- exp(colSums(diff*diff)/(-2*(sigma)^2))
            tem2  <- tem2+ylabel[order[i-r]]*crossprod(kt_,beta_ast)[1,1]
          }
          tem2 <- tem2/M
        }
        sum_t_delta <- sum_t_delta+max(tem1-2*ylabel[order[i]]*tem2,0)
        tem3        <- crossprod(tem,svpara)[1,1]
        svpara      <- svpara + lambda_t*ylabel[order[i]]*as.vector(beta_ast)
        Norm        <- sqrt(Norm^2+2*ylabel[order[i]]*lambda_t*tem3+lambda_t^2*tem1)
        if(Norm >U)
        {
          svpara <- svpara*U/Norm
          Norm <- U
        }
      }else{
        k           <- k+1
        svmat       <- cbind(svmat,trdata[order[i],])
#        sv_index[k] <- order[i]
        svpara[k]   <- lambda_t*ylabel[order[i]]
        
#        update the inverse kernel matrix 
        tem_d       <- beta_ast
        tem_d[k]    <- -1
        incre       <- tem_d %*% t(tem_d)/alpha_t
        incre[1:(k-1),1:(k-1)] <- incre[1:(k-1),1:(k-1)]+Inver_K
        
        Inver_K     <- incre
        Gram        <- cbind(Gram,as.vector(kt))
        Gram        <- rbind(Gram,c(kt,1))
        if(k<=(B/2))
        {
          Gram_B <- cbind(Gram_B,as.vector(kt))
          Gram_B <- rbind(Gram_B,c(kt,1))
        }
        Norm     <- sqrt(Norm^2+2*ylabel[order[i]]*lambda_t*f_t_1+lambda_t^2)
        if(Norm >U)
        {
          svpara <- svpara*U/Norm
          Norm   <- U
        }
        sum_t_delta <- sum_t_delta+max(1-2*ylabel[order[i]]*k_t_M,0)
      }
    }
    if(k==B0)
      bar_t[re] <- i
  }
  flag <- 0
  if(bar_t[re]<length_tr)
  {
  #  sum_t_delta<- 0
    for(t in (i+1):length_tr)
    {
      diff  <- t(trdata[order[(t-M):(t-1)], ]) - trdata[order[t], ]
      kt_1  <- exp(colSums(diff*diff)/(-2*(sigma)^2))
      k_t_M <- 1/M*crossprod(ylabel[order[(t-M):(t-1)]],kt_1)[1,1]
      if(flag == 0)
        lambda_t <- U/sqrt(3+sum_t_delta)
      else{
        lambda_t <- coe*U/sqrt(3+sum_t_delta)
      }
#      lambda_t   <- 0.1
      
      diff  <- svmat- trdata[order[t], ]
      kt    <- exp(colSums(diff*diff)/(-2*(sigma)^2))
      f_t_1 <- crossprod(svpara[1:k],kt)[1,1]
      fx    <- f_t_1 + lambda_t*k_t_M
      hatyi <- 1
      if(fx < 0)
        hatyi  <- -1
      if(hatyi != ylabel[order[t]])
      {
        error <- error + 1
      }
      if(ylabel[order[t]]*fx<1)
      {
        #### compute alpha_t
        k<- k+1
        svmat     <- cbind(svmat,trdata[order[t],])
#        sv_index[k] <- order[t]
        svpara[k] <- lambda_t*ylabel[order[t]]
        Norm      <- sqrt(Norm^2+2*ylabel[order[t]]*lambda_t*f_t_1+lambda_t^2)
        if(Norm > U)
        {
          svpara  <- svpara*U/Norm
          Norm    <- U
        }
        sum_t_delta  <- sum_t_delta+max(1-2*ylabel[order[t]]*k_t_M,0)
        sumtdelta[re] <- sumtdelta[re] + max(1-2*ylabel[order[t]]*k_t_M,0)
        if(k<=(B/2))
        {
          Gram_B <- cbind(Gram_B,as.vector(kt))
          Gram_B <- rbind(Gram_B,c(kt,1))
        }else{
          kt_    <- kt[1:(B/2)]
          max_kt   <- max(kt_)
          tem_     <- which(kt_==max_kt,arr.ind = TRUE)
          indx[k-(B/2)]<- tem_[1]
        }
        if(k==B)
        {
          flag <- 1
          tem_svmat  <- svmat
          svmat      <- matrix(0,nrow = feature_tr,ncol=B/2)
          svmat      <- tem_svmat[,1:(B/2)]
          
          tem_svpara <- svpara
          svpara     <- array(0,B/2)
          svpara     <- tem_svpara[1:(B/2)]
          
          sum_t_delta<- 0
          k          <- B/2
          Norm2      <- 0 
          for(r in (k+1):B)
          {
            ind <- indx[r-k]
            svpara[ind]   <- svpara[ind] + tem_svpara[r]
          }
          tem    <- Gram_B%*%svpara
          Norm2  <- crossprod(tem,svpara)[1,1]
          
          svpara <- svpara*U/sqrt(Norm2)
          Norm   <- U
        }
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
  alg_name = c("POMDR-M-"),
  dataname = paste0(Dataset[d_index], ".train"),
  ker_para = sigma,
  sv_num   = sum(All_bud)/re,
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time      <- sd(runtime),
  sd_err       <-sd(errorrate)
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
sprintf("average A_T is %.5f in dataset", mean(sumtdelta))
sprintf("average bar_t is %.5f in dataset", mean(bar_t))
