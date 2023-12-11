setwd("D:/experiment/Conference Paper/AAAI/AAAI2024/code")
rm(list = ls())

d_index <- 4

dpath           <- file.path("D:/experiment/AAAI/AAAI2024/dataset") 

Dataset         <- c("w8a", "ijcnn1_all","a9a_all","phishing","cod-rna","SUSY50000")

savepath1      <- paste0("D:/experiment/Conference Paper/AAAI/AAAI2024/Result/",
                         paste0("BOGD++-",Dataset[d_index],".txt")) 

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1]                                             

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

para1_setting <- list( 
  #  eta    = 2^{-2}, ## 2^0, 2^1, 2^2, 2^3,2^4
  eta    = 1000/sqrt(length_tr),
  lambda = 2^(3)/length_tr^2,
#  lambda = 0.001/sqrt(length_tr),
  gamma  = 2^{0},   ## 2^0, 2^1, 2^2, 2^3,2^4
  B      = 400
)
# -4 -3 -2 -1 0 1 2 3 4
sigma    <- 2^(1)

reptimes <- 2
runtime  <- c(rep(0, reptimes))
errorrate<- c(rep(0, reptimes))
N        <- c(rep(0, reptimes))

for(re in 1:reptimes)
{
  order    <- sample(1:length_tr,length_tr,replace = F)   #dis
  svpara   <- array(0,1)    # store the parameter for each support vector                                 # store the selected times of each kernel parameter
  svmat    <- matrix(0,nrow = feature_tr,ncol=1)
  t1       <- proc.time()  #proc.time()
  
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  svpara[1]  <- para1_setting$eta*ylabel[order[1]]
  k          <- 1
  
  for (i in 2:length_tr)
  {
    err <- 0
    diff_S_i  <- svmat - trdata[order[i], ]
    tem     <- colSums(diff_S_i*diff_S_i)
    fx        <- (svpara[1:k] %*% exp(tem/(-2*(sigma^2))))[1,1]
    hatyi     <- 1     
    if(fx < 0 )  
      hatyi <- -1
    if(hatyi != ylabel[order[i]])
      err  <- 1
    svpara <- (1-para1_setting$eta*para1_setting$lambda)*svpara
    if(fx*ylabel[order[i]] < 1)
    {
      N[re] <- N[re] +1
      if(k>= para1_setting$B)
      {
        abs_alpha   <- abs(svpara)
        s           <- (para1_setting$B-1)/sum(abs_alpha)
        samp_p      <- 1-s*abs_alpha
        dis         <- sample(1:para1_setting$B, 1, replace = T, prob = samp_p) #non-uniform sample
        
        svmat[,dis] <- trdata[order[i],]
        svpara      <- (1-para1_setting$eta*para1_setting$lambda)*svpara/(1-samp_p[dis])                 # svpara/(1/para1_setting$B)
        svpara[dis] <- para1_setting$eta*ylabel[order[i]]
        
#        projection
        threshold  <- para1_setting$eta*para1_setting$gamma
        abs_alpha  <- abs(svpara)
        vio_index  <- which(abs_alpha>threshold,arr.ind = TRUE)
        svpara[vio_index] <- threshold*sign(svpara[vio_index])
      }
      else
      {
        svmat       <- cbind(svmat,trdata[order[i],])
        svpara[k+1] <- para1_setting$eta*ylabel[order[i]]
        k  <- k+1
      }
    }
    error <- error + err   #record the err of selected superarm
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- sum(error)/length_tr
}

save_result <- list(
  note     = c(" the next term are:alg_name--dataname--sam_num--sigma--sv_num--gamma--eta--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_error"),
  alg_name = c("BOGD++"),
  dataname = paste0(Dataset[d_index], ".train"),
  Update   = mean(N),
  sv_num   = para1_setting$B,
  eta      = para1_setting$eta,
  lambda   = para1_setting$lambda,  
  gamma    = para1_setting$gamma,
  run_time = as.character(runtime),
  err_num  = as.character(errorrate), 
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter is %f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the mistake ratio is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of mistake ratio is %.5f in dataset", sd(errorrate))
sprintf("average updating time is %.5f in dataset", mean(N))