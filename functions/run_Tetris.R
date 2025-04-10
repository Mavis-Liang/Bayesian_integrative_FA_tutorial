<<<<<<< HEAD
source("Tetris.R")

# Define customized new_control() or getLambda() here if needed

run_Tetris <- function(sim_data, big_T=NULL, nrun=10000, burn=8000){
  S <- length(sim_data$N_s)
  set_alpha <- ceiling(1.25*S)

  if(is.null(big_T)){
    fit<- tetris(sim_data$Y_list, alpha=set_alpha, beta=1, nprint = 200, 
                 nrun=nrun, burn=burn)
    print("Start to choose big_T. It might take a long time.")
    big_T <- choose.A(fit, alpha_IBP=set_alpha, S=S)
  }
  
  big_T <- big_T
  run_fixed <- tetris(sim_data$Y_list, alpha=set_alpha, beta=1, 
                      fixed=TRUE, A_fixed=big_T, nprint = 200, nrun=nrun, burn=burn)
  
  return(run_fixed)
}
=======
source("Tetris.R")
run_Tetris <- function(sim_data, big_T=NULL){
  S <- length(sim_data$N_s)
  set_alpha <- ceiling(1.25*S)
  
  if(is.null(big_T)){
    fit<- tetris(sim_data$Y_list, alpha=set_alpha, beta=1)
    print("Start to choose big_T. It might take a long time.")
    big_T <- choose.A(fit, alpha_IBP=set_alpha, S=S)
  }
  
  big_T <- big_T
  run_fixed <- tetris(sim_data$Y_list, alpha=set_alpha, beta=1, 
                      fixed=TRUE, A_fixed=big_T)
  
  return(run_fixed)
}

>>>>>>> 5581b0f45e94b799518d5702b6c2dc3f4bbf6443
