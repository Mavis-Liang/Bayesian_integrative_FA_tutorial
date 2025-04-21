if (!require("caret", quietly = TRUE)) {
  install.packages("caret")
}
if (!require("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
if (!require("matlab", quietly = TRUE)) {
  install.packages("matlab")
}
if (!require("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
# simulate data for scenario 4 in the manuscript: mimicing nutritional data, based on Tetris
gen_scenario4 <- function(S = 12, N_s = c(1362, 217, 417, 1012, 2241, 205, 2403, 3775, 1790, 761, 373, 465), P = 42) {
  N <- sum(N_s)
  # Define probabilities for categorical variables based on the table
  age_probs <- c("18-44" = 0.593, "45-55" = 0.190, "55-74" = 0.216)
  sex_probs <- c("Female" = 0.525, "Male" = 0.475)
  immigrant_gen_probs <- c("First" = 0.772, "Second or Later" = 0.228)
  years_lived_probs <- c(
    "Born in mainland United States" = 0.283, 
    "< 10 y" = 0.502, 
    "≥10 y" = 0.215
  )
    employment_probs <- c(
      "Retired and not currently employed" = 0.082,
      "Not retired and not currently employed" = 0.411,
      "Part-time (≤35 hrs)" = 0.168,
      "Full-time (35+ hrs)" = 0.339
    )
    marital_status_probs <- c(
      "Single" = 0.338,
      "Married/living with partner" = 0.499,
      "Separated/divorced/widowed" = 0.163
    )
    income_probs <- c("< $30k" = 0.617, "≥ $30k" = 0.318, "Not reported" = 0.064)
    education_probs <- c(
      "Less than high school" = 0.33,
      "High school or equivalent" = 0.284,
      "Greater than high school" = 0.386
    )
    physical_activity_probs <- c("No" = 0.576, "Yes" = 0.424)
    dietary_supplement_probs <- c("No" = 0.593, "Yes" = 0.407)
    bmi_probs <- c(
      "Underweight" = 0.012,
      "Normal" = 0.216,
      "Overweight" = 0.374,
      "Obesity" = 0.398
    )

    # Simulate energy intake as a continuous variable
    # Mean = 1910.1 kcal/day, SD = 10.8 kcal/day (based on the table)
    energy_intake <- rnorm(N, mean = 1910.1, sd = 10.8)

    # Generate data for each variable
    X <- data.frame(
      Age_Group = sample(names(age_probs), size = N, replace = TRUE, prob = age_probs),
      Sex = sample(names(sex_probs), size = N, replace = TRUE, prob = sex_probs),
      Immigrant_Generation = sample(names(immigrant_gen_probs), size = N, replace = TRUE, prob = immigrant_gen_probs),
      Years_Lived = sample(names(years_lived_probs), size = N, replace = TRUE, prob = years_lived_probs),
      Employment = sample(names(employment_probs), size = N, replace = TRUE, prob = employment_probs),
      Marital_Status = sample(names(marital_status_probs), size = N, replace = TRUE, prob = marital_status_probs),
      Yearly_Income = sample(names(income_probs), size = N, replace = TRUE, prob = income_probs),
      Education_Status = sample(names(education_probs), size = N, replace = TRUE, prob = education_probs),
      Physical_Activity = sample(names(physical_activity_probs), size = N, replace = TRUE, prob = physical_activity_probs),
      Dietary_Supplement = sample(names(dietary_supplement_probs), size = N, replace = TRUE, prob = dietary_supplement_probs),
      BMI = sample(names(bmi_probs), size = N, replace = TRUE, prob = bmi_probs),
      Energy_Intake = energy_intake
    )

    # One-hot encode categorical variables
    one_hot_X <- caret::dummyVars(" ~ .", data = X)
    X <- predict(one_hot_X, newdata = X) %>% as.data.frame()
    # Normalize continuous variables
    X$Energy_Intake <- scale(X$Energy_Intake) 
    cumulative_sizes <- c(0, cumsum(N_s))
    X_list <- lapply(1:12, function(s) {
      X_s <- X[(cumulative_sizes[s] + 1):cumulative_sizes[s + 1], ]
      X_s
    })

    # Define beta coefficients
    Q <- ncol(X)
    beta <- matrix(runif(P * Q, -1, 1), nrow = P, ncol = Q)

    # Define "big-T" matrix
    K <- 4 + 7 + 12
    big_T <- matrix(0, nrow = S, ncol = K)
    # select common factors
    big_T[, 1:4] <- 1
    # select partial factors
    BX <- c(1,2,6)
    CA <- c(2, 3, 4)
    CHI <- c(3, 7, 10, 11)
    MIA <- c(4, 5, 12)
    M <- c(6, 7, 8)
    PR <- c(9, 10)
    SA <- c(11, 12)
    big_T[BX,5] <- big_T[CA, 6] <- big_T[CHI, 7] <- big_T[MIA, 8] <- big_T[M, 9] <- big_T[PR, 10] <- big_T[SA, 11] <- 1

    # Select study-specific factors
    for (s in 1:S) {
      big_T[s, 11 + s] <- 1
    }

    # Define T_s
    T_s <- lapply(1:S, function(s) {
      T_s <- diag(big_T[s, ])
      T_s
    })


    # Loading matrix
    sparsity <- 0.4
    Phi_long <- as.vector(matlab::zeros(P, K))
    noZERO_count <- P * K * (1 - sparsity)
    noZEROs <- runif(noZERO_count, 0.6, 1)
    sign <- sample(x = length(noZEROs), 
                    size = (length(noZEROs) / 2))# Randomly assign negative sign
    noZEROs[sign] <- noZEROs[sign] * (-1)
    position_noZERO <- sample(x = K * P, size = length(noZEROs))
    Phi_long[position_noZERO] <- noZEROs
    Phi_star <- matrix(Phi_long, P, K)

    # Generate Sigma_s
    Psi_list <- list()
    Psi_list <- lapply(1:S, function(s) {
      Psi_s <- diag(runif(P, 0, 1), P)
      Psi_s
    })
    Sigma_list <- list()
    Sigma_list <- lapply(1:S, function(s) {
      Sigma_s <- Phi_star %*% T_s[[s]] %*% t(Phi_star) + Psi_list[[s]]
      Sigma_s
    })

    # Generate Y_s
    Y_list <- lapply(1:S, function(s) {
      mu_s <- as.matrix(X_list[[s]]) %*% t(beta)
      Sigma_s <- Sigma_list[[s]]
      E_s <- MASS::mvrnorm(n = N_s[[s]], mu = rep(0, P), Sigma = Sigma_s)
      Y_s <- mu_s + E_s
    })

    Y_mat <- do.call(rbind, Y_list)
    M_list <- lapply(1:S, function(s) {
      matrix(1, nrow = N_s[s], ncol = 1)
    })
    M <- Matrix::bdiag(M_list) %>% as.matrix()
    Phi <- as.matrix(Phi_star[,colSums(big_T)==S])
    SigmaPhi = tcrossprod(Phi)
    P = diag((colSums(big_T) == S)*1)
    LambdaList <- lapply(1:S, function(s){
        T_s <- diag(big_T[s,])
        Lambda_s <- Phi_star %*% (T_s - P)
        Lambda_s <- Lambda_s[,-which(colSums(Lambda_s == 0) == nrow(Lambda_s))]
        Lambda_s <- matrix(Lambda_s, nrow=nrow(Phi_star))}
        )
    SigmaLambdaList <- lapply(1:S, function(s){
        tcrossprod(LambdaList[[s]])})
    Psi_mat <- lapply(Psi_list,diag) %>% do.call(cbind, .)

    return(list(Y_mat=Y_mat, Y_list=Y_list, N_s=N_s, M=M, 
                  X=as.matrix(X), Beta=beta,
                  Phi=Phi, SigmaPhi = tcrossprod(Phi), 
                  LambdaList=LambdaList,
                  SigmaLambdaList = SigmaLambdaList,
                  SigmaMarginal = Sigma_list,
                  Psi_list=Psi_list, 
                  Psi_mat=Psi_mat, big_T=big_T))
}
