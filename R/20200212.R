#------------------------------------------------------------#
# Code for "Optimal prediction in an additive 
#          functional model"
#------------------------------------------------------------#
library(ggplot2)

#----------------------------------------------------#
# Some basic functions
#----------------------------------------------------#
Xfun = function(t, zeta, Z){
  # Input : 
  # t            - knobs of the axis
  # zeta         - vector, zeta = [zeta_1, zeta_2, ..., zeta_K]
  # Z            - vector, Z = [Z_1, Z_2, ..., Z_K]
  # Output: 
  # X            - X := X(t), t in [0, 1], covariate X, functional data 
  
  NumT = length(t)
  K = length(zeta)
  Xvec = numeric(NumT)
  for(i in 1:NumT){
    Xvec[i] = zeta[1] * Z[1]
    for(k in 2:K){
      Xvec[i] = Xvec[i] + sqrt(2) * zeta[k] * 
        Z[k] * cos(k * pi * t[i])
    }
  }
  return(Xvec)
}

XYdata = function(n = 100, K = 50, NumT = 1000, nv = 1.1){
  # Input :
  # n          - sample size
  # K          - order of nonzero coefficient of X(t)
  # NumT       - number of knobs of X = (X(t_1), ..., X(t_NumT))
  # Output:
  # XYdata     - list of Xmat and Y, where Xmat = [X1, X2, ..., Xn]
  #              Y = [Y1, Y2, ..., Yn] 
  
  # zeta[k]^2 are eignvalues
  zeta = numeric(K)
  for(k in 1:K){zeta[k] = (-1)^(k+1) * k^(-nv / 2)}
  
  # Z matrix, Z[k, ] = (Z1, Z2, ..., ZK)^T
  Zmat = matrix(runif(n * K, -sqrt(3), sqrt(3)), nrow = K, ncol = n)
  
  # X matrix, Xmat = [X1(t), X2(t), ..., Xn(t)]
  t = seq(from = 0, to = 1, by = 1 / (NumT - 1))
  Xmat = matrix(NA, nrow = NumT, ncol = n)
  for(i in 1:n){Xmat[, i] = Xfun(t, zeta, Zmat[, i])}
  
  # beta := beta(t)
  beta = numeric(NumT)
  for(i in 1:NumT){
    beta[i] = 0.3
    for(k in 2:K){
      beta[i] = (beta[i] + 4 * sqrt(2) * 
        (-1)^(k + 1) / k^2 * cos(k * pi * t[i]))
    }
  }
  
  # Y vector
  Y = t(Xmat) %*% beta / (NumT - 1) + rnorm(n, 0, 0.5)
  return(list(Xmat, Y))
}

Yhat = function(n, K, NumT, d1hat){
  
  data = XYdata(n, K, NumT)
  X = data[[1]]
  Y = data[[2]]
  
  t = seq(from = 0, to = 1, by = 1 / (NumT -1))
  s = seq(from = 0, to = 1, by = 1 / (NumT -1))
  
  # Matrix M
  Mmat = matrix(c(rep(0.5, n), t(X) %*% matrix(c(rep(1, NumT)), nrow = NumT, ncol = 1) / (NumT - 1)), nrow = n, ncol = 2)
  
  # Function J(x)
  J = numeric(NumT)
  for(x in 1:NumT){
    J[x] = t[x]^2 * log(t[x])
  }
  
  # Matrix Sigma
  Sigmamat = matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      for(p in 1:NumT){
        for(q in 1:NumT){
          Sigmamat[i, j] = 0
          Sigmamat[i, j] = Sigmamat[i, j] + J[p] * sqrt((t[p] - s[q])^2 + (X[p, i] - X[q, j])^2) * (1 / NumT)^2
        }
      }
    }
  }
  
  # Vector Ytilde
  Ytilde = Y - d1hat
  
  #lambda
  lambda = n^(-2 / 3)
  
  # Vector Thetahat
  Thetahat = solve(rbind(cbind(n * lambda * diag(n) + Sigmamat, Mmat), cbind(t(Mmat) %*% Sigmamat, t(Mmat) %*% Mmat))) %*% rbind(Ytilde, t(Mmat) %*% Ytilde)
  
  # Vector Yhat
  Yhat = cbind(Sigmamat, Mmat) %*% Thetahat + d1hat
  
  return(Yhat)
}

n = 20
K = 50
NumT = 100
d1hat = 0.000048

Yhat1 = Yhat(n, K, NumT, d1hat)
