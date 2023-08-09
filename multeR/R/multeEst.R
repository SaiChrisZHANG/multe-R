## estimation functions called by multe.R

## main estimation function
multeEst.main <- function(d,base_val) {

  # set up ===============
  ## preparing design matrices
  Tm <- d$Tm # treatment, baseline value as first column
  Tlevels <- d$Tlevels # levels of treatment
  Cm <- d$Cm # control groups

  n <- length(d$Y)
  k <- length(col(Tm)[1,])
  T0 <- Tm[,1]
  alpha0 <- multeEst.ols(d$Y[T0==1],Cm[T0==1,])
  psi_alpha0 <- (T0 * Cm * as.vector(d$Y - Cm %*% alpha0)) %*%
    qr.solve(crossprod(T0 * Cm, Cm)/n)

  est = se_or = se_po = matrix(NA, nrow = k-1, ncol = 3)

  ## common weights
  ps <- Cm %*% multeEst.ols(Tm,Cm) # propensity score
  Tt <- Tm - ps # residuals
  lam <- 1 / rowSums(1/ps)
  rcalpha <- multeEst.olsw(d$Y,Tm,lam/rowSums(ps*Tm))
  rcres <- d$Y - Tm %*% rcalpha
  est[,3] <- rcalpha[2:length(rcalpha)] - rcalpha[1]

  ## fitted t values
  ts <- Cm %*% multeEst.ols(as.vector(rcres)*Tm/(ps^2)*lam, Cm)
  s0 <- Cm %*% multeEst.ols(as.vector(rcres)*T0/ps[,1]*(lam^2)/(ps^2), Cm)

  # TE estimation ===============
  var_po_onem = var_or_onem = matrix(0, nrow = k, ncol = 1)
  psi_pom = psi_orm = matrix(0, nrow = n, ncol = k*3)
  Cmean <- as.vector(colMeans(Cm))

  ## looping through treatment level
  for (j in 2:k){
    alphak <- multeEst.ols(d$Y[Tm[,j]==1],Cm[Tm[,j]==1,])
    psi_alphak <- ( Tm[,j] * as.vector(d$Y- Cm%*%alphak) * Cm) %*%
      qr.solve(crossprod(Tm[,j]*Cm, Cm)/n)

    ### estimation 1: ATE
    ### WARNING: here psi_pom are numerically off a bit
    est[j-1,1] <- Cmean %*% (alphak - alpha0)
    psi_or <- (psi_alphak - psi_alpha0) %*% Cmean
    psi_orm[,j] <- psi_or
    psi_pom[,j] <- psi_or + ((Cm-Cmean)%*%(alphak-alpha0))
    se_or[j-1,1] <- sqrt(var(psi_or) * (n-1)/n^2)
    se_po[j-1,1] <- sqrt(var(psi_pom[,j]) * (n-1)/n^2)

    ### estimation 2: One treatment at a time
    s <- which(T0==1|Tm[,j]==1) # select one treatment and control
    Tm_s <- Tm[s,j]
    Cm_s <- Cm[s,]
    T0_s <- T0[s]
    Y_s  <- d$Y[s]
    Xdot <- Tm_s - Cm_s %*% multeEst.ols(Tm_s, Cm_s)
    rk <- multeEst.olsres(Y_s, cbind(Xdot, Cm_s))
    est[j-1,2] <- rk$coef[1]
    se_po[j-1,2] <- sqrt(sum((rk$res^2)*(Xdot^2)) / sum(Xdot^2)^2)

    eps <- Tm_s * (Y_s - Cm_s %*% alphak) + T0_s * (Y_s - Cm_s %*% alpha0)
    se_or[j-1,2] <- sqrt(sum((eps^2)*(Xdot^2))/sum(Xdot^2)^2)
    var_po_onem[j,1] <- sum((rk$res^2) * (Xdot^2)) / sum(Xdot^2)^2
    var_or_onem[j,1] <- sum((eps^2) * (Xdot^2)) / sum(Xdot^2)^2

    ### estimation 3: common weight
    psi_or <- lam * (Tm[,j] * (d$Y - Cm %*% alphak) / ps[,j] -
                       T0 * (d$Y - Cm %*% alpha0) / ps[,1]) / mean(lam)
    sk <- Cm %*% multeEst.ols(as.vector(rcres * Tm[,j] / ps[,j] * (lam^2)
                                        ) / (ps^2), Cm)
    psi_po <- (lam*rcres*(Tm[,j]/ps[,j] - T0/ps[,1]) + Tt[,1] * ts[,1]
               - Tt[,j] * ts[,j] + rowSums(Tt * (sk - s0))) / mean(lam)
    se_po[j-1,3] <- sqrt(var(psi_po)*(n-1)/n^2)
    se_or[j-1,3] <- sqrt(var(psi_or)*(n-1)/n^2)
    psi_pom[,j+(k*2)] <- psi_po
    psi_orm[,j+(k*2)] <- psi_or
  }

  # compute variance-covariance matrix ===========
  ### NOTE: Var(beta) = Var(psi)/n. That's why po_vcov has n^2 in the denominator and not n
  psi_pom_tl <- sweep(psi_pom,2,(colSums(psi_pom) / n))
  psi_orm_tl <- sweep(psi_orm,2,(colSums(psi_orm) / n))
  po_vcov <- crossprod(psi_pom_tl) / (n^2)
  or_vcov <- crossprod(psi_orm_tl) / (n^2)

  for (j in 2:k) {
    po_vcov[j+k,j+k] = var_po_onem[j,1]
    or_vcov[j+k,j+k] = var_or_onem[j,1]
  }

  ## build a block matrix out of po_vcov and or_vcov
  po_vcov <- as.matrix(Matrix::bdiag((po_vcov[1:k,1:k]),
                                     (po_vcov[(k+1):(k*2),(k+1):(k*2)]),
                                     (po_vcov[(2*k+1):(k*3),(2*k+1):(k*3)])))
  or_vcov <- as.matrix(Matrix::bdiag((or_vcov[1:k,1:k]),
                                     (or_vcov[(k+1):(k*2),(k+1):(k*2)]),
                                     (or_vcov[(2*k+1):(k*3),(2*k+1):(k*3)])))


  # store estimation results in a list ============
  est_result <- list(est = est, se_po = se_po, se_or = se_or,
                     n_obs = n, n_trt = k, Tlevels = Tlevels,
                     po_vcov = po_vcov, or_vcov= or_vcov)
  est_result
}


## auxiliary functions
multeEst.ols <- function(Y,X) {
  ## OLS estimation
  ols_est <- solve(crossprod(X, X)) %*% crossprod(X, Y)
  ols_est
}

multeEst.olsw <- function(Y,X,W) {
  ## weighted OLS estimation
  olsw_est <- qr.solve(crossprod(X*W,X)) %*% crossprod(X*W, Y)
  olsw_est
}

multeEst.olsres <- function(Y,X) {
  ## OLS residuals
  coef <- multeEst.ols(Y,X)
  res <- Y - X%*%coef
  ols_res <- list(coef=coef, res=res)
  ols_res
}

multeEst.olswres <- function(Y,X,W) {
  ## OLS residuals
  coef <- multeEst.olsw(Y,X,W)
  res <- Y - X%*%coef
  ols_res <- list(coef=coef, res=res)
  ols_res
}
