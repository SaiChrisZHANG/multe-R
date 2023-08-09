## variance decomposition functions called by multe.R

## main decomposition function
multeDecomp.main <- function(d,save_lambda_as,save_tau_as) {

  Tm <- d$Tm # treatment, baseline value as first column
  Tlevels <- d$Tlevels # levels of treatment
  Cm <- d$Cm # control groups

  # set up =================
  n <- length(d$Y)
  k <- length(col(Tm)[1,])
  kc <- length(col(Cm)[1,])

  ktc <- kc + (k-1) * kc
  M <- kronecker(diag(1,k-1),rep(1,kc))
  Tm <- Tm[,2:k]
  CmTm <- matrix(0,nrow=n,ncol=ktc)

  # elements of CmTm
  CmTm[,1:kc] <- Cm
  for (j in 1:(k-1)){
    col_start <- kc + (j-1) * kc + 1
    col_end <- kc + j * kc
    CmTm[,col_start:col_end] <- Tm[,j]*Cm
  }

  ri <- multeEst.olsres(d$Y, CmTm)
  gamma <- ri$coef[(kc+1):ktc]
  rd <- multeEst.olsres(CmTm[,(kc+1):ktc], cbind(Tm,Cm))
  psi_gamma <- ((as.vector(ri$res) * CmTm) %*% solve(crossprod(CmTm)))[,(kc+1):ktc]

  # sort columns by size
  ghelper <- c(matrix(rep(matrix(1:(k-1),ncol=(k-1)),kc),ncol=(k-1),byrow=TRUE))
  gi <- order(ghelper, gamma) # ascending
  gd <- order(ghelper,-gamma) # descending
  est <- se <- matrix(NA, nrow = k-1, ncol = 5)
  ggi <- gamma[gi] * M
  ggd <- gamma[gd] * M
  psigi <- psi_gamma[,gi]
  psigd <- psi_gamma[,gd]

  # standard errors
  for (j in 1:(k-1)){
    rddX <- multeEst.olsres(Tm[,j], cbind(Tm[,-j], Cm))

    deltak <- rd$coef[j,]
    psi_deltak <- as.vector(rddX$res) * (rd$res) / sum(as.vector(rddX$res)^2)
    psi <- psi_deltak %*% (gamma * M) + psi_gamma %*% (deltak * M)
    estk <- crossprod(gamma, (deltak*M))

    di <- order(ghelper, deltak)

    psi_deltak <- psi_deltak[,di]
    deltak = deltak[di] * M
    psimax = psi_deltak %*% ggi + psigi %*% deltak
    psimin = psi_deltak %*% ggd + psigd %*% deltak

    est[j,] = c(sum(estk), estk[j], sum(estk[-j]),
                sum(crossprod(ggd,deltak[,-j])),
                sum(crossprod(ggi,deltak[,-j]))
                )
    se[j,] = sqrt(c(sum(rowSums(psi)^2), sum(psi[,j]^2),
                    sum(rowSums(psi[,-j,drop=F])^2),
                    sum(rowSums(psimin[,-j,drop=F])^2),
                    sum(rowSums(psimax[,-j,drop=F])^2)
                    )
                  )

  }

  # finalizing output
  tauhat_names <- sprintf("%i",c(1:(k-1)))
  lambda_names <- sprintf("%i",kronecker(c(1:(k-1)),rep(10,(k-1))) + rep(c(1:(k-1)),(k-1)))

  delta_kl <- matrix(as.vector(t(rd$coef[1:(k-1),])),ncol=(k-1)^2)
  delta_pr <- delta_kl / as.vector(colMeans(Cm))
  gammam <- t(matrix(gamma,ncol=(k-1)))

  # estimation results
  result <- list(est = est, se = se,
                 n_obs = n, n_trt = k, Tlevels = Tlevels)

  # also return lambda and tau if specified
  if (!missing(save_lambda_as)){
    lambda_saved <- Cm %*% delta_pr
    colnames(lambda_saved) <- paste0(save_lambda_as, lambda_names)
    rownames(lambda_saved) <- as.numeric(d$index)
    result <- append(result,list(lambda_saved = as.data.frame(lambda_saved)))
  }
  if (!missing(save_tau_as)){
    tauhat_saved <- tcrossprod(Cm,gammam)
    colnames(tauhat_saved) <- paste0(save_tau_as, tauhat_names)
    rownames(tauhat_saved) <- as.numeric(d$index)
    result <- append(result,list(tauhat_saved = as.data.frame(tauhat_saved)))
  }

  result
}
