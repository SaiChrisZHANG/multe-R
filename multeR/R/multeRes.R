## Print results

multeRes.print <- function(res,res_decomp,minmax,vce,alpha) {

  ## print estimation results
  cat("Treatment Effect Estimates (with ",vce, " standard errors):\n", sep="")
  cat("ATE\n")
  print(multeRes.EsttabConst(res,1,vce,alpha),quote=F,right=T)
  cat("One-at-a-time\n")
  print(multeRes.EsttabConst(res,2,vce,alpha),quote=F,right=T)
  cat("Common weights\n")
  print(multeRes.EsttabConst(res,3,vce,alpha),quote=F,right=T)

  if (!missing(res_decomp)) {
    cat("\n")
    cat("Contamination Bias Decomposition:\n")
    outtab <- multeRes.EsttabDecomp(res_decomp,minmax)
    print(outtab,quote=F,right=T)
  }
}


## construct the results table for treatment effect estimation, to be aligned with the Stata output
multeRes.EsttabConst <- function(res,n,vce,alpha) {
  ## column names
  colnm_or <- c("Coef.","Oracle S.E.","z","P>|z|",
                paste0((1-alpha)*100,"% Confidence Interval"))
  colnm_po <- c("Coef.","Robust S.E.","z","P>|z|",
                paste0((1-alpha)*100,"% Confidence Interval"))

  if (vce=="oracle"){
    outtab <-cbind(sprintf("%.7f",res$est[,n]),
                   sprintf("%.7f",res$se_or[,n]),
                   sprintf("%.2f",res$est[,n]/res$se_or[,n]),
                   sprintf("%.3f",2*pnorm(q=abs(res$est[,n]/res$se_or[,n]),lower.tail=FALSE)),
                   paste0("(",
                          sprintf("%.7f",res$est[,n] - stats::qnorm(1-alpha/2)*res$se_or[,n]),
                          ", ",
                          sprintf("%.7f",res$est[,n] + stats::qnorm(1-alpha/2)*res$se_or[,n]),
                          ")")
    )
    colnames(outtab) <- colnm_or
    rownames(outtab) <- res$Tlevels[-1]
  }
  else {
    outtab <- cbind(sprintf("%.7f",res$est[,n]),
                    sprintf("%.7f",res$se_po[,n]),
                    sprintf("%.2f",res$est[,n]/res$se_po[,n]),
                    sprintf("%.3f",2*pnorm(q=abs(res$est[,n]/res$se_po[,n]),lower.tail=FALSE)),
                    paste0("(",
                           sprintf("%.7f",res$est[,n] - stats::qnorm(1-alpha/2)*res$se_po[,n]),
                           ", ",
                           sprintf("%.7f",res$est[,n] + stats::qnorm(1-alpha/2)*res$se_po[,n]),
                           ")")
    )
    colnames(outtab) <- colnm_po
    rownames(outtab) <- res$Tlevels[-1]
  }
  outtab
}

## construct the results table for contamination bias decomposition
multeRes.EsttabDecomp <- function(res_decomp,minmax) {
  outest <- matrix(sprintf("%.7f",res_decomp$est),ncol=5)
  rownames(outest) <- res_decomp$Tlevels[-1]
  outse <- matrix(paste0("(",sprintf("%.7f",res_decomp$se),")"),ncol=5)
  rownames(outse) <- rep(NULL,length(res_decomp$Tlevels[-1]))
  outtab <- rbind(outest,outse)[rep(seq_len(nrow(outest)),each=2)+c(0,nrow(outest)),]

  colnames(outtab) <- c("Coef.","Own Effect","Bias","Min Bias","Max Bias")

  if (minmax) {
    outtab
  }
  else {
    outtab[,c("Coef.","Own Effect","Bias")]
  }
}
