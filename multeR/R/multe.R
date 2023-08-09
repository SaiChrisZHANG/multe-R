#' Multiple Treatment Effects Regression
#'
#' Calculate estimators adjusted for the contamination bias in multiple
#' treatment setting with heterogeneous treatment effects.
#'
#'
#' @template multeFormula
#' @param vce Specifies the type of standard errors to report, it can be a
#' string equal to
#' \describe{
#'    \item{"robust"}{Heteroskedasticity-robust standard errors (default).}
#'
#'    \item{"oracle"}{Heteroskedasticity-robust standard errors, treating
#'    propensity scores as known.}
#'
#' }
#' @param base_val Specifies the baseline value of the treatment variable,
#' the first unique value of the treatment variable is used as default.
#' @param decomp Compute contamination bias decomposition.
#' @param minmax Save saturated group-specific treatment effects as tau
#' and/or implicit ATE regression weights as lambda; computes bias
#' decomposition internally.
#' @param alpha Determines the confidence level, \code{1-alpha} for
#' constructing/optimizing confidence intervals.
#' @param save_lambda_as Saves the set of implicit ATE regression weights,
#' with the input as the prefix of variable names. The saved data frame
#' can be merged back to the data set via index.
#' @param save_tau_as Saves the saturated group-specific treatment effects,
#' with the input as the prefix of variable names. The saved data frame
#' can be merged back to the data set via index.
#' @param print print a summary of the results.
#' @return Returns a summary of the results, containing
#'    \describe{
#'      \item{\code{estimation}}{For treatment effect estimation results, where
#'         \describe{
#'          \item{\code{est}}{Point estimates of the treatment effect, via 3
#'               methods: ATE (\code{est[,1]}), one-at-a-time (\code{est[,2]}),
#'              and common weights (\code{est[,3]}).}
#'
#'          \item{\code{se_po}}{Heteroskedasticity-robust standard errors.}
#'
#'          \item{\code{se_or}}{Oracle standard errors: heteroskedasticity-robust,
#'              treating propensity scores as known.}
#'
#'          \item{\code{po_vcov}}{Variance-covariance matrix, heterogeneity-robust.}
#'
#'          \item{\code{or_vcov}}{Variance-covariance matrix, oracle.}
#'}}
#'      \item{\code{decomposition}}{For contamination bias decomposition
#'        results, where
#'        \describe{
#'          \item{\code{est}}{Point estimates of the contamination bias decomposition,
#'              including: coefficients (\code{est[,1]}), own effect
#'              (\code{est[,2]}), contamination bias (\code{est[,3]}), and its
#'              maximum (\code{est[,4]}) and minimum (\code{est[,5]}).}
#'
#'          \item{\code{se}}{Standard errors.}
#'
#'          \item{\code{lambda_saved}}{Saved lambdas, the set of implicit ATE
#'              regression weights. Only appear when \code{save_lambda_as}
#'              is specified.}
#'
#'          \item{\code{tauhat_saved}}{Saved lambdas, the saturated
#'              group-specific treatment effects. Only appear when
#'              \code{save_tau_as} is specified.}
#' }}
#'}
#'    And in both \code{estimation} and \code{decomposition},
#'     \describe{
#'   \item{\code{n_obs}}{Number of effect observations used for analysis (with
#'              weak overlapping strata dropped.}
#'
#'   \item{\code{n_trt}}{Number of treatment arms.}
#'
#'   \item{\code{Tlevels}}{The label vector treatment arms.}
#'
#' }
#'    If \code{decomp} is not specified, only returns \code{estimation} as a
#'    list.
#' @references{
#'
#' \cite{Goldsmith-Pinkham, Paul, Peter Hull, and Michal Kolesár. Contamination
#' bias in linear regressions. No. w30108. National Bureau of Economic
#' Research, 2022.
#' \doi{10.3386/w30108}}
#'
#' }
#' @examples
#'
#' # Project STAR dataset
#' multe(score ~ treatment | school, data = star)
#' multe(score ~ treatment | school, data = star, decomp = TRUE, minmax = TRUE)
#' multe(score ~ treatment | school, data = star, decomp = TRUE, save_lambda_as = "lambda")
#'
#' @export
multe <- function(formula, data, subset, na.action,
                  vce="robust", base_val,
                  decomp=FALSE, minmax=FALSE, alpha=0.05,
                  save_lambda_as, save_tau_as, print=TRUE) {

  ## argument logics: minmax, save_lambda_as, save_tau_as all conditional on decomp=TRUE
  if (!decomp){
    if (minmax) {
      stop("Must set decomp=TRUE as well to set minmax=TRUE")
    }
    if (!missing(save_lambda_as)) {
      stop("Must set decomp=TRUE as well to save lambda")
    }
    if (!missing(save_tau_as)) {
      stop("Must set decomp=TRUE as well to save tau")
    }
  }

  ## construct model frame
  cl <- mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  formula <- Formula::as.Formula(formula)

  ## one dependent variable (LHS), one treatment assignment indicator (RHS),
  ## and one group of controls to determine the strata
  stopifnot(length(formula)[1] == 1L,
            "Please specify control variables!" =
              length(formula)[2] == 2)
  mf$formula <- formula

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # load data
  d <- multeData(mf,base_val)
  # main estimation
  res <- multeEst.main(d,base_val)

  # contamination bias decomposition
  if (decomp) {
    res_decomp <- multeDecomp.main(d,save_lambda_as,save_tau_as)
  }

  ## print results
  if (print) {
    if (decomp) {
      multeRes.print(res,res_decomp,minmax,vce,alpha)
    }
    else {
      multeRes.print(res,minmax=minmax,vce=vce,alpha=alpha)
    }
  }

  ## save results
  if (decomp) {
    ret <- list(estimation=res,decomposition=res_decomp)
  }
  else {
    ret <- res
  }

  invisible(ret)
}
