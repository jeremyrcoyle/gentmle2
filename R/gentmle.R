library(boot)
############## somewhat general TMLE framework takes initial data, and a tmle parameter spec todo:
############## reimplement truncation (can just be a separate submodel) implement gradient for
############## optimization approaches handle parameters which require fluctuating g

#' @title gentmle
#' @description General TMLE function that takes care of the bookkeeping of estimation and update steps.
#'
#' @param paramlist, named list of parameters to estimate. See define_param for details
#' @param submodel, submodel along which to fluctuate
#' @param loss, loss function to optimize
#' @param approach, One of initial, recursive (small delta), line, full
#' @param max_iter, Maximum number of iteration steps
#'
#' @export
##' @example /inst/examples/ey1_example.R
##' @example /inst/examples/eysi_example.R
gentmle <- function(initdata, params, submodel = submodel_logit, loss = loss_loglik,
                    approach = "recursive", max_iter = 100, ...) {
  tmleenv <- list2env(initdata)
  converge <- F
  last_risk <- Inf

  initests <- NULL
  n <- length(tmledata$Y)
  for (j in seq_len(max_iter)) {
    # print('Qks') print(mean(tmleenv$Qk)) print(mean(tmleenv$Q1k))
    # print(mean(tmleenv$Q0k))

    # calculate parameters and associated values
    evals <- lapply(params, eval_param, tmleenv)
    Dstar <- sapply(evals, `[[`, "IC")
    ED <- apply(Dstar, 2, mean)
    sig <- sapply(evals, function(param) sd(param$IC))
    psi <- sapply(evals, `[[`, "psi")
    if (j == 1) {
      initests <- psi
    }
    # get HAs (for each param)
    HAs <- sapply(evals, `[[`, "HA")
    H1s <- sapply(evals, `[[`, "H1")
    H0s <- sapply(evals, `[[`, "H0")

    # print('HAs') print(mean(HAs)) print(mean(H0s)) print(mean(H1s))

    if (approach == "full") {
      # optimize across multidimensional HA
      HA <- HAs
      H1 <- H1s
      H0 <- H0s

    } else {
      # combine into single HA
      EDnormed <- ED/norm(ED, type = "2")
      HA <- HAs %*% EDnormed
      H1 <- H1s %*% EDnormed
      H0 <- H0s %*% EDnormed
    }

    if (approach == "initial") {
      eps <- 0
    } else if (approach == "recursive") {
      eps <- 0.0001
    } else if (approach == "line" || approach == "full") {
      init_eps <- rep(0, ncol(HA))
      opt <- optim(par = init_eps, risk_eps, HA = HA, tmleenv = tmleenv, method = "L-BFGS-B")
      eps <- opt$par
    }

    risk <- risk_eps(eps, HA, tmleenv)

    # check for improvement
    if (risk >= last_risk) {
      # we failed to improve, so give up
      j <- j - 1
      converge <- T
      break
    }
    if (all(abs(ED) < sig/n)) {
      # update Qk
      tmleenv$Qk <- as.vector(eval(submodel, list(HA = HA, eps = eps, Qk = tmleenv$Qk)))
      tmleenv$Q1k <- as.vector(eval(submodel, list(HA = H1, eps = eps, Qk = tmleenv$Q1k)))
      tmleenv$Q0k <- as.vector(eval(submodel, list(HA = H0, eps = eps, Qk = tmleenv$Q0k)))
      evals <- lapply(params, eval_param, tmleenv)
      Dstar <- sapply(evals, `[[`, "IC")
      ED <- apply(Dstar, 2, mean)
      psi <- sapply(evals, `[[`, "psi")
      converge <- T
      break
    } else {
      tmleenv$Qk <- as.vector(eval(submodel, list(HA = HA, eps = eps, Qk = tmleenv$Qk)))
      tmleenv$Q1k <- as.vector(eval(submodel, list(HA = H1, eps = eps, Qk = tmleenv$Q1k)))
      tmleenv$Q0k <- as.vector(eval(submodel, list(HA = H0, eps = eps, Qk = tmleenv$Q0k)))
    }

  }

  ED2 <- sapply(evals, function(param) mean(param$IC^2))
  ED3 <- sapply(evals, function(param) mean(param$IC^3))
  result <- list(initdata = initdata, tmledata = data.frame(as.list(tmleenv)), initests = initests,
                 tmleests = psi, steps = j, Dstar = Dstar, ED = ED, ED2 = ED2, ED3 = ED3)

  class(result) <- "gentmle"

  return(result)
}

#' @export
ci_gentmle <- function(gentmle_obj, level = 0.95) {

    n <- nrow(gentmle_obj$initdata)
    n_ests <- length(gentmle_obj$tmleests)
    ldply(seq_len(n_ests), function(i) {

        est <- gentmle_obj$tmleests[i]
        se <- sqrt(gentmle_obj$ED2[i])/sqrt(n)
        z <- (1 + level)/2
        lower <- est - qnorm(z) * se
        upper <- est + qnorm(z) * se
        data.frame(parameter = names(est), est = est, se = se, lower = lower, upper = upper)
    })
}

#' @export
print.gentmle <- function(gentmle_obj) {
    cat(sprintf("TMLE ran for %d step(s)\n", gentmle_obj$steps))
    EDtext <- sprintf("E[%s]=%1.2e", names(gentmle_obj$ED), gentmle_obj$ED)
    cat(sprintf("The mean of the IC is %s\n", paste(EDtext, collapse = ", ")))

    cat("\n\n")
    print(ci_gentmle(gentmle_obj))

    cat("\n")

}
