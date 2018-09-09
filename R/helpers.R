# truncation function for Q so logistic regression doesn't break on Y close to 0 or 1
#' @export
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
    pmin(pmax(x, lower), upper)
}

# evaluate parameter specification using data from an environment
eval_param <- function(param, tmleenv) {
    res <- list(psi = NULL, HA = NULL, H1 = NULL, H0 = NULL, IC = NULL)

    # this ordering allows IC to be a function of the other elements
    res$psi <- eval(param$psi, tmleenv)
    res$HA <- eval(param$HA, tmleenv)
    res$H1 <- eval(param$HA, list(A = 1), tmleenv)
    res$H0 <- eval(param$HA, list(A = 0), tmleenv)
    res$IC <- eval(param$IC, res, tmleenv)

    return(res)
}

# evaluate risk for a given Qk
risk_Qk <- function(Qk, tmleenv, loss, HA) {
    mean(eval(loss, list(Qk = Qk, HA = HA), tmleenv))
}

# fluctuate on a submodel and evaluate risk for a given epsilon
risk_eps <- function(eps, HA, tmleenv, submodel, loss) {
    Qk_eps <- eval(submodel, list(HA = HA, eps = eps), tmleenv)
    risk_Qk(Qk_eps, tmleenv, loss, HA = HA)
}
