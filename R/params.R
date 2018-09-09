
# E_W[E[Y|A=1,W]]
#' @export
param_EY1 <- define_param(psi = "mean(Q1k)", HA = "A/gk", IC = "HA*(Y-Qk)+Q1k-psi")

# E_W[E[Y|A=1,W]-E[Y|A=0,W]]
#' @export
param_ATE <- define_param(psi = "mean(Q1k-Q0k)", HA = "(A/gk - (1 - A)/(1 - gk))", IC = "HA*(Y-Qk)+Q1k-Q0k-psi")

# V_W[E[Y|A=1,W]-E[Y|A=0,W]]
#' @export
param_sigmaATE <- define_param(psi = "mean((Q1k-Q0k-mean(Q1k-Q0k))^2)", HA = "2 * (Q1k - Q0k - mean(Q1k-Q0k)) * (A/gk - (1 - A)/(1 - gk))",
    IC = "HA*(Y-Qk)+(Q1k - Q0k - mean(Q1k-Q0k))^2 - psi")

#' @export
param_ATT <- define_param(psi = "mean((Q1k-Q0k)*gk/mean(A))", HA = "A/mean(A) - (1-A) * gk / (mean(A) * (1-gk))",
                          CA = "(Q1k - Q0k - psi)/mean(A)", IC = "HA*(Y-Qk) + CA*(A - gk) + (Q1k - Q0k - psi)*gk/mean(A)")

#' @export
param_EB2 <- define_param(psi = "mean((Q1k-Q0k)^2)", HA = "2 * (Q1k - Q0k) * (A/gk - (1 - A)/(1 - gk))",
                          IC = "HA*(Y-Qk)+(Q1k - Q0k)^2 - psi")



