# logistic model/log likelihood loss

#' @export
Q_submodel_logit <- quote(plogis(qlogis(Qk) + HA %*% eps))

#' @export
Q_loss_loglik <- quote(-1 * (Y * log(Qk) + (1 - Y) * log(1 - Qk)))

#' @export
g_submodel_logit <- quote(plogis(qlogis(gk) + CA %*% eps))


#' @export
g_loss_loglik <- quote(-1 * (A * log(gk) + (1 - A) * log(1 - gk)))

#' @export
Q_submodel_logit_weighted <- quote(plogis(qlogis(Qk) + sign(HA) %*% eps))

#' @export
Q_loss_loglik_weighted <- quote(-abs(HA) * (Y * log(Qk) + (1 - Y) * log(1 - Qk)))
