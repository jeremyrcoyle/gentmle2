# logistic model/log likelihood loss

#' @export
submodel_logit <- quote(plogis(qlogis(Qk) + HA %*% eps))

#' @export
loss_loglik <- quote(-1 * (Y * log(Qk) + (1 - Y) * log(1 - Qk)))

#' @export
submodel_logit_intercept <- quote(plogis(qlogis(Qk) +
                                           as.matrix((HA<0)*-1 + (HA>0) + (HA==0)*0) %*% eps))

#' @export
loss_loglik_wts <- quote(abs(HA) * (-1 * (Y * log(Qk) + (1 - Y) * log(1 - Qk))))




