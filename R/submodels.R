# logistic model/log likelihood loss

#' @export
submodel_logit <- quote(plogis(qlogis(Qk) + HA %*% eps))

#' @export
loss_loglik <- quote(-1 * (Y * log(Qk) + (1 - Y) * log(1 - Qk))) 
