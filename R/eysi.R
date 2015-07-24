library(tmle)
library(foreach)
library(doSNOW)
library(origami)
library(glmnet)
library(ggplot2)
library(plyr)
library(dplyr)
registerDoSNOW(makeCluster(4, type = "SOCK"))

Qbar0 <- function(A, W) {
  
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]
  Qbar <- plogis(ifelse(W4 > 0, 1 - W1^2 + 3 * W2 + (A==1) * (5 * W3^2 - 4.45) + (A==2), 
                        -0.5 - W3 + 2 * W1 * W2 + (A==1) * (3 * abs(W2) - 1.5) + 5 * (A==3) * W3))
  return(Qbar)
}

g0 <- function(W) {
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]
  
  # rep(0.5, nrow(W))
  A1 <- plogis(W1 - W2)
  A2 <- plogis(W2 - W1)
  A3 <- plogis(W2 - W3)
  A <- cbind(A1, A2, A3)
  
  # make sure A sums to 1
  A <- t(apply(A, 1, function(x) x/sum(x)))
}

gen_data <- function(n = 1000, p = 4) {
  W <- matrix(rnorm(n * p), nrow = n)
  colnames(W) <- paste("W", seq_len(p), sep = "")
  pA <- g0(W)
  A <- factor(apply(pA, 1, function(pAi) which(rmultinom(1, 1, pAi) == 1)))
  A_vals <- levels(A)
  #A <- t(apply(pA, 1, function(pAi) rmultinom(1, 1, pAi) == 1))
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Ya=sapply(A_vals,Qbar0,W)
  d0 <- apply(Ya,1,which.max)
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  data.frame(W, A, Y,Ya=Ya, d0,Yd0)
}

data <- gen_data(1000, 5)

Anode <- "A"
Wnodes <- grep("^W", names(data), value = T)
# Q_library <- c("SL.glm", "SL.glmem", "SL.glmnetprob", "SL.step.forward", "SL.gam", 
#     "SL.rpart", "SL.rpartPrune", "SL.mean")
Q_library <- c("SL.glm","SL.step.forward", "SL.gam","SL.rpart", "SL.rpartPrune", "SL.mean")
# Q_fit=origami_SuperLearner(data$Y,data[,c("A",Wnodes)],SL.library=Q_library,family=binomial())
Q_fit=glm(data$Y~.,data[,c("A",Wnodes)],family=binomial(link="logit"))
g_fit=multinomial_SuperLearner(data$A,data[,Wnodes])

A_vals=sort(unique(data$A))

Q_a <- sapply(A_vals, function(A_val) {
  newdata <- data[, c(Anode, Wnodes)]
  newdata[, Anode] <- A_val
  predict(Q_fit, newdata,type="response")
})

pA=predict(g_fit,newdata=data[,Wnodes])$pred

# A sample gstar--treat with A=1, if the patient received A=1 or 2, otherwise leave alone
gstar = function(gk,A) {
  ifelse(A==1,gk[1]+gk[2],ifelse(A==3,gk[3],0))
}

# Initial data set-up
initdata=data
initdata$Q_a=Q_a
initdata$pA=pA

# stochastic intervention update for the mean
eysi_update <- function(tmledata, Q.trunc = 0.001) {
  
  subset <- with(tmledata, which(Qk>0 & Qk<1))
  eps_q <- 0
  
  # fluctuate Q
  tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
  qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qktrunc)))
  #     qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(Qk))
  eps_q <- qfluc$eps
  tmledata$Q_a <- with(tmledata, plogis(qlogis(Q_a) + H * eps_q))
  # tmledata$Qk=qfluc$update
  
  list(tmledata = tmledata, coefs = c(eps_q))
  
}

#' @export

# Stochastic intervention estimate for the mean
eysi_estimate <- function(tmledata, ...) {
  
  # assign probs under stochastic intervention
  pAstar = mapply(A_vals,FUN=function(x) apply(tmledata$pA,1,...,x))
  tmledata$pAstar = pAstar
  
  # compute the parameter under stochastic intervention for fixed g
  psi <- mean(apply(tmledata$Q_a*tmledata$pAstar,1,sum))
  
  # which A popped up in reality
  ind = mapply(A_vals,FUN=function(x) as.numeric(tmledata$A==x))
  
  # clever coordinate for each treatment
  tmledata$H = with(tmledata,pAstar/pA)
  
  # clever coordinate for the fitted QAW
  tmledata$HA = apply(ind*tmledata$H,1,sum)
  
  # create QAW--our predictions over the treatments that occurred 
  tmledata$Qk = apply(ind*tmledata$Q_a,1,sum)
  
  # influence curves
  tmledata$empirical = apply(with(tmledata,pAstar*Q_a),1,sum)
  Dstar_psi <- with(tmledata, HA * (Y - Qk) + empirical - psi)
  
  list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
  
}

testing = gentmle(initdata,eysi_estimate,eysi_update,max_iter=100,gstar)
