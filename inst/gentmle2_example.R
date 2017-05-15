library(devtools)
install_github("jeremyrcoyle/gentmle")

library(gentmle2)
################################## generate some data
Qbar0 <- function(A, W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    Qbar <- plogis(A + A * W1 + W2)
    return(Qbar)
}

g0 <- function(W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    # rep(0.5, nrow(W))
    plogis(0.25 * W1 - 0.1 * W2)
}

gen_data <- function(n = 1000, p = 2) {
    W <- matrix(rnorm(n * p), nrow = n)
    colnames(W) <- paste("W", seq_len(p), sep = "")
    A <- rbinom(n, 1, g0(W))
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Y0 <- as.numeric(u < Qbar0(0, W))
    Y1 <- as.numeric(u < Qbar0(1, W))
    data.frame(W, A, Y, Y0, Y1)
}

n <- 1000
data <- gen_data(n)
Wnodes <- grep("^W", names(data), value = T)

test_data <- gen_data(1e+05)
test_data$Q1k <- Qbar0(1, test_data[, Wnodes])
test_data$Q0k <- Qbar0(0, test_data[, Wnodes])
psi0_sigma <- eval(param_sigmaATE$psi, test_data)
psi0_ATE <- eval(param_ATE$psi, test_data)
gk <- g0(data[, Wnodes])
Qk <- Qbar0(data$A, data[, Wnodes])
Q1k <- Qbar0(1, data[, Wnodes])
Q0k <- Qbar0(0, data[, Wnodes])
tmledata <- data.frame(A = data$A, Y = data$Y, gk = gk, Q1k = Q1k, Q0k = Q0k, Qk = Qk)

# tmledata with crappy Q predictions (to make things at least a little interesting)

bad_tmledata <- data.frame(A = data$A, Y = data$Y, gk = gk, Q1k = runif(n), Q0k = runif(n), Qk = mean(data$Y))

####################################################### check out Jonathan's code
source("~/Downloads/onestep_blipvar.R")
Q <- bad_tmledata[, c("Qk", "Q0k", "Q1k")]
names(Q) <- c("QAW", "Q0W", "Q1W")
Q <- as.matrix(Q)
bounds <- c(0, 1)
stiles_res <- onestep_atesig(bad_tmledata$Y, bad_tmledata$A, Q, g1W = bad_tmledata$gk, depsilon = 1e-04, max_iter = 10000, 
    gbounds = bounds, Qbounds = bounds)
stiles_res$psi
stiles_res$ICmeans
mean(eval(loss_loglik, list(Y = bad_tmledata$Y, Qk = stiles_res$Q[, "QAW"])))

# recursive/one-step TMLE (should produce the same results -- it's close)
params <- list(ATE = param_ATE, sigma = param_sigmaATE)
res <- gentmle(bad_tmledata, params = params, max_iter = 10000, approach = "recursive")
res$initests
res$tmleests
res$ED
mean(eval(loss_loglik, res$tmledata))
res$steps

# full search TMLE (kinda like fitting a GLM)
res <- gentmle(bad_tmledata, params = params, max_iter = 100, approach = "full")
res$tmleests
res$ED
mean(eval(loss_loglik, res$tmledata))
res$steps

# one dimensional line search TMLE (somewhere in between)
res <- gentmle(bad_tmledata, params = params, approach = "line", max_iter = 100)
res$tmleests
res$ED
mean(eval(loss_loglik, res$tmledata))
res$steps 
