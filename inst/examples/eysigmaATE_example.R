###### Example of TMLE for variance of conditional average treatment effect or blip variance
###### var(E[Y|A=1, W] - E[Y|A=0, W])

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
    data.frame(W, A, Y)
}

data <- gen_data(1000)
Wnodes <- grep("^W", names(data), value = T)
gk <- g0(data[, Wnodes])
Qk <- Qbar0(data$A, data[, Wnodes])
Q1k <- Qbar0(1, data[, Wnodes])
Q0k <- Qbar0(0, data[, Wnodes])

# notice here we specify two parameters to be simultaneously estimated
initdata <- data.frame(A = data$A, Y = data$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "recursive", max_iter = 10000)
print(result)

# for iterative TMLE, choose full--different approaches sometimes give slightly
# different answers
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "full")
print(result)

# One can also form simultaneous confidence bounds for numerous params using the
# influence curves by specifying simultaneous.inference = TRUE
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "recursive", simultaneous.inference = TRUE)
print(result)


