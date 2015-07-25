###### Example of TMLE for the treatment-specific mean E[Y_1]

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

data <- gen_data(1000)
Wnodes <- grep("^W", names(data), value = T)
gk <- g0(data[, Wnodes])
Qk <- Qbar0(data$A, data[, Wnodes])
Q1k <- Qbar0(1, data[, Wnodes])
Q0k <- Qbar0(1, data[, Wnodes])
tmledata <- data.frame(A = data$A, Y = data$Y, gk = gk, Qk = Q1k)
result <- gentmle(tmledata, ey1_estimate, ey1_update)
print(result)

tmledata2 <- data.frame(A = data$A, Y = data$Y, gk = gk, QAk = Qk, Q1k = Q1k, Q0k = Q0k)
result2 <- gentmle(tmledata2, ey1_estimate2, ey1_update2)
print(result2) 
