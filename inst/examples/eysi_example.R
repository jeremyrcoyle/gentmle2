######
# Example of TMLE for a stochastic intervention mean E[Y_gstar]

Qbar0 <- function(A, W) {
    
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    Qbar <- plogis(ifelse(W4 > 0, (A == 1) + (A == 1) * (5 * W1^2 - 4.45), (A == 
        2) + (A == 3) + (A == 2) * (4 * W2) + (A == 3) * (5 * W3)))
    return(Qbar)
}

g0 <- function(W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    
    # rep(0.5, nrow(W))
    A1 <- plogis(W1)
    A2 <- plogis(W2)
    A3 <- plogis(W3)
    A <- cbind(A1, A2, A3)
    
    # make sure A sums to 1
    A <- normalize_rows(A)
}

gen_data <- function(n = 1000, p = 4) {
    W <- matrix(rnorm(n * p), nrow = n)
    colnames(W) <- paste("W", seq_len(p), sep = "")
    pA <- g0(W)
    A <- factor(apply(pA, 1, function(pAi) which(rmultinom(1, 1, pAi) == 1)))
    A_vals <- vals_from_factor(A)
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Q0aW <- sapply(A_vals, Qbar0, W)
    d0 <- apply(Q0aW, 1, which.max)
    Yd0 <- as.numeric(u < Qbar0(d0, W))
    data.frame(W, A, Y, Q0aW, d0, Yd0)
}

data <- gen_data(1000, 5)

Anode <- "A"
Wnodes <- grep("^W", names(data), value = T)

Q_fit <- glm(data$Y ~ ., data[, c("A", Wnodes)], family = binomial(link = "logit"))
g_fit <- multinomial_SuperLearner(data$A, data[, Wnodes])

A_vals <- vals_from_factor(data$A)

Q_a <- sapply(A_vals, function(A_val) {
    newdata <- data[, c(Anode, Wnodes)]
    newdata[, Anode] <- A_val
    predict(Q_fit, newdata, type = "response")
})

pA <- predict(g_fit, newdata = data[, Wnodes])$pred

# A sample gstar--treat with A=1, if the patient received A=1 or 2, otherwise
# leave alone
gstar <- function(A, gk) {
    ifelse(A == 1, gk[, 1] + gk[, 2], ifelse(A == 3, gk[, 3], 0))
}

# Initial data set-up
initdata <- data
initdata$Q_a <- Q_a
initdata$pA <- pA

# estimate using Qn, gn
result <- gentmle(initdata, eysi_estimate, eysi_update, max_iter = 100, gstar = gstar)
print(result)

# estimate using Q0, g0
initdata <- data
initdata$Q_a <- sapply(A_vals,Qbar0,data[,Wnodes])
initdata2 <- initdata
initdata2$pA <- g0(data[, Wnodes])
result2 <- gentmle(initdata2, eysi_estimate, eysi_update, max_iter = 100, gstar = gstar)
print(result2) 
