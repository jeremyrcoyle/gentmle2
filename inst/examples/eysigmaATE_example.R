gendata <- function(n) {
    W1 <- rnorm(n)
    W2 <- rnorm(n)
    W3 <- rnorm(n)
    W4 <- rnorm(n)
    A <- rbinom(n, 1, g0(W1, W2, W3, W4))
    Y <- rbinom(n, 1, Q0(A, W1, W2, W3, W4))
    data.frame(W1, W2, W3, W4, A, Y)
}

g0 <- function(W1, W2, W3, W4) {
    plogis(-0.28 * W1 + 0.5 * W2 + 0.08 * W3 - 2 * W4)
}
Q0 <- function(A, W1, W2, W3, W4) {
    plogis(4 * A + 20 * A * W2 - 50 * W4^2 - 10 * W3^2 + W1 * W2)
}

initdata <- gendata(1000)
data1 <- initdata
data1$A <- 1
data0 <- initdata
data0$A <- 0

Qfit <- glm(Y ~ A + W1 + W2 + W3 + W4, data = initdata, family = binomial(link = "logit"))
gfit <- glm(A ~ W1 + W2 + W3 + W4, data = initdata, family = binomial(link = "logit"))

initdata$gk <- predict(gfit, type = "response")
initdata$Qk <- predict(Qfit, type = "response")
initdata$Q1k <- predict(Qfit, newdata = data1, type = "response")
initdata$Q0k <- predict(Qfit, newdata = data0, type = "response")

results <- gentmle(initdata, eysigmaATE_estimate, eysigmaATE_update, max_iter = 100)

CI <- ci_gentmle(results, level = 0.95)
CI 
