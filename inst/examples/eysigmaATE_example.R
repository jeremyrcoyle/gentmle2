###### Example of TMLE for simultaneous estimation of two parameters,
###### ATE and blip variance

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


