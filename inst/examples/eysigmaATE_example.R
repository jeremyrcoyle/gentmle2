###### Example of TMLE for simultaneous estimation of two parameters,
###### ATE and blip variance

# notice here we specify two parameters to be simultaneously estimated
initdata <- data.frame(A = data$A, Y = data$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "recursive", max_iter = 10000)
print(result)

# for iterative TMLE, choose full--different approaches sometimes give slightly
# different answers

# full approach computes a separate epsilon for each parameter
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "full")
print(result)

# line approach is the iterative analog of the 1 step TMLE
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "line")
print(result)

# recursive is the 1 step tmle as in Mark van der Laan's recent work
result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "recursive", max_iter = 10000)

# One can also form simultaneous confidence bounds for numerous params using the
# influence curves by specifying simultaneous.inference = TRUE

result <- gentmle(initdata = initdata, params = list(param_ATE, param_sigmaATE),
                  approach = "recursive", max_iter = 10000, simultaneous.inference = TRUE)
print(result)


