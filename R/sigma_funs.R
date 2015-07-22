
############## TMLE targeting Sigma^2(Ey1)

#' @export
sigma_update <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # logit_fluctuate Q
        qfluc <- logit_fluctuate(tmledata, Y ~ -1 + CyA + offset(qlogis(Qk)), subset)
        eps_q <- qfluc$eps
        tmledata$Qk <- with(tmledata, plogis(qlogis(Qk) + Cy1 * eps_q))
    }
    
    # logit_fluctuate g tmledata$Ca=with(tmledata,Qk*(1-Qk)/(gk^2))
    gfluc <- logit_fluctuate(tmledata, A ~ -1 + Ca + offset(qlogis(gk)))
    eps_g <- gfluc$eps
    tmledata$gk <- with(tmledata, plogis(qlogis(gk) + Ca * eps_g))
    
    list(tmledata = tmledata, coefs = c(eps_q, eps_g))
    
}

#' @export
sigma_estimate <- function(tmledata) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$Cy1 <- with(tmledata, (1/gk) * ((1 - 2 * Qk)/gk + 2 * (Qk - psi)))
    tmledata$CyA <- with(tmledata, A * Cy1)
    tmledata$Ca <- with(tmledata, Qk * (1 - Qk)/(gk^2))
    
    # influence curve
    varterm <- with(tmledata, Qk * (1 - Qk)/gk + (Qk - psi)^2)
    sigma2 <- mean(varterm)
    D_Qw <- varterm - sigma2
    D_Qbar <- with(tmledata, CyA * (Y - Qk))
    D_gbar <- with(tmledata, -Ca * (A - gk))
    Dstar_sigma <- D_Qw + D_Qbar + D_gbar
    
    
    list(tmledata = tmledata, ests = c(sigma2 = sigma2), Dstar = list(Dstar_sigma = Dstar_sigma))
} 
