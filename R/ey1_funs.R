
############## TMLE targeting EY1
#' @export
ey1_update <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # luctuate Q
        qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qk)), subset)
        eps_q <- qfluc$eps
        tmledata$Qk <- with(tmledata, plogis(qlogis(Qk) + H1 * eps_q))
        # tmledata$Qk=qfluc$update
    }
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

#' @export
ey1_estimate <- function(tmledata) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
}

############## TMLE targeting EY1 (more canonical than above, but should be equivalent)

#' @export
ey1_update2 <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < QAk & QAk < 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # logit_fluctuate Q
        qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(QAk)), subset)
        eps_q <- qfluc$eps
        tmledata$Q1k <- with(tmledata, plogis(qlogis(Q1k) + H1 * eps_q))
        tmledata$QAk <- with(tmledata, plogis(qlogis(QAk) + HA * eps_q))
        # tmledata$Qk=qfluc$update
    }
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

#' @export
ey1_estimate2 <- function(tmledata) {
    
    psi <- mean(tmledata$Q1k)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - QAk) + Q1k - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
} 
