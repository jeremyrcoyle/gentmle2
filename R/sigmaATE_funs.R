sigmaATE_update=function(tmledata,Q.trunc = 0.001, ...){
  #fix points where Q is already 0 or 1 - perfect prediction
  subset=with(tmledata,which(0<Qk & Qk<1))  
  eps_q=0
  tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
  if(length(subset)>0){
    #fluctuate Q
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qktrunc)))
    eps_q=qfluc$eps
    tmledata$Qk=with(tmledata,plogis(qlogis(Qktrunc)+HA*eps_q))
    tmledata$Q1k=with(tmledata,plogis(qlogis(Q1k)+H1*eps_q))
    tmledata$Q0k=with(tmledata,plogis(qlogis(Q0k)+H0*eps_q))
  } 
  list(tmledata=tmledata,coef=eps_q)  
  
}

sigmaATE_estimate=function(tmledata,Q.trunc = 0.001, ...){
  
  ATE=mean(with(tmledata,Q1k-Q0k))
  tmledata$HA=with(tmledata,2*(Q1k-Q0k-ATE)*(A/gk-(1-A)/(1-gk)))
  tmledata$H1=with(tmledata,2*(Q1k-Q0k-ATE)/gk)
  tmledata$H0=with(tmledata,-2*(Q1k-Q0k-ATE)/(1-gk))
  
  #influence curve
  sigma <- var(with(tmledata,Q1k-Q0k))
  Dstar_sigma <- with(tmledata,(Q1k-Q0k-ATE)^2-sigma+HA*(Y-Qk))
  list(tmledata=tmledata,ests=c(sigma=sigma),Dstar=list(Dstar_sigma=Dstar_sigma))
}

