library(markedTMB)
reps=100
N=500
reals=vector("list",length=reps)
for(i in 1:reps)
{
  xx=sim_msjs(N,pent=c(.1,.2,.1,.2,.1,.1,.1,.1),pi=c(0.1,0.9),
              phi=0.7,p=c(0.6,0.4),strata.labels=c("A","B"))
  
  dp=process.data(xx,model="MSJS",strata.labels=c("A","B"))
  ddl=make.design.data(dp)
  mod_js=crm(dp,ddl,model.parameters=list(S=list(formula=~1),
                                          p=list(formula=~-1+stratum),pi=list(formula=~-1+stratum),Psi=list(formula=~1),
                                          pent=list(formula=~-1+time)))
  reals[[i]]=sapply(mod_js$results$beta,plogis)
  reals[[i]]$pi=exp(c(mod_js$results$beta$pi))/(sum(exp(mod_js$results$beta$pi)))
  reals[[i]]$pent=exp(c(0,mod_js$results$beta$pent))/(1+sum(exp(mod_js$results$beta$pent)))
  
}

mean(sapply(reals,function(x)x$S))
apply(sapply(reals,function(x)x$p),1,mean)
apply(sapply(reals,function(x)x$pent),1,mean)
apply(sapply(reals,function(x)x$pi),1,mean)

