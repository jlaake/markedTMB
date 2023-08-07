#' Compute HMM matrices
#' 
#' Computes gamma, dmat and delta (initial) matrices(arrays) and returns them in a list.
#' 
#' @param object fitted crm model
#' @param state.names if not NULL use these labels instead of strata.labels in model
#' @param par if not NULL use these values for matrix calculations
#' @return list with gamma, dmat and delta arrays
#' @author Jeff Laake
#' @export compute_matrices
#' @keywords utility
compute_matrices=function(object,state.names=NULL,par=NULL)
{
  load_lib(object$model)
  # use last.par or specified parameters
  if(is.null(par))par=object$results$f$env$last.par
  mat=object$results$f$report(par)
  # depending on model adjust matrices
	if(object$model=="MSLD")
	{
	  names(mat)=c("dmat","gamma")
	  mat$delta=matrix(0,nrow=nrow(object$data$data),ncol=length(object$data$strata.labels)*2+1)
	  for(i in 1:nrow(object$data$data))
	    mat$delta[i,object$data$start[i,1]]=1
	  # Adjust dmat to have nocc+1 entries for time to work with global decode
	  temp=array(0,dim=c(dim(mat$dmat)[1],dim(mat$dmat)[2]+1,dim(mat$dmat)[3],dim(mat$dmat)[4]))
	  temp[,2:dim(temp)[2],,]=mat$dmat
	  for(i in 1:nrow(temp))
	    for(j in 1:length(object$data$strata.labels))
	    {
	      temp[i,object$data$start[i,2],j+1,j]=1
	      temp[i,object$data$start[i,2],1,j]=0
	    }
	  mat$dmat=temp
	  mat=list(delta=mat$delta,dmat=mat$dmat,gamma=mat$gamma)
	}
  if(object$model=="MSCJS")
  {
    names(mat)=c("dmat","gamma")
    mat$delta=matrix(0,nrow=nrow(object$data$data),ncol=length(object$data$strata.labels)+1)
    for(i in 1:nrow(object$data$data))
      mat$delta[i,object$data$start[i,1]]=1
    temp=array(0,dim=c(dim(mat$dmat)[1],dim(mat$dmat)[2]+1,dim(mat$dmat)[3],dim(mat$dmat)[4]))
    temp[,2:dim(temp)[2],,]=mat$dmat
    for(i in 1:nrow(temp))
      for(j in 1:length(object$data$strata.labels))
      {
        temp[i,object$data$start[i,2],j+1,j]=1
        temp[i,object$data$start[i,2],1,j]=0
      }
    mat$dmat=temp
    mat=list(delta=mat$delta,dmat=mat$dmat,gamma=mat$gamma)
  }
  if(object$model%in%c("MSJS","MSJSU"))
  {
    names(mat)=c("N","p0","dmat","gamma")
    mat$delta=matrix(0,nrow=nrow(object$data$data),ncol=length(object$data$strata.labels)+1)
    mat$delta[,1]=1
    object$data$start=rbind(c(1,1),object$data$start)
    temp=array(0,dim=c(dim(mat$dmat)[1],dim(mat$dmat)[2]+1,dim(mat$dmat)[3],dim(mat$dmat)[4]))
    temp[,2:dim(temp)[2],,]=mat$dmat
    for(i in 1:nrow(temp))
      for(j in 1:length(object$data$strata.labels))
      {
        temp[i,object$data$start[i,2],j+1,j]=1
        temp[i,object$data$start[i,2],1,j]=0
      }
    mat$dmat=temp
    #mat=list(delta=mat$delta,dmat=mat$dmat,gamma=mat$gamma)
  }
  if(object$model=="MVMSCJS")
  {    
    names(mat)=c("delta","dmat","gamma")
    mat$delta=cbind(mat$delta,rep(0,nrow(mat$delta)))
  }
   if(is.null(state.names))
  {
    if(object$model=="MSLD")
      state.names=c(object$data$strata.labels,paste(object$data$strata.labels,":NewDead",sep=""),"Dead")
    else
      state.names=c(object$data$strata.labels,"Dead")
  }
	if(length(state.names)!=dim(mat$gamma)[3]) stop("\n length of state.names not correct. Should be ",dim(mat$gamma)[3])
	obs.names=object$data$ObsLevels
	dimnames(mat$gamma)[3:4]=c(list(state.names),list(state.names))
	names(dimnames(mat$gamma))=c("Id","Occasion","From_state","To_state")
	dimnames(mat$dmat)[3:4]=c(list(obs.names),list(state.names))
	names(dimnames(mat$dmat))=c("Id","Occasion","Observation","State")
	colnames(mat$delta)=state.names
	return(mat)
}
#' Computes backward probabilities 
#' 
#' Computes backward probability sequence for a set of capture histories
#' 
#' @param object fitted crm model 
#' @param state.names if not NULL use these labels instead of strata.labels in model
#' @param par if not NULL use these values for matrix calculations
#' @author Jeff Laake
#' @return array of backward probabilities (one for each id, state, occasion)
#' @export backward_prob
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 61.
backward_prob=function(object,state.names=NULL,par=NULL)
{  	
	if(!is.null(object$results$mat))
	{
		dmat=object$results$mat$dmat
		gamma=object$results$mat$gamma
	}else
	{
	  load_lib(object$model)
	  matlist=compute_matrices(object=object,state.names=state.names,par=par)
		dmat=matlist$dmat
		gamma=matlist$gamma
	}
	x=object$data$ehmat
	T=object$data$nocc
	first=object$data$start[,2]
	m=object$data$m
	beta=array(NA,dim=c(nrow(x),ncol(x),m))
	# Loop over capture histories
	for(i in 1:nrow(x))
	{
		occ=T
		beta[i,occ,]=rep(1,m)
		# Loop over occasions for this encounter history (x)
		for(t in T:(first[i]+1))
		{
			occ=occ-1
			# Compute backward probability for this occasion
			beta[i,occ,]=gamma[i,t-1,,]%*%diag(dmat[i,t,x[i,t],])%*%beta[i,occ+1,]  
		}
	}
	return(beta)
} 
#' Local decoding of HMM 
#' 
#' Computes state predictions one at a time for each occasion for each individual
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param state.names names for states if over-riding defaults
#' @param par if not NULL use these values for matrix calculations
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export local_decode
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 80.

local_decode=function(object,state.names=NULL,par=NULL)
{  	
  load_lib(object$model)
  if(is.null(par)) 
    result=list(lnl=object$results$neg2lnl)
  else
    result=list(lnl=object$results$f$fn(par))
	result$beta=backward_prob(object,state.names=state.names,par=par)
	stateprob=result$alpha*result$beta/exp(result$lnl)
	states=apply(stateprob,c(1,2),function(x){ if(any(is.na(x))) return(NA) else return(state.names[which.max(x)])})
	return(states)
} 
#' Global decoding of HMM 
#' 
#' Computes sequence of state predictions for each individual
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param state.names names to over-ride default names for states
#' @param par if not NULL use these values for matrix calculations
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export global_decode
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 82.
global_decode=function(object,state.names=NULL,par=NULL)
{  	
  load_lib(object$model)
  if(is.null(object$results$mat))
      parmlist=compute_matrices(object,state.names=state.names,par=par)
  else
      parmlist=object$results$mat
  dmat=parmlist$dmat
  gamma=parmlist$gamma
  state.names=dimnames(gamma)[3][[1]]
  delta=parmlist$delta
	x=object$data$ehmat
	T=object$data$nocc
	if(is.null(object$data$m))
	  m=dim(dmat)[4]
	else
	  m=object$data$m
	first=object$data$start[,2]
	states=matrix(NA,nrow=nrow(x),ncol=T)
	for(i in 1:nrow(x))
	{
		psi=matrix(NA,nrow=T,ncol=m)
		# Assign psi value at first occasion
		psi[first[i],]=delta[i,]%*%diag(dmat[i,first[i],x[i,first[i]],]) 
		psi[first[i],]=psi[first[i],]/sum(psi[first[i],])
		for(t in (first[i]+1):T)
		{
			psi[t,]=apply(psi[t-1,]*gamma[i,t-1,,],2,max)%*%diag(dmat[i,t,x[i,t],])
			psi[t,]=psi[t,]/sum(psi[t,])
		}
		states[i,T]=which.max(psi[T,])
		for(t in (T-1):first[i])
			states[i,t]=which.max(psi[t,]*gamma[i,t,,states[i,t+1]])
		# check to make sure predicted state matches any observed state
		observed=object$data$ehmat[i,]-1
		if(any(states[i,]!=observed & observed!=0))
		{
		  warning("\nFor the following observation the predicted states conflict with the observed states")
		  cat("\n",object$data$data[i,])
		  cat("\npredicted")
		  cat("\n",states[i,])
		}
	}
	states=t(apply(states,1,function(x) state.names[x]))
	return(states)
}

#' Estimate of abundance
#' 
#' Computes total abundance by occasion and abundance by stratum and occasion
#' for fitted MSJS and MSJSU models
#' 
#' @param object fitted crm model (must be an MSJS or MSJSU model)
#' @param par if not NULL use these values for matrix calculations
#' @param state.names names to over-ride default names for states
#' @author Jeff Laake
#' @return dataframe of abundance estimates
#' @export abundance_estimate
#' @keywords utility
abundance_estimate=function(object,par=NULL,state.names=NULL)
{
  # if par is NULL use the fitted parameter values otherwise it uses those
  # specified in par which is used for bootstrapping
  if(!object$model%in%c("MSJS","MSJSU")) stop("Only used for JS type models")
  # compute dmat and gamma matrices - then get predicted states
  parmlist=compute_matrices(object,state.names=state.names,par=par)
  object$results$mat=parmlist
  gd=global_decode(object,par=par)
  # get frequencies of each history ignoring first one which is all 0's
  freq=object$data$freq[-1]
  # replicate each set of predicted states by frequency
  gd=gd[rep(1:nrow(gd),freq),]
  # get count of predicted states by strata over columns(occasions)
  strata=object$data$strata.labels
  strata=strata[strata!="N"]
  xx=apply(gd,2,function(x)table(factor(x,levels=strata)))
  # compute abundance by multiplying total abundance by proportion in each state for each occasion
  # this excludes any that are dead or still in state N
  est=rep(object$results$mat$N,object$data$nocc)*xx/sum(freq)
  # remove first occasion which is prior to sampling occasions (all in N state)
  est=est[,-1]
  # get totals over strata and create dataframe
  est=rbind(est,colSums(est))
  est=floor(est)
  rownames(est)[3]="Total"
  N=data.frame(occasion=rep(1:(object$data$nocc-1),nrow(est)),strata=rep(rownames(est),each=10),N=as.vector(t(est)))
  return(N)
}

#' Bootstrap limits for abundance estimation 
#' 
#' Computes bootstrap std error and confidence limits for estimated abundance in a MSJS or MSJSU model
#' 
#' @param object fitted crm model of type MSJS or MSJSU
#' @param nreps number of bootstrap replicates
#' @param alpha value for confidence interval 
#' @author Jeff Laake
#' @return dataframe of abundance estimates and std error and confidence intervals
#' @export abundance_estimate_bs
#' @keywords utility
abundance_estimate_bs=function(object,nreps=1000,alpha=0.05)
{
  if(!object$model%in%c("MSJS","MSJSU")) stop("Only used for JS type models")
  if(alpha>=1 | alpha <=0)stop("Invalid value of alpha")
  # get abundance estimate dataframe with fitted parameter values
  N=abundance_estimate(object)
  #generate nreps vectors of parameters from multivariate-normal
  par_reps=rmvnorm(nreps,mean=unlist(object$results$beta),sigma=object$results$beta.vcv)
  # generate abundance estimate for first replicate and then create matrix with
  # a column for each estimate
  est=abundance_estimate(object,par=par_reps[1,])
  boot_reps=matrix(NA,ncol=nrow(est),nrow=nreps)
  colnames(boot_reps)=paste(est$occasion,est$strata,sep="-")
  boot_reps[1,]=est$N
  # loop over remaining replicates and store in matrix boot_reps
  for(i in 2:nreps)
    boot_reps[i,]=abundance_estimate(object,par=par_reps[i,])$N
  # construct integer values for lower and upper limit of sorted replicate values
  Llimit=floor(nreps*alpha/2)
  Ulimit=floor(nreps*(1-alpha/2))
  # get confidence limits for each estimate
  CLimits=apply(boot_reps,2,function(x) {
    xsort=sort(x)
    c(xsort[Llimit],xsort[Ulimit])
  })
  # get std error for each estimate
  se=apply(boot_reps,2,function(x) sqrt(var(x)))
  # return estimates with se and limits
  N=cbind(N,se,t(CLimits))
  colnames(N)[4:6]=c("SE","LCL","UCL")
  return(N)
}

