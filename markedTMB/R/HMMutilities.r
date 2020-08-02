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
		for(t in (first[i]+1):T)
		{
			psi[t,]=apply(psi[t-1,]*gamma[i,t-1,,],2,max)%*%diag(dmat[i,t,x[i,t],])
		}
		states[i,T]=which.max(psi[T,])
		for(t in (T-1):first[i])
			states[i,t]=which.max(psi[t,]*gamma[i,t,,states[i,t+1]])
	}
	states=t(apply(states,1,function(x) state.names[x]))
	return(states)
}
