sim_msjs=function(N,pent,pi,phi,p,delta=NULL,strata.labels=NULL)
{
  #N is a single value
  #pent is a vector of probabilities that sum to 1 over occasions (nocc)
  #pi is a vector of probabilities that sum to 1; the number of elements defines the number of strata (ns)
  #phi can be a matrix (nS by nocc-1) to have strata-time varying survival or a vector (nS for stratum varying) or single constant value
  #p can be a matrix (nS by nocc) to have strata-time varying capture probabilities or a vector (nS for stratum varying) or single constant value
  #delta - if NULL, no uncertainty is included. If non-NULL it should be a vector whose length is nS and values are probability of observing strata
  # value (getting it right) given it was seen. A "u" is put in place of strata if uncertain.
  #strata.labels only needed to specify non-numeric single character labels
  #
  # number of capture occasions
  nocc=length(pent)
  # number of strata
  nS=length(pi)
  # check values
  if(sum(pent)!=1)stop("sum of pent probabilities must be 1")
  if(sum(pi)!=1)stop("sum of pi probabilities must be 1")
  if(!is.null(delta))
    if(length(delta)!=nS) stop("length of delta must match ns")
  if(is.matrix(phi))
  {
    if(nrow(phi)!=nS)stop("number of rows of phi must equal number of strata")
    if(ncol(phi)!=nocc-1)stop("number of columns of phi must equal number of occasions-1")
  } else
  {
    if(length(phi)==1)
      phi=matrix(phi,nrow=nS,ncol=nocc-1)
    else
      if(length(phi)==nS)
        phi=matrix(phi,nrow=nS,ncol=nocc-1)
      else
        stop("If phi is a vector it must be a single value or equal to the number of strata")
  }
  if(is.matrix(p))
  {
    if(nrow(p)!=nS)stop("number of rows of p must equal number of strata")
    if(ncol(p)!=nocc)stop("number of columns of p must equal number of occasions")
  } else
  {
    if(length(p)==1)
      p=matrix(p,nrow=nS,ncol=nocc)
    else
      if(length(p)==nS)
        p=matrix(p,nrow=nS,ncol=nocc)
      else
        stop("If p is a vector it must be a single value or equal to the number of strata")
  }
  if(!is.null(strata.labels))
  {
    if(length(strata.labels)!=nS)stop("length of strata.labels must match length of pi")
    
  }
  # Generate random entry occasion for N individuals
  ch=rmultinom(1,N,pent)
  nocc=length(pent)
  # Generate stratum for entry
  # ch (nS by nocc matrix) will be the number of entrants in each stratum (row) before each occasion
  ch=apply(ch,1,function(x) rmultinom(1,x,pi))
  # Loop over each cohort and construct random survival and capture process to create capture history matrix (chmat)
  chmat=NULL
  for(i in 1:ncol(ch))
  {
    # Survival process - once dead stays dead which is what cumprod function does
    if(i==1)
    {
      xmat=matrix(rep(1:nS,ch[,i]),nrow=sum(ch[,i]),ncol=nocc)
      smat=cbind(rep(1,sum(ch[,i])),matrix(rbinom(sum(ch[,i])*(nocc-1),1,as.vector(phi[xmat[,i],])),nrow=sum(ch[,i]),ncol=nocc-1))
      #smat=cbind(rep(1,sum(ch[,i])),matrix(rbinom(sum(ch[,i])*(nocc-1),1,phi),nrow=sum(ch[,i]),ncol=nocc-1))
      smat=t(apply(smat,1,cumprod))
      xmat=xmat*smat
    }
    else
    {
      xmat=cbind(matrix(0,nrow=sum(ch[,i]),ncol=i-1),matrix(rep(1:nS,ch[,i]),nrow=sum(ch[,i]),ncol=nocc-i+1))
      if(i<nocc)
      {
        smat=cbind(matrix(1,nrow=sum(ch[,i]),ncol=i),matrix(rbinom(sum(ch[,i])*(nocc-i),1,as.vector(phi[xmat[,i],i:(nocc-1)])),nrow=sum(ch[,i]),ncol=nocc-i))
        smat=t(apply(smat,1,cumprod))
        xmat=xmat*smat
      }
    }
    # Capture process
    pmat=matrix(rbinom(nrow(xmat)*nocc,1,as.vector(p[xmat[,i],])),nrow=nrow(xmat),ncol=nocc)
    pmat=pmat*xmat
    # state uncertainty
    if(!is.null(delta))
    {
      dmat=matrix(rbinom(nrow(xmat)*nocc,1,1-as.vector(delta[xmat[,i]])),nrow=nrow(xmat),ncol=nocc)
      pmat[pmat>0&dmat==0]=nS+1
    }
    chmat=rbind(chmat,pmat)
    
  }
  # add labels if any
  if(!is.null(delta))strata.labels=c(strata.labels,"u")
  if(!is.null(strata.labels))
    chmat[chmat!=0]=strata.labels[chmat[chmat!=0]]
  ch=apply(chmat,1,paste,collapse="")
  ch=data.frame(table(ch))
  ch$ch=as.character(ch$ch)
  colnames(ch)=c("ch","freq")
  return(ch[-1,])
}

