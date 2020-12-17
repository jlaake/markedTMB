#' Fitting function for Multistate CJS models with TMB
#' 
#' A function for computing MLEs for a Multi-state Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{mscjs_tmb} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of simplified dataframes for design data
#' @param fullddl list of complete dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, S.dm,p.dm,Psi.dm,S.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix all parameters to speed up execution
#' @param initial list of initial values for parameters if desired; if each is a named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optim}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param re if TRUE creates random effect model admbcjsre.tpl and runs admb optimizer
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to tmb
#' @param clean if TRUE, deletes the dll and recompiles 
#' @param getreals if TRUE, compute real values and std errors for TMB models; may want to set as FALSE until model selection is complete
#' @param useHess if TRUE, the TMB hessian function is used for optimization; using hessian is typically slower with many parameters but can result in a better solution
#' @param savef if TRUE, save optimization function in model for reporting
#' @param ... not currently used
#' @export
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optim}; if 0 \code{optim} thinks it
#' converged} \item{count}{\code{optim} results of number of function
#' evaluations} \item{reals}{dataframe of data and real S and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake
#' @references Ford, J. H., M. V. Bravington, and J. Robbins. 2012. Incorporating individual variability into mark-recapture models. Methods in Ecology and Evolution 3:1047-1054.
#' @examples 
#' \donttest{
#' # this example requires TMB
#' # The same example is in the RMark package and it is included here to
#' # illustrate the differences in the handling of mlogit parameters between RMark 
#' # and marked.  The MARK software handles parameters like Psi which must sum to 1
#' # by excluding one of the cells that is used as a reference cell and is computed by
#' # subtracting the other cell values from 1 so the total sums to 1.  This is often
#' # handled with an mlogit parameter in which the cell values are exp(beta) and the
#' # reference cell is set to 1 and the values are divided by the sum across the cells
#' # so the resulting values are probabilities that sum to 1. In marked, instead of removing
#' # one of the cells, all are included and the user must select which should be the
#' # reference cell by setting the value fix=1 for that cell and others are NA so they are
#' # estimated. For transition parameters like Psi, the default design data is setup so 
#' # that the probability of remaining in the cell (stratum=tostratum) is the reference cell
#' # and fix set to 1.  Thus, this means 2 changes are needed to the script in RMark.
#' # The first is to remove the statement skagit.ddl$Psi$fix=NA because that over-rides
#' # the default fix values.  The other is to add
#' # skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
#' #  skagit.ddl$Psi$time==5]=0
#' # to change the value from 1 to 0 which forces movement from B to A in the interval 5 to 6. If
#' # this is not done then Psi B to B=Psi B to A=0.5 because each is 1 and when they are normalized
#' # they are divided by the sum which is 2 (1/2).
#' data(skagit)
#' skagit.processed=process.data(skagit,model="MSCJS",groups=c("tag"),strata.labels=c("A","B"))
#'skagit.ddl=make.design.data(skagit.processed)
#'#
#'# p
#'#
#'# Can't be seen at 5A or 2B,6B (the latter 2 don't exist)
#'skagit.ddl$p$fix=ifelse((skagit.ddl$p$stratum=="A"&skagit.ddl$p$time==5) | 
#' (skagit.ddl$p$stratum=="B"&skagit.ddl$p$time%in%c(2,6)),0,NA)
#'# Estimated externally from current data to allow estimation of survival at last interval
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v7"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.687
#'skagit.ddl$p$fix[skagit.ddl$p$tag=="v9"&skagit.ddl$p$time==6&skagit.ddl$p$stratum=="A"]=0.975
#'#
#'# Psi
#'#
#'# only 3 possible transitions are A to B at time interval 2 to 3 and 
#'# for time interval 3 to 4 from A to B and from B to A
#'# rest are fixed values
#'############ change for RMark to marked; remove next line
#'#skagit.ddl$Psi$fix=NA
#'# stay in A for intervals 1-2, 4-5 and 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="A"&
#'  skagit.ddl$Psi$tostratum=="B"&skagit.ddl$Psi$time%in%c(1,4,5)]=0
#'# stay in B for interval 4-5
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"
#'  &skagit.ddl$Psi$time==4]=0
#'# leave B to go to A for interval 5-6
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#' skagit.ddl$Psi$time==5]=1
#'############ change for RMark to marked; add next line to set B to B to 0 otherwise it has
#'############ been set to 1 by default which would make psi B to B = psi B to A = 0.5
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="B"&
#' skagit.ddl$Psi$time==5]=0
#'# "stay" in B for interval 1-2 and 2-3 because none will be in B
#'skagit.ddl$Psi$fix[skagit.ddl$Psi$stratum=="B"&skagit.ddl$Psi$tostratum=="A"&
#' skagit.ddl$Psi$time%in%1:2]=0
#'# 
#'# S
#'#
#'# None in B, so fixing S to 1
#'skagit.ddl$S$fix=ifelse(skagit.ddl$S$stratum=="B"&skagit.ddl$S$time%in%c(1,2),1,NA)
#'skagit.ddl$S$fix[skagit.ddl$S$stratum=="A"&skagit.ddl$S$time==4]=1
#'# fit model
#'p.timexstratum.tag=list(formula=~time:stratum+tag,remove.intercept=TRUE)
#'Psi.sxtime=list(formula=~-1+stratum:time)
#'S.stratumxtime=list(formula=~-1+stratum:time)
#'#
#'mod1=crm(skagit.processed,skagit.ddl,
#' model.parameters=list(S=S.stratumxtime,p= p.timexstratum.tag,Psi=Psi.sxtime),hessian=TRUE)
#'}
mscjs_tmb=function(x,ddl,fullddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		re=FALSE,compile=FALSE,extra.args="",clean=TRUE,getreals=FALSE, useHess=FALSE,savef=TRUE, ...)
{
	accumulate=FALSE
	nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use x$intervals if none are in ddl$Phi
	if(!is.null(ddl$S$time.interval))		
	  time.intervals=matrix(fullddl$S$time.interval[fullddl$S$stratum==x$strata.labels[1]],nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
	if(is.vector(x$time.intervals))
		time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
		time.intervals=x$time.intervals

#  Store data from x$data into x
	strata.labels=x$strata.labels
	uS=x$unobserved
#  set default frequencies if not used
	freq=NULL
	if(!is.null(x$data$freq))freq=x$data$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
	ch=x$data$ch
	imat=process.ch(ch,freq,all=FALSE)
	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=length(ch))
	for(nlabel in 1:length(strata.labels))
		chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
	chmat=t(apply(chmat,1,function(x) as.numeric(x)))
#  Use specified initial values or create if null
	if(is.null(initial))
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
				p=rep(0,ncol(dml$p$fe)),
				S=rep(0,ncol(dml$S$fe)))
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(S.dm=dml$S$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,imat=imat,S.fixed=parameters$S$fixed,
			p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
	if(accumulate)
	{
		cat("Accumulating capture histories based on design. This can take awhile.\n")
		flush.console()
		model_data.save=model_data   
		#model_data=mscjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#  Scale the design matrices and parameters with either input scale or computed scale
	scale=1
	scale=set.scale(names(dml),model_data,scale)
	model_data=scale.dm(model_data,scale)
	setup_tmb("multistate_tmb",clean=clean)

	# S design matrix
	phidm=as.matrix(model_data$S.dm)
	phifix=rep(-1,nrow(phidm))
	if(!is.null(ddl$S$fix))
		phifix[!is.na(ddl$S$fix)]=ddl$S$fix[!is.na(ddl$S$fix)]
	phi_slist=simplify_indices(cbind(phidm,phifix))
	phi_relist=setup_re(fullddl$S,parameters$S$formula)
	
	# p design matrix
	pdm=as.matrix(model_data$p.dm)
	pfix=rep(-1,nrow(pdm))
	if(!is.null(ddl$p$fix))
		pfix[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
	p_slist=simplify_indices(cbind(pdm,pfix))
	p_relist=setup_re(fullddl$p,parameters$p$formula)
	
	#Psi design matrix
	psidm=as.matrix(model_data$Psi.dm)
	psifix=rep(-1,nrow(psidm))
	if(!is.null(ddl$Psi$fix))
		psifix[!is.na(ddl$Psi$fix)]=ddl$Psi$fix[!is.na(ddl$Psi$fix)]
	psi_slist=simplify_indices(cbind(psidm,psifix))
	psi_relist=setup_re(fullddl$Psi,parameters$Psi$formula)
	

	f = MakeADFun(data=list(n=length(model_data$imat$freq),m=model_data$imat$nocc,nS=length(strata.labels),
					ch=chmat,frst=model_data$imat$first,freq=model_data$imat$freq,tint=model_data$time.intervals,
					nrowphi=length(phi_slist$set),	phidm=phidm[phi_slist$set,,drop=FALSE],
					phifix=phifix[phi_slist$set],phiindex=phi_slist$indices[ddl$S.indices],
					phi_nre=phi_relist$nre,phi_krand=phi_relist$krand,phi_randDM=phi_relist$randDM,phi_randDM_i=phi_relist$randDM_i,
					phi_randIndex=phi_relist$randIndex,phi_randIndex_i=phi_relist$randIndex_i,phi_counts=phi_relist$counts,phi_idIndex=phi_relist$idIndex,
					phi_idIndex_i=phi_relist$idIndex_i,
					nrowp=length(p_slist$set),pdm=pdm[p_slist$set,,drop=FALSE],pfix=pfix[p_slist$set],pindex=p_slist$indices[ddl$p.indices],
					p_nre=p_relist$nre,p_krand=p_relist$krand,p_randDM=p_relist$randDM, p_randDM_i=p_relist$randDM_i, 
					p_randIndex=p_relist$randIndex,p_randIndex_i=p_relist$randIndex_i,p_counts=p_relist$counts,
					p_idIndex=p_relist$idIndex,p_idIndex_i=p_relist$idIndex_i,
					nrowpsi=length(psi_slist$set),	psidm=psidm[psi_slist$set,,drop=FALSE],
					psifix=psifix[psi_slist$set],psiindex=psi_slist$indices[ddl$Psi.indices],
					psi_nre=psi_relist$nre,psi_krand=psi_relist$krand,psi_randDM=psi_relist$randDM,psi_randDM_i=psi_relist$randDM_i,
					psi_randIndex=psi_relist$randIndex,psi_randIndex_i=psi_relist$randIndex_i,psi_counts=psi_relist$counts,psi_idIndex=psi_relist$idIndex,
					psi_idIndex_i=psi_relist$idIndex_i,
					getreals=as.integer(getreals)),
			        parameters=list(phibeta=par$S,pbeta=par$p,psibeta=par$Psi,log_sigma_phi=rep(-1,phi_relist$nsigma),
					log_sigma_p=rep(-1,p_relist$nsigma),log_sigma_psi=rep(-1,psi_relist$nsigma),u_phi=rep(0,phi_relist$nre),
					u_p=rep(0,p_relist$nre),u_psi=rep(0,psi_relist$nre)),random=c("u_phi","u_p","u_psi"),DLL="multistate_tmb")
	cat("\nrunning TMB program\n")                         
	if(method=="nlminb")
	{
	  if(!useHess)
		   mod=nlminb(f$par,f$fn,f$gr,control=control,...)
	  else
	    mod=nlminb(f$par,f$fn,f$gr,f$he,control=control,...)
		lnl=mod$objective
		par=mod$par
		convergence=mod$convergence
	} else
	{
	  if(method=="SANN")
	  {		  
	    control$maxit=itnmax
	    mod=optim(f$par,f$fn,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
	    par=mod$par
	    convergence=mod$convergence
	  } else
	  {
	    control$starttests=FALSE
 		  if(!useHess)
		     mod=optimx(f$par,f$fn,f$gr,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
		  else
		     mod=optimx(f$par,f$fn,f$gr,f$he,hessian=FALSE,control=control,itnmax=itnmax,method=method,...)
		  par <- coef(mod, order="value")[1, ]
		  mod=as.list(summary(mod, order="value")[1, ])
	    convergence=mod$convcode
	  }
	  lnl=mod$value
	}
	fixed.npar=ncol(phidm)+ncol(pdm)+ncol(psidm)
	par_summary=sdreport(f,getJointPrecision=getreals&hessian)
	S_re=NULL
	p_re=NULL
	Psi_re=NULL
	if(p_relist$nre+phi_relist$nre+psi_relist$nre>0)
	{
	  random_values=par_summary$par.random
	 	par=par_summary$par.fixed[1:fixed.npar]
		cjs.beta.fixed=unscale.par(par,scale)
		cjs.beta.sigma=par_summary$par.fixed[-(1:fixed.npar)]
		sigma=NULL
		if(phi_relist$krand>0)
		{
			Phi_sigma=cjs.beta.sigma[1:phi_relist$krand]
			names(Phi_sigma)=colnames(phi_relist$randDM)
			sigma=list(S_logsigma=Phi_sigma)
			S_re=split(random_values[names(random_values)%in%"u_phi"],rep(1:length(phi_relist$nre_byeffect),phi_relist$nre_byeffect))
		} 
		if(p_relist$krand>0)
		{
			p_sigma=cjs.beta.sigma[(phi_relist$krand+1):(phi_relist$krand+p_relist$krand)]
			names(p_sigma)=colnames(p_relist$randDM)
			sigma=c(sigma,list(p_logsigma=p_sigma))
			p_re=split(random_values[names(random_values)%in%"u_p"],rep(1:length(p_relist$nre_byeffect),p_relist$nre_byeffect))
		}
		if(psi_relist$krand>0)
		{
			Psi_sigma=cjs.beta.sigma[(phi_relist$krand+p_relist$krand+1):(phi_relist$krand+p_relist$krand+psi_relist$krand)]
			names(Psi_sigma)=colnames(psi_relist$randDM)
			sigma=c(sigma,list(Psi_logsigma=Psi_sigma))
			Psi_re=split(random_values[names(random_values)%in%"u_psi"],rep(1:length(psi_relist$nre_byeffect),psi_relist$nre_byeffect))
		} 
		cjs.beta=c(cjs.beta.fixed,sigma)
		beta.vcv=par_summary$cov.fixed
		colnames(beta.vcv)=names(unlist(cjs.beta))
		rownames(beta.vcv)=colnames(beta.vcv)
	}else
	{	
		cjs.beta=unscale.par(par,scale)
		if(hessian) 
		{
			message("Computing hessian...")
		  beta.vcv=solvecov(f$he(par))$inv
			colnames(beta.vcv)=names(unlist(cjs.beta))
			rownames(beta.vcv)=colnames(beta.vcv)
		} else
			beta.vcv=NULL
	}	
	# Create results list 
	if(getreals)
	{
		reals=split(par_summary$value,names(par_summary$value))
		reals.se=split(par_summary$sd,names(par_summary$value))	
		names(reals)=c("p","S","Psi")
		names(reals.se)=c("p","S","Psi")
		reals$S[fullddl$S$Time<fullddl$S$Cohort]=NA
		reals.se$S[fullddl$S$Time<fullddl$S$Cohort]=NA
		reals$p[fullddl$p$Time<fullddl$p$Cohort]=NA
		reals.se$p[fullddl$p$Time<fullddl$p$Cohort]=NA
		reals$Psi=as.vector(aperm(array(reals$Psi,dim=c(model_data$imat$nocc-1,length(strata.labels),length(strata.labels),length(model_data$imat$freq))),c(3,2,1,4)))
		reals.se$Psi=as.vector(aperm(array(reals.se$Psi,dim=c(model_data$imat$nocc-1,length(strata.labels),length(strata.labels),length(model_data$imat$freq))),c(3,2,1,4)))
		reals$Psi[fullddl$Psi$Time<fullddl$Psi$Cohort]=NA
		reals.se$Psi[fullddl$Psi$Time<fullddl$Psi$Cohort]=NA
	}
	else
	{
		reals=NULL
		reals.se=NULL
	}
	res=list(beta=cjs.beta,re=list(S=S_re,p=p_re,Psi=Psi_re),neg2lnl=2*lnl,AIC=2*lnl+2*sum(sapply(cjs.beta,length)),
			beta.vcv=beta.vcv,TMBreals=reals,TMBreals.se=reals.se,convergence=convergence,optim.details=mod,
			model_data=model_data,options=list(scale=scale,accumulate=accumulate,initial=initial,method=method,
			chunk_size=chunk_size,itnmax=itnmax,control=control))		
		
#  Restore non-accumulated, non-scaled dm's etc
	res$model_data=model_data.save
	# if savef add it to the model
	if(savef)res$f=f
#  Assign S3 class values and return
	class(res)=c("crm","mscjs")
	return(res)
}




