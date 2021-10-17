#' Capture-recapture model fitting function
#'
#' Fits user specified models to some types of capture-recapture wholly in R
#' and not with MARK.  A single function that processes data, creates the
#' design data, makes the crm model and runs it
#'
#' This function is operationally similar to the function \code{mark in RMark}
#' in that is is a shell that calls several other functions to perform the
#' following steps: 1) \code{\link{process.data}} to setup data and parameters
#' and package them into a list (processed data),2)
#' \code{\link{make.design.data}} to create the design data for each parameter
#' in the specified model, 3) \code{\link{create.dm}} to create the design
#' matrices for each parameter based on the formula provided for each
#' parameter, 4) call to the specific function for model fitting. As with \code{mark} the calling
#' arguments for \code{crm} are a compilation of the calling arguments for each
#' of the functions it calls (with some arguments renamed to avoid conflicts).
#' expects to find a value for \code{ddl}.  Likewise, if the data have not been
#' processed, then \code{ddl} should be NULL.  This dual calling structure
#' allows either a single call approach for each model or alternatively the preferred method
#' where the data area processed and the design data (\code{ddl}) created once and
#' then a whole series of models can be analyzed without repeating those steps.
#'
#' There are some optional arguments that can be used to set initial values and
#' control other aspects of the optimization.  The optimization is done with
#' the R package/function \code{optimx} and the arguments \code{method} and
#' \code{hessian} are described with the help for that function.  
#'
#' In the current implementation, a logit link is used to constrain the
#' parameters in the unit interval (0,1) except for probability of entry which
#' uses an mlogit and N which uses a log link. For the probitCJS model, a probit link is
#' used for the parameters. These could be generalized to
#' use other link functions. Following the notation of MARK, the parameters in
#' the link space are referred to as \code{beta} and those in the actual
#' parameter space of \code{Phi} and \code{p} as reals.
#'
#' Initial values can be set in 2 ways.  1) Define a list of named vectors with 
#' the initial beta parameter values (eg logit link) in \code{initial}. 
#' The names of the vectors should be the parameter names in the model. Any unspecified
#' values are set to 0. 2) Specify a previously run model for initial. The code will match the names of the
#' current design matrix to the names in \code{beta} and use the appropriate
#' initial values. Any non-specified values are set to 0.  If no value is specified for initial,
#' all beta are started at a value of 0.
#'
#' If you have a study with unequal time intervals between capture occasions,
#' then these can be specified with the argument \code{time.intervals}.
#'
#' The argument \code{accumulate} defaults to \code{TRUE}.  When it is
#' \code{TRUE} it will accumulate common capture histories that also have
#' common design and common fixed values (see below) for the parameters.  This
#' will speed up the analysis because in the calculation of the likelihood
#' it loops over the unique values. In general the
#' default will be the best unless you have many capture histories and are
#' using many individual covariate(s) in the formula that would make each entry
#' unique. In that case there will be no effect of accumulation but the code
#' will still try to accumulate. In that particular case by setting
#' \code{accumulate=FALSE} you can skip the code run for accumulation.
#'
#' Most of the arguments controlling the fitted model are contained in lists in
#' the arguments \code{model.parameters} and \code{design.parameters} which are
#' similar to their counterparts in \code{mark inb RMark}. Each is a named list
#' with the names being the parameters in the model (e.g., Phi and p in "cjs"
#' and "Phi","p","pent","N" in "js"). Each named element is also a list
#' containing various values defining the design data and model for the
#' parameter. The elements of \code{model.parameters} can include
#' \code{formula} which is an R formula to create the design matrix for the
#' parameter and \code{fixed} is a matrix of fixed values as described below.
#' The elements of \code{design.parameters} can include \code{time.varying},
#' \code{fields}, \code{time.bins},\code{age.bins}, and \code{cohort.bins}. See
#' \code{\link{create.dmdf}} for a description of the first 2 and
#' \code{\link{create.dm}} for a description of the last 3.
#'
#' Real parameters can be set to fixed values by specifying values for a field called fix in the design data for a parameter.
#' If the value of fix is NA the parameter is estimated and if it is not NA then the real
#' parameter is fixed at that value.  If you also specify fixed as decribed above, they will over-ride any
#' values you have also set with fix in the design data. To set all of the real values for a
#' particular occasion you can use the following example with the dipper data
#' as a template:
#'
#' \code{model.parameters=list(Phi=list(formula=~1,}
#' \code{fixed=cbind(1:nrow(dipper),rep(2,nrow(dipper)),rep(1,nrow(dipper)))))}
#'
#' The above sets \code{Phi} to 1 for the interval between occasions 2 and 3
#' for all animals.
#'
#' Alternatively, you could do as follows:
#'
#' data(dipper)
#' dp=process.data(dipper)
#' ddl=make.design.data(dp)
#' ddl$Phi$fix=ifelse(ddl$Phi$time==2,1,NA)
#'
#' At present there is no modification of the parameter count
#' to address fixing of real parameters except that if by fixing reals, a beta is not needed in the design it will be dropped.
#' For example, if you were to use ~time for Phi with survival fixed to 1 for time 2, then then beta for that time would not
#' be included.
#'
#'
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model (eg, "cjs", "js")
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional vector of initial values for beta parameters; if
#' named from previous analysis only relevant values are used
#' @param groups Vector of names factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param method optimization method
#' @param debug if TRUE, shows optimization output
#' @param hessian if TRUE, computes v-c matrix using hessian
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce
#' computation
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories and design matrices; amount used is 8*chunk_size/1e6 MB (default
#' 80MB)
#' @param control control string for optimization functions
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations for optimization
#' @param scale vector of scale values for parameters
#' @param run if TRUE, it runs model; otherwise if FALSE can be used to test model build components
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param save.matrices for HMM models this option controls whether the gamma,dmat and delta matrices are saved in the model object
#' @param savef if TRUE, save the makeAdFun result from TMB to report real values and matrices
#' @param getreals if TRUE, compute real values and std errors for TMB models; may want to set as FALSE until model selection is complete
#' @param real.ids vector of id values for which real parameters should be output with std error information for TMB models; if NULL all ids used
#' @param check if TRUE values of gamma, dmat and delta are checked to make sure the values are valid with initial parameter values.
#' @param useHess if TRUE, the TMB hessian function is used for optimization; using hessian is typically slower with many parameters but can result in a better solution
#' @param optimize if TRUE, optimizes to get parameter estimates; set to FALSE to extract estimates of ADREPORTed values only
#' @param vcv if TRUE, computes var-covariance matrix of ADREPORTed values
#' @param unit_scale default TRUE, if FALSE any time scaled parameter (e.g. Phi,S) is scaled when computing real value such that it represents the length of the interval rather than a unit interval
#' @param ... optional arguments passed to js or cjs and optimx
#' @importFrom graphics boxplot par
#' @importFrom stats as.formula binomial coef density
#'             glm.fit median model.frame model.matrix optim
#'              plogis pnorm predict rgamma rmultinom
#'              rnorm sd nlminb
#' @importFrom utils capture.output flush.console
#'             read.delim
#' @import data.table
#' @return crm model object with class=("crm",submodel)
#' @author Jeff Laake
#' @export crm
#' @import optimx Matrix numDeriv
#' @seealso \code{\link{make.design.data}},\code{\link{process.data}}
#' @keywords models
crm <- function(data,ddl=NULL,begin.time=1,model="MSCJS",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
 groups = NULL, time.intervals = NULL,debug=FALSE, method="nlminb", hessian=FALSE, accumulate=TRUE,chunk_size=1e7,
 control=list(),refit=1,itnmax=5000,scale=NULL,run=TRUE,compile=FALSE,extra.args=NULL,
 strata.labels=NULL,clean=FALSE,save.matrices=FALSE,savef=TRUE,getreals=FALSE,real.ids=NULL,check=FALSE,useHess=FALSE,optimize=TRUE,
 vcv=FALSE,unit_scale=TRUE,...)
{
model=toupper(model)
ptm=proc.time()
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
#
if(is.null(data$data))
{
   if(!is.null(ddl))
   {
      warning("Warning: specification of ddl ignored, as data have not been processed")
      ddl=NULL
   }
   if(debug)message("Model: ",model,"\n")
   if(debug)message("Processing data...\n")
   flush.console()
   data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=1,
	   groups = groups, age.var = NULL, initial.ages = NULL,
	   time.intervals = time.intervals,nocc=NULL,accumulate=accumulate,strata.labels=strata.labels)
}
else
{
	data.proc=data
	model=data$model
}
#
# Setup parameter list
#
number.of.groups=1
if(!is.null(data.proc$group.covariates))number.of.groups=nrow(data.proc$group.covariates)
par.list=setup.parameters(data.proc$model,check=TRUE)
#
# Check validity of parameter list; stop if not valid
#
if(!valid.parameters(model,model.parameters)) stop()
parameters=setup.parameters(data.proc$model,model.parameters,data$nocc,number.of.groups=number.of.groups)
parameters=parameters[par.list]
# See if any formula contain random effects and set re
re=FALSE
for (i in 1:length(parameters))
{
	if(is.null(parameters[[i]]$formula)) parameters[[i]]$formula=~1
	mlist=proc.form(parameters[[i]]$formula)

	if(!is.null(mlist$re.model))
	{
		re_names=sub("^\\s+", "",sapply(strsplit(names(mlist$re.model),"\\|"),function(x)x[2]))
		if(length(re_names)>1 | !"id" %in% re_names) crossed=TRUE
		if((length(re_names)> 1 || re_names[1]!="time" ) & any(data.proc$freq>1))
			stop("\n data cannot be accumulated (freq>1) except with temporal random effects only; set accumulate=FALSE\n")
		re=TRUE
	}
	if(parameters[[i]]$nointercept)parameters[[i]]$remove.intercept=TRUE
}
#
# If the design data have not been constructed, do so now
#
external.ddl=FALSE
if(is.null(ddl))
{
  if(debug)message("Creating design data...\n")
	flush.console()
	ddl=make.design.data(data.proc,design.parameters)
} else
{
#
# if ddl is external check to make sure ddl.rda exists and set logical
#
  if(is.character(ddl)&&toupper(ddl)=="EXTERNAL")
  {
    external.ddl=TRUE
    if(file.exists("ddl.rda"))
    {
      load(file="ddl.rda")
      if(!exists("ddl")) stop("\nexternal ddl.rda file must contain object named ddl")
    } else
    {
      stop("\nCannot find external file named ddl.rda")
    }
  }
#
# check to make sure ddl dataframes are in order
#
	for (i in 1:length(parameters))
	{
		if(!is.null(ddl[[i]]$order))
		   if(any(ddl[[i]]$order!=1:nrow(ddl[[i]])))
			   stop(paste("Design data for parameter",names(parameters)[i],"is out of order."))
	}
    if(!is.null(design.parameters))
		  for(parname in names(design.parameters))
			   ddl$design.parameters[[parname]]=c(ddl$design.parameters[[parname]],design.parameters[[parname]])
	  design.parameters=ddl$design.parameters
}
#
#   setup fixed values if old way used in parameter specification
#
ddl=set.fixed(ddl,parameters)
#
# check to make sure ddl values are setup correctly for each parameter
#
if(check) check_ddl_values(model,ddl)
#
# simplify the ddl
if(model=="SMSLD")
{
  ddl$S=ddl$S[ddl$S$Time>=ddl$S$Cohort,]
  ddl$r=ddl$r[ddl$r$Time>=ddl$r$Cohort,]
  ddl$p=ddl$p[ddl$p$Time>=ddl$p$Cohort,]
  ddl$Psi=ddl$Psi[ddl$Psi$Time>=ddl$Psi$Cohort,]
}
fullddl=ddl
if(debug)message("Simplifying design data\n")
ddl=simplify_ddl(ddl,parameters) # add indices to ddl and reduce ddl to unique values used
#
# check to see if all values for a parameter have been fixed.  If so, then set formula to ~0
#
for (i in 1:length(parameters))
{
	if(!is.null(ddl[[i]]$fix))
	{
		if(!is.null(ddl[[i]]$fix) && all(!is.na(ddl[[i]]$fix)))
		{
			message(paste("All values for",names(parameters)[i],"have been fixed. Setting formula to ~0\n"))
			parameters[[i]]$formula=~0
		} else {
			if(parameters[[i]]$formula==~0)
				stop(paste("Cannot use formula ~0 for",names(parameters)[i],"when some of the parameters must be estimated.\n"))
		}
	} else
	   if(parameters[[i]]$formula==~0)
		   stop(paste("Cannot use formula ~0 for",names(parameters)[i],"when some of the parameters must be estimated.\n"))
}
# Create design matrices for each parameter
if(debug)message("Creating design matrices\n")
dml=create.dml(ddl,model.parameters=parameters,design.parameters=design.parameters,chunk_size=chunk_size)
# store fulldml
fulldml=dml
for(parx in names(dml))
{
  fulldml[[parx]]$fe=dml[[parx]]$fe[ddl[[paste(parx,".indices",sep="")]],,drop=FALSE]
  #parameters[[parx]]$indices=ddl[[paste(parx,".indices",sep="")]]
}
# NEED TO REVISIT THIS
# For HMM call set.initial to get ptype and set initial values
if(nchar(model)>=4 &substr(model,1,4)=="MVMS")
	initial.list=set.initial(names(dml),dml,initial)
else
	initial.list=NULL
# if not running, return object with data,ddl,dml etc
if(!run) return(list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,ddl=ddl,dml=dml,results=initial.list))
#
# Depending on method set some values
#
if("SANN"%in%method)
{
	if(length(method)>1)
		warning("***SANN can only be used by itself; other methods ignored.")
	method="SANN"
  control$maxit=itnmax
}
if("nlminb"%in%method)
{
	control$eval.max=itnmax
	control$iter.max=itnmax
}
#
# Call estimation function which depends on the model
#
if(debug)message("Fitting model\n")
#
# MSCJS model
#
if(model=="MSCJS")
		runmodel=mscjs_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,getreals=getreals,useHess=useHess,savef=savef,...)
#
# MSJS model
#
if(model=="MSJS")
  runmodel=msjs_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
                     refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,getreals=getreals,useHess=useHess,savef=savef,...)
#
# MSJSU model
#
if(model=="MSJSU")
  runmodel=msjsu_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
                    refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,getreals=getreals,useHess=useHess,savef=savef,...)
#
# MSLD model
#
if(model=="MSLD" || model=="SMSLD")
{
  save(fullddl,file="tmp.rda")
  rm(fullddl)
  if(model=="MSLD")
  runmodel=msld_tmb(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
                    refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,
                    clean=clean,getreals=getreals,useHess=useHess,savef=savef,...)
  else
    runmodel=smsld_tmb(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
                      refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,
                      clean=clean,getreals=getreals,useHess=useHess,savef=savef,...)
  load("tmp.rda")
}
#
# MVMS model
#
if(model=="MVMSCJS")
{
	sup=data.proc$fct_sup(list(obslevels=data.proc$ObsLevels))
	mx=data.proc$m
# call HMMlikelihood to check for problems in setup
#  if(check)
#  {
#    message("Checking for problems in design data setup\n")
#    xx=HMMLikelihood(par=unlist(initial.list$par),xx=data.proc$ehmat,mx=mx,
#                     type=initial.list$ptype,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
#                     fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=fullddl,dml=fulldml,
#                     parameters=parameters,sup=sup,check=TRUE)
#  }
  runmodel=mvmscjs_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
	                   refit=refit,control=control,itnmax=itnmax,re=FALSE,compile=compile,clean=clean,sup=sup,getreals=getreals,real.ids=real.ids,useHess=useHess,optimize=optimize,
	                   vcv=vcv,savef=savef,...)
	par=coef(runmodel)[,1]
	runmodel$options=c(runmodel$options,list(accumulate=accumulate,initial=initial.list$par,method=method,
			chunk_size=chunk_size,itnmax=itnmax,control=control))
	parlist=split(par,initial.list$ptype)
	par=vector("list",length=length(names(initial.list$par)))
	names(par)=names(initial.list$par)
	for(p in names(parlist))
	{
		par[[p]]=parlist[[p]]
		names(par[[p]])=colnames(dml[[p]]$fe)
	}
	runmodel$beta=par
	runmodel$par=NULL
	if(is.null(runmodel$neg2lnl))
		runmodel$neg2lnl=2*runmodel$optim.details$value
	runmodel$AIC=runmodel$neg2lnl+2*sum(sapply(runmodel$beta,length))
	if(!is.null(runmodel$hessian))
	{
		runmodel$beta.vcv=solvecov(runmodel$hessian)$inv
		colnames(runmodel$beta.vcv)=names(unlist(runmodel$beta))
		rownames(runmodel$beta.vcv)=colnames(runmodel$beta.vcv)
	}
	class(runmodel)=c("crm","mle",model)
}
# Return fitted MARK model object or if external, return character string with same class and save file
if(!is.null(runmodel$convergence) && runmodel$convergence!=0)
{
	warning("******Model did not converge******")
	msg=attr(runmodel$optim.details,"details")$message
	if(is.null(msg)) msg="Exceeded maximum number of iterations"
	warning(msg)
}
object=list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,results=runmodel)
class(object)=class(runmodel)
#
# if save.matrices call here and store
#
if(save.matrices)object$results$mat=compute_matrices(object)
#
if(getreals)object$results$reals=predict(object,ddl=fullddl,unique=TRUE,se=hessian,unit_scale=unit_scale)
if(file.exists("tmp.rda"))unlink("tmp.rda")
message(paste("\nElapsed time in minutes: ",round((proc.time()[3]-ptm[3])/60,digits=4),"\n"))
return(object)
}
# solvecov code was taken from package fpc: Christian
# Hennig <chrish@@stats.ucl.ac.uk> \url{http://www.homepages.ucl.ac.uk/~ucakche/}
solvecov=function (m, cmax = 1e+10)
# from package fpc
{
	options(show.error.messages = FALSE)
	covinv <- try(solve(m))
	if (!is(covinv,"try-error"))
		coll = FALSE
	else {
		p <- nrow(m)
		cove <- eigen(m, symmetric = TRUE)
		coll <- TRUE
		if (min(cove$values) < 1/cmax) {
			covewi <- diag(p)
			for (i in 1:p) if (cove$values[i] < 1/cmax)
					covewi[i, i] <- cmax
				else covewi[i, i] <- 1/cove$values[i]
		}
		else covewi <- diag(1/cove$values, nrow = length(cove$values))
		covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
	}
	options(show.error.messages = TRUE)
	out <- list(inv = covinv, coll = coll)
	out
}
