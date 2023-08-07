#' Defines model specific parameters (internal use)
#'
#' Compares \code{model}, the name of the type of model (eg "CJS") to the list
#' of acceptable models to determine if it is supported and then creates some
#' global fields specific to that type of model that are used to modify the
#' operation of the code.
#'
#' In general, the structure of the different types of models (e.g.,
#' "CJS","Recovery",...etc) are very similar with some minor exceptions.  This
#' function is not intended to be called directly by the user but it is
#' documented to enable other models to be added.  This function is called by
#' other functions to validate and setup model specific parameters.  For
#' example, for live/dead models, the length of the capture history is twice
#' the number of capture occasions and the number of time intervals equals the
#' number of capture occasions because the final interval is included with dead
#' recoveries.  Whereas, for recapture models, the length of the capture
#' history is the number of capture occasions and the number of time intervals
#' is 1 less than the number of occasions.  This function validates that the
#' model is valid and sets up some parameters specific to the model that are
#' used in the code.
#'
#' @param model name of model type (must be in vector \code{valid.models})
#' @param nocc length of capture history string
#' @param mixtures number of mixtures
#' @export
#' @aliases setup.model setupHMM
#' @return model.list - a list with following elements \item{etype}{encounter
#' type string for MARK input; typically same as model} \item{nocc}{number of
#' capture occasions} \item{num}{number of time intervals relative to number of
#' occasions (0 or -1)} \item{mixtures}{number of mixtures if any}
#' \item{derived}{logical; TRUE if model produces derived estimates}
#' @author Jeff Laake
#' @seealso \code{\link{setup.parameters}}, \code{\link{valid.parameters}}
#' @keywords utility
setup.model <-
		function(model,nocc,mixtures=1)
{
    # Read in parameter definitions
	fdir=system.file(package="markedTMB")
	fdir=file.path(fdir,"models.txt")
	model_definitions=read.delim(fdir,header=TRUE,
			colClasses=c("character",rep("numeric",1),rep("logical",4),rep("numeric",1)))
	model_def=model_definitions[model_definitions$model==model,]
	if(nrow(model_def)==0)
		stop("Invalid type of model = ",model," Valid types are\n", paste(model_definitions$model,collapse="\n"))
    # model_def$nocc=nocc/model_def$divisor; not used at present
	model_def$nocc=nocc
	model_def=as.list(model_def)
	return(model_def)
}
setupHMM=function(model_def,model,strata.labels)
{
  if(toupper(model)=="MSLD")
  {
    model_def$hmm$strata.labels=strata.labels
    model_def$hmm$m=2*length(strata.labels)+1
    model_def$hmm$ObsLevels=c("0",strata.labels,"1")
  }
  if(toupper(model)=="MSCJS")
  {
    model_def$hmm$strata.labels=strata.labels
    model_def$hmm$m=length(strata.labels)+1
    model_def$hmm$ObsLevels=c(0,strata.labels)
  }
  if(toupper(model)=="MSJS")
  {
    model_def$hmm$strata.labels=c("N",strata.labels)
    model_def$hmm$m=length(strata.labels)+2
    model_def$hmm$ObsLevels=c(0,model_def$hmm$strata.labels)
  }
  if(toupper(model)=="MSJSU")
  {
    model_def$hmm$strata.labels=c("N",strata.labels)
    model_def$hmm$m=length(strata.labels)+2
    model_def$hmm$ObsLevels=c(0,model_def$hmm$strata.labels,"U")
  }
  if(toupper(model)=="MVMSCJS")
	{
		model_def$hmm$fct_sup=mvms_sup
		model_def$hmm$strata.list=set_mvms(strata.labels)
		model_def$hmm$strata.labels=apply(model_def$hmm$strata.list$df.states,1,paste,collapse="")
		model_def$hmm$m=nrow(model_def$hmm$strata.list$df.states)+1
		model_def$hmm$ObsLevels=c(0,apply(model_def$hmm$strata.list$df,1,paste,collapse=""))
	}
	if(toupper(model)=="MVMS")
	{
		model_def$hmm$strata.list=set_mvms(strata.labels)
		model_def$hmm$strata.labels=apply(model_def$hmm$strata.list$df.states,1,paste,collapse="")
		model_def$hmm$m=nrow(model_def$hmm$strata.list$df.states)+1
		model_def$hmm$ObsLevels=c(0,apply(model_def$hmm$strata.list$df,1,paste,collapse=""))
	}
	if(toupper(model)=="ATTEND")
	{
		model_def$hmm$strata.labels=strata.labels
		model_def$hmm$m=length(strata.labels)+1
		model_def$hmm$ObsLevels=c(0,1)
	}
	return(model_def)
}
mvms_sup=function(x) 
{
  #   supplemental function to provide information to dmat function 
  #   it is only run once per model fit because it doesn't change
  obslevels=x$obslevels
  #   Function to identify indices that should be filled in matrix
  identify=function(x,y) { 
    col=grep(y,s)
    if(length(col)>0)
      return(data.frame(row=x,col=col))
    else
      return(NULL)
  }
  #   state names
  unknown=grep("u",obslevels,fixed=TRUE)
  if(length(unknown)==0)
    s=obslevels[-1]
  else
    s=obslevels[-grep("u",obslevels,fixed=TRUE)][-1]
  #   indices for p and delta for non-zero values
  indices_forp=as.matrix(do.call("rbind",mapply(identify,1:length(obslevels),gsub("\\+","\\\\+",gsub("u",".",obslevels)))))
  s=c(s,"Dead")
  np=nrow(indices_forp)
  pcounts=table(indices_forp[,2])
  indices_forp=indices_forp[order(indices_forp[,2]),]
  return(list(indices_forp=indices_forp,np=np,pcounts=pcounts,obslevels=obslevels))
}



