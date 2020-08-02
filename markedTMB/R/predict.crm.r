#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular subset of 
#' parameters. You can compute real estimates for a subset of values or a new set of values is to specify a limited 
#' range of the values in ddl for each parameter. Make sure to include a complete set of values that spans
#' the factor levels and individual covariates used in the formulas for the model object or you will receive an
#' error that the number of columns in the design matrix does not match the number of beta parameters.  You cannot 
#' change the levels of any factor variable or modify the design data in anyway that changes the design matrix.
#' 
#' If the real estimates are in the model object and se and vcv are FALSE and ddl not specified, 
#' the code will simply pull the values from the model object.
#' 
#' 
#' @usage \method{predict}{crm}(object,ddl=NULL,parameter=NULL,unique=TRUE,
#'                    vcv=FALSE,se=FALSE,chat=1,subset=NULL,select=NULL,
#'                    real.ids=NULL,merge=FALSE,unit_scale=TRUE,...)
#' @param object model object;
#' @param ddl list of simplified dataframes for design data
#' @param ddl list of dataframes for design data (not simplified)
#' @param parameter name of real parameter to be computed (eg "Phi")
#' @param unique TRUE if only unique values should be returned; if TRUE forces computation even if estimates are in model object
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates; if TRUE forces computation even if estimates are in model object
#' @param se logical; if TRUE, computes std errors and conf itervals of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include
#' @param real.ids animal ids passed to TMB code for computation of real parameter values
#' @param merge default FALSE but if TRUE, the ddl for the parameter is merged (cbind) to the estimates
#' @param unit_scale default TRUE, if FALSE any time scaled parameter (e.g. Phi,S) is scaled when computing real value such that it represents the length of the interval rather than a unit interval
#' @param ... generic arguments not used here
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if vcv=TRUE it also contains
#' standard errors and confidence intervals} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @export
#' @keywords utility
predict.crm <-function(object,ddl=NULL,parameter=NULL,unique=TRUE,vcv=FALSE,se=FALSE,
                       chat=1,subset=NULL,select=NULL,real.ids=NULL,merge=FALSE,unit_scale=TRUE,...)
{
  if(is.null(ddl))
  {
    if(!(se||vcv))
    {
      if(!is.null(object$results$reals))
          return(object$results$reals)
      else
         stop("No real values available in object and ddl not specified")
    }
    else
    {
      stop("No ddl specified")
    }
  }
  if(is.null(parameter))
  {
    results=NULL
    for (parameter in names(object$model.parameters))
      results[[parameter]]=compute_real(object,parameter,ddl,unique,vcv,se,chat,subset=substitute(subset),select,include=object$model.parameters[[parameter]]$include,merge=merge,unit_scale=unit_scale)
    return(results)
  } else
    return(compute_real(object,parameter,ddl,unique,vcv,se,chat,subset=substitute(subset),select,include=object$model.parameters[[parameter]]$include,merge=merge,unit_scale=unit_scale))	
}
