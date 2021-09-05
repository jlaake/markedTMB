#' Parametric bootstrap of state predictions
#' 
#' Bootstraps the state predictions for a model by assuming a multivariate normal distribution for
#' the model parameter estimates. For each bootstrap replicate it replaces model parameter estimates with a 
#' random sample from a multivariate normal where the mean vector is the fitted model parameter estimates and
#' the var-cov matrix is the estimated matrix from the model fit. 
#' 
#' @param model fitted model object
#' @param nboot number of bootstrap replicates; if 0 it returns the model state predictions with global_decode at the estimated parameter values
#' @param beta_only if TRUE, generates beta_reps from model and returns it. This is for cases where the same values of beta_reps are used for another bootstrap (eg projections using survival)
#' @param beta_reps  if NULL the values are generated; otherwise it is a matrix of values with rows being the reps
#' @export
#' @importFrom mvtnorm rmvnorm
#' @return If nboot=0, a matrix of model state predictions where rows are the individual histories and columns are the occasions.
#' If nboot>0, a list of length nboot with state prediction matrices is returned.
#' @author Jeff Laake
#' @keywords utility
boot_pred_states=function(model,nboot=0,beta_only=FALSE,beta_reps=NULL)
{
  if(nboot==0 &is.null(beta_reps))
  {
    return(global_decode(model))
  } else
  {
    if(is.null(beta_reps))
    {
      # get beta (estimates) and vcv matrix from model
      beta=unlist(model$results$beta)
      beta.vcv=model$results$beta.vcv
      # get multivariate normal random beta values
      beta_reps<-rmvnorm(nboot,mean=beta,sigma=beta.vcv)
    } else
    {
      nboot=nrow(beta_reps)
      if(ncol(beta_reps)!=length(unlist(model$results$beta)))stop("beta_reps doesn't match the model")
    }
    if(beta_only)
      return(beta_reps)
    pred_states_list=vector("list",length=nboot)
    # compute state predictions for each bootstrap replicate
    for(i in 1:nboot)
    {
      cat("\nrep = ",i)
      pred_states_list[[i]]=global_decode(model,par=beta_reps[i,])
    }
    return(pred_states_list)
  }
}

#' Summarize statistic values across parametric bootstraps of state predictions
#' 
#' From bootstraps replicates of the state predictions, compute the value of a statistic summarizing some
#' attribute. Computes the statistic value for the fitted model, standard error of the statistic value (std dev across bootstrap replicates),
#' lower and upper 95% confidence intervals assuming an asymptotic normal distribution and percentile 
#' confidence intervals across replicates.
#' 
#' @param pp state prediction matrix at fitted model estimates
#' @param bs list of bootstrap replicates of state prediction matrices
#' @param f function to compute some desired value for each state prediction matrix
#' @param index if NULL, f only produces a single value; otherwise index specifies which of the values should be summarized
#' @export
#' @return A dataframe containing the estimate, std error, 95% normal confidence interval and 95% percentile intervals.
#' @author Jeff Laake
#' @keywords utility
Values=function(pp,bs,f,index=NULL)
{
  if(is.null(index))
  {
    estimate=f(pp)
    estimates=sapply(bs,f)
    se=sqrt(var(estimates))
  } else
  {
    estimate=f(pp)[[index]]
    estimates=unlist(sapply(bs,f)[index,])
    se=sqrt(var(estimates))
  }
  nboot=length(bs)
  lcl.pct=sort(estimates)[ceiling(nboot*.025)]
  ucl.pct=sort(estimates)[floor(nboot*.975)]
  return(data.frame(estimate=estimate,se=se,lcl=estimate-1.96*se,ucl=estimate+1.96*se,lcl.pct=lcl.pct,ucl.pct=ucl.pct))
}
