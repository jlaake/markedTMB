check_ddl_values=function(model,ddl)
{
if(substr(model,1,4)=="MVMS")
{
  check_mlogit=function(x){ifelse(sum(x,na.rm=TRUE)==0||(any(is.na(x))&!(any(x[!is.na(x)]==1))),TRUE,FALSE)}
  if(is.null(ddl$pi$fix))
    message("\n Warning: No values provided for fix for pi. Must have a reference cell via formula.")
  else
  {
    bad_pi=sapply(split(ddl$pi$fix,ddl$pi$id),check_mlogit)
    if(any(bad_pi))
    {
      message("\n Check values of fix for pi with id:")
      cat(names(which(bad_pi)))
    }
  }
  if(is.null(ddl$delta$fix))
  {
    message("\n Warning: No values provided for fix for delta. Must have a reference cell via formula.")
  }else
  {
    xx=list(id=ddl$delta$id,occ=ddl$delta$occ,stratum=ddl$delta$stratum)
    bad_delta=sapply(split(ddl$delta$fix,xx),check_mlogit)
    if(any(bad_delta))
    {
      message("\n Warning: Check values of fix for delta for the following records with id.occ.stratum.")
      cat(names(which(bad_delta)))
    }
  }
  if(is.null(ddl$Psi$fix))
  {
    message("\n Warning: No values provided for fix for Psi. Must have a reference cell via formula.")
  }else
  {
    xx=list(id=ddl$Psi$id,occ=ddl$Psi$occ,stratum=ddl$Psi$stratum)
    bad_Psi=sapply(split(ddl$Psi$fix,xx),check_mlogit)
    if(any(bad_Psi))
    {
      message("\n Warning: Check values of fix for Psi for the following records with id.occ.stratum.")
      cat(names(which(bad_Psi)))
    }
  }
}
return(NULL)
}
