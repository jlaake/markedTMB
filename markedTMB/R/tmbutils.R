#' TMB setup
#'
#' Sets up executable for the .cpp file (tpl) by looking for exe in package directory or
#' compiles cpp file in local directory (clean=FALSE) of from package directory.
#'
#' @param tpl character string for admb template file
#' @param clean if TRUE, deletes the tpl (.cpp) and executable files in local directory and copies from package directory
#' @param debug if TRUE, writes out debug line values in code
#' @import optimx 
#' @export
setup_tmb=function(tpl,clean=FALSE,debug=FALSE)
{
  if(debug)message("TMB program setup\n")
  sdir=system.file(package="markedTMB")
  #   if argument clean is TRUE, delete dll and cpp files as well
  if(clean)
  {
    if(debug)message("deleting old TMB program and rebuilding\n")
    if(file.exists(paste(tpl,".cpp",sep=""))) unlink(paste(tpl,".cpp",sep=""))
    if(tpl%in%names(getLoadedDLLs()))dyn.unload(dynlib(tpl))
    if(file.exists(dynlib(tpl))) unlink(dynlib(tpl))
    file.copy(file.path(sdir,paste(tpl,".cpp",sep="")),file.path(getwd(),paste(tpl,".cpp",sep="")),overwrite=TRUE)
    if(debug)message("compiling and linking TMB program\n")
    compile(paste(tpl,".cpp",sep=""), flags = '-Wno-ignored-attributes')    # Compile the C++ file
    dyn.load(dynlib(tpl))                           # Dynamically link the C++ code
  } else
  {
    # if dll is available load it
    if(file.exists(dynlib(tpl)))
    {
      if(debug)message("loading existing TMB program\n")
      if(!tpl%in%names(getLoadedDLLs())) {
        dyn.load(dynlib(tpl))
      }
    } else
      # check for cpp file in directory; if it doesn't exist then copy from package directory
      if(!file.exists(paste(tpl,".cpp",sep="")))
      {
        file.copy(file.path(sdir,paste(tpl,".cpp",sep="")),file.path(getwd(),paste(tpl,".cpp",sep="")),overwrite=TRUE)
        if(tpl%in%names(getLoadedDLLs()))dyn.unload(dynlib(tpl))
        if(file.exists(dynlib(tpl))) unlink(dynlib(tpl))
        if(debug)message("compiling and linking TMB program\n")
        compile(paste(tpl,".cpp",sep=""), flags = '-Wno-ignored-attributes')    # Compile the C++ file
        dyn.load(dynlib(tpl))          # Dynamically link the C++ code
      } else
      {
        if(file.exists(paste(tpl,".o",sep=""))) unlink(paste(tpl,".o",sep=""))
        if(debug)message("compiling and linking TMB program\n")
        compile(paste(tpl,".cpp",sep=""), flags = '-Wno-ignored-attributes')    # Compile the C++ file
        if(is.loaded(dynlib(tpl)))dyn.unload(dynlib(tpl))
        dyn.load(dynlib(tpl))          # Dynamically link the C++ code
      }
  }
  invisible()
}
#' Load library for TMB model 
#' 
#' Loads the model specific library (eg dll) for a model
#' 
#' @param model type of crm model 
#' @return NULL
#' @export load_lib
#' @keywords utility
load_lib <- function(model)
{  	
  mod_list <- c("MSCJS","MVMSCJS","MSLD")
  if(model%in%mod_list){
    if(!"markedTMB_TMBExports"%in%names(getLoadedDLLs())){
      fff <- paste0(system.file(package='markedTMB'), "/libs/markedTMB_TMBExports")
      dyn.load(dynlib(fff))
    }
  } else {
    # Use the following template to test new models. Once development is complete
    # the new model should be incorporated into the `mod_list` above
    # if(model=="MSCJS"){
    #   if(!"multistate_tmb"%in%names(getLoadedDLLs())) dyn.load(dynlib("multistate_tmb"))
    # }
  }
  return(NULL)
}

