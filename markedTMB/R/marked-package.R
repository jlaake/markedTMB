#' Dipper capture-recapture data
#'
#' A capture-recapture data set on European dippers from France that
#' accompanies MARK as an example analysis using the CJS and POPAN models.  The
#' dipper data set was orginally described as an example by Lebreton et al
#' (1992).
#'
#' There isn't a specific CJS model in markedTMB but you can use the model MSCJS with
#' a dummy state that is never observed.  This example shows how that can be done to 
#' fit standard JS models.
#'
#'
#' Note that the covariate "sex" defined in dipper has values "Male" and
#' "Female".  It cannot be used directly in a formula for MARK without using it
#' do define groups because MARK.EXE will be unable to read in a covariate with
#' non-numeric values.  By using \code{groups="sex"} in the call the
#' \code{\link{process.data}} a factor "sex" field is created that can be used
#' in the formula.  Alternatively, a new covariate could be defined in the data
#' with say values 0 for Female and 1 for Male and this could be used without
#' defining groups because it is numeric.  This can be done easily by
#' translating the values of the coded variables to a numeric variable.  Factor
#' variables are numbered 1..k for k levels in alphabetic order.  Since Female
#' < Male in alphabetic order then it is level 1 and Male is level 2.  So the
#' following will create a numeric sex covariate.
#'
#' \preformatted{ dipper$numeric.sex=as.numeric(dipper$sex)-1 }
#'
#' @name dipper
#' @docType data
#' @format A data frame with 294 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird} \item{sex}{the sex of the bird: a factor with levels
#' \code{Female} \code{Male}} }
#' @source Lebreton, J.-D., K. P. Burnham, J. Clobert, and D. R. Anderson.
#' 1992. Modeling survival and testing biological hypotheses using marked
#' animals: case studies and recent advances. Ecol. Monogr. 62:67-118.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(dipper)
#' ms1=process.data(dipper,model="MSCJS",strata.labels=c("1","2"))
#' ms1.ddl=make.design.data(ms1)
#' ms1.ddl$Psi$fix[ms1.ddl$Psi$stratum!=ms1.ddl$Psi$tostratum]=0
#' ms1.ddl$p$fix=NA
#' ms1.ddl$p$fix[ms1.ddl$p$stratum==2]=0
#' mod_cjs=crm(ms1,ms1.ddl,model.parameters=list(S=list(formula=~1),
#'                                                p=list(formula=~time)),hessian=TRUE)
#' hmm_mat=compute_matrices(mod_cjs,state.names=c("Alive","NULL","Dead"))
#' }                                                

NULL


#' Multistrata example data
#'
#' An example data set which appears to be simulated data that accompanies MARK
#' as an example analysis using the Multistrata model.
#'
#' This is a data set that accompanies program MARK as an example for the
#' Multistrata model and is also in the RMark pacakge. Here I use it to show the
#' 3 ways models can be fitted to multistrata data. The model MSCJS is not run because it
#' requires ADMB or the exe constructed from ADMB which is not available if downloaded from CRAN.
#'
#' @name mstrata
#' @docType data
#' @format A data frame with 255 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird with strata} \item{freq}{the number of birds with that capture
#' history} }
#' @keywords datasets
#' @examples
#' \donttest{
#' data(mstrata)
#' ms1=process.data(mstrata,model="MSCJS",strata.labels=c("A","B","C"))
#' ms1.ddl=make.design.data(ms1)
#' mod_mscjs=crm(ms1,ms1.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#'                                                p=list(formula=~time)),hessian=TRUE,
#'                                                save.matrices=TRUE)
#' mod_mscjs
#' mod_mscjsr=crm(ms1,ms1.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#'                                                p=list(formula=~1+(1|time))),hessian=TRUE)
#' mod_mscjsr
#' # strata.labels for MVMS models must be specified as a list because more than one variable
#' # can be used
#' ms2=process.data(mstrata,model="MVMSCJS",strata.labels=list(state=c("A","B","C")))
#' ms2.ddl=make.design.data(ms2)
#' ms2.ddl$delta$fix=1
#' # uses TMB with mvmscjs
#' mod_mvmscjs=crm(ms2,ms2.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#'                                            p=list(formula=~time)),hessian=TRUE)
#' mod_mvmscjs
#' }
NULL

#' Multivariate State example data
#'
#' An example data set using California sea lions to demonstrate the Multivariate
#' state model with 3 variables defining the states : area, ltag and rtag.  The left tag (ltag) and
#' right tag (rtag) variables can be unknown.
#' #'
#' @name sealions
#' @docType data
#' @format A data frame with 485 observations on the following 3 variables.
#' \describe{ \item{ch}{a character vector of 3-character values for each occasion separated by commas.}
#'  \item{sex}{the sex of the sea lion designated as F or M}
#'  \item{weight}{the anomaly of the weight from the sex-specific mean}
#' }
#' @keywords datasets
#' @examples
#' \donttest{
#' ### Load packages ###
#' # The splines package is only necessary for fitting b-spline curves used in the paper
#' # It is not required for the multivate models in the marked package
#' library(splines)
#'
#' # Get data
#' data(sealions)
#'
#' # Process data for multivariate models in marked
#' dp=process.data(sealions,model="mvmscjs",
#'   strata.labels=list(area=c("A","S"),ltag=c("+","-","u"),rtag=c("+","-","u")))
#'
#' ### Make design data
#' ddl=make.design.data(dp)
#'
#' # Create pup variable for Phi
#' ddl$Phi$pup=ifelse(ddl$Phi$Age==0,1,0)
#' ddl$Phi$sex=factor(ddl$Phi$sex)
#'
#' # Detection model
#' # Set final year (2014)  p=0 (no resight data) for ANI
#' ddl$p$fix = ifelse(ddl$p$Time==17 & ddl$p$area=="A", 0, ddl$p$fix)
#'
#' # Delta model
#' # create indicator variables for 'unknown' tag observations
#' ddl$delta$obs.ltag.u = ifelse(ddl$delta$obs.ltag=="u", 1, 0)
#' ddl$delta$obs.rtag.u = ifelse(ddl$delta$obs.rtag=="u", 1, 0)
#'
#' # Psi model
#' # Set Psi to 0 for cases which are not possible - missing tag to having tag
#' ddl$Psi$fix[as.character(ddl$Psi$ltag)=="-"&as.character(ddl$Psi$toltag)=="+"]=0
#' ddl$Psi$fix[as.character(ddl$Psi$rtag)=="-"&as.character(ddl$Psi$tortag)=="+"]=0
#' # Create indicator variables for transitioning between states
#' ddl$Psi$AtoS=ifelse(ddl$Psi$area=="A"&ddl$Psi$toarea=="S",1,0)  # ANI to SMI movement
#' ddl$Psi$StoA=ifelse(ddl$Psi$area=="S"&ddl$Psi$toarea=="A",1,0)  # SMI to ANI movement
#' ddl$Psi$lpm=ifelse(ddl$Psi$ltag=="+"&ddl$Psi$toltag=="-",1,0)   # Losing left tag
#' ddl$Psi$rpm=ifelse(ddl$Psi$rtag=="+"&ddl$Psi$tortag=="-",1,0)   # Losing right tag
#' ddl$Psi$sex=factor(ddl$Psi$sex)
#'
#' # formulas
#' Psi.1=list(formula=~-1+ AtoS:sex + AtoS:sex:bs(Age) + StoA:sex + StoA:sex:bs(Age) +
#'                      I(lpm+rpm) +I(lpm+rpm):Age + lpm:rpm)
#' p.1=list(formula=~-1+time:area)
#' delta.1=list(formula= ~ -1 + obs.ltag.u + obs.rtag.u + obs.ltag.u:obs.rtag.u)
#' Phi.1=list(formula=~sex*bs(Age)+pup:weight+area)
#'
#' # Fit model with TMB
#'  mod_sealions=crm(dp,ddl,model.parameters=list(Psi=Psi.1,p=p.1,delta=delta.1,Phi=Phi.1),
#'  method="nlminb",hessian=TRUE,save.matrices=TRUE)
#' }
NULL

#' An example of the Mulstistrata (multi-state) model in which states are routes taken by migrating fish.
#'
#' @name skagit
#' @docType data
#' @format A data frame with 100 observations on the following 2 variables.
#' \describe{ \item{ch}{capture history}
#'  \item{tag}{tag type v7 or v9}
#' }
#' @keywords datasets
#' @author Megan Moore <megan.moore at noaa.gov>
#' @examples
#'# There are just two states which correspond to route A and route B. There are 6 occasions
#'# which are the locations rather than times. After release at 1=A there is no movement
#'# between states for the first segment, fish are migrating downriver together and all pass 2A.
#'# Then after occasion 2, migrants go down the North Fork (3A) or the South Fork (3B),
#'# which both empty into Skagit Bay. Once in saltwater, they can go north to Deception Pass (4A)
#'# or South to a receiver array exiting South Skagit Bay (4B). Fish in route A can then only go
#'# to the Strait of Juan de Fuca, while fish in route B must pass by Admiralty Inlet (5B).
#'# Then both routes end with the array at the Strait of Juan de Fuca.
#'#
#'#       1A
#'#        |
#'#        2A
#'#      /     \
#'#    3A        3B
#'#   /  \      /  \
#'# 4A   4B    4A  4B
#'#  |     \    /   |
#'#   5A    5B  5A   5B
#'#      \   \   /    /
#'#            6
#'#
#'# from 3A and 3B they can branch to either 4A or 4B; branches merge at 6
#'# 5A does not exist so p=0; only survival from 4A to 6 can be
#'# estimated which is done by setting survival from 4A to 5A to 1 and
#'# estimating survival from 5A to 6 which is then total survival from 4A to 6.
#'
#'# See help for mscjs_tmb for an example that explains difference between marked and RMark
#'# with regard to treatment of mlogit parameters like Psi.
NULL

#'  Mulstistate Live-Dead Paradise Shelduck Data
#'
#' @name Paradise_shelduck
#' @docType data
#' @aliases ps
#' @description Paradise shelduck recapture and recovery data in multistrata provided by Richard Barker and Gary White.
#' @format  A data frame with 1704 observations of 3 variables
#'  \describe{
#'  \item{ch}{a character vector containing the capture history (each is 2 character positions LD) for 6 occasions}
#'  \item{freq}{capture history frequency}
#'  \item{sex}{Male or Female}
#'  }
#' @keywords datasets
#' @author Jeff Laake
#' @references Barker, R.J, White,G.C, and M. McDougall. 2005. MOVEMENT OF PARADISE SHELDUCK BETWEEN MOLT SITES:
#' A JOINT MULTISTATE-DEAD RECOVERY MARK–RECAPTURE MODEL. JOURNAL OF WILDLIFE MANAGEMENT 69(3):1194–1201.
#' @examples
#' \donttest{
#' # In the referenced article, there are 3 observable strata (A,B,C) and 3 unobservable
#' # strata (D,E,F). This example is setup by default to use only the 3 observable strata
#' # to avoid problems with multiple modes in the likelihood.
#' # Code that uses all 6 strata are provided but commented out.
#' # With unobservable strata, simulated annealing should
#' # be used (options="SIMANNEAL")
#' data("Paradise_shelduck")
#' # change sex reference level to Male to match design matrix used in MARK
#' ps$sex=relevel(ps$sex,"Male")
#' # Process data with MSLiveDead model using sex groups and specify only observable strata
#' ps_dp=process.data(ps,model="MSLD",groups="sex",strata.labels=c("A","B","C"))
#' # Process data with MSLiveDead model using sex groups and specify observable and
#' # unboservable strata
#' # ps_dp=process.data(ps,model="MSLD",groups="sex",strata.labels=c("A","B","C","D","E","F"))
#' # Make design data and specify constant PIM for Psi to reduce parameter space. No time
#' #variation was allowed in Psi in the article.
#' ddl=make.design.data(ps_dp)
#' # Fix p to 0 for unobservable strata (only needed if they are included)
#' ddl$p$fix=NA
#' ddl$p$fix[ddl$p$stratum%in%c("D","E","F")]=0
#' # Fix p to 0 for last occasion
#' ddl$p$fix[ddl$p$time%in%6:7]=0.0
#' # Fix survival to 0.5 for last interval to match MARK file (to avoid confounding)
#' ddl$S$fix=NA
#' ddl$S$fix[ddl$S$time==6]=0.5
#' # create site variable for survival which matches A with D, B with E and C with F
#' ddl$S$site="A"
#' ddl$S$site[ddl$S$stratum%in%c("B","C")]=as.character(ddl$S$stratum[ddl$S$stratum%in%c("B","C")])
#' ddl$S$site[ddl$S$stratum%in%c("E")]="B"
#' ddl$S$site[ddl$S$stratum%in%c("F")]="C"
#' ddl$S$site=as.factor(ddl$S$site)
#' # create same site variable for recovery probability (r)
#' ddl$r$site="A"
#' ddl$r$site[ddl$r$stratum%in%c("B","C")]=as.character(ddl$r$stratum[ddl$r$stratum%in%c("B","C")])
#' ddl$r$site[ddl$r$stratum%in%c("E")]="B"
#' ddl$r$site[ddl$r$stratum%in%c("F")]="C"
#' ddl$r$site=as.factor(ddl$r$site)
#' # Specify formula used in MARK model
#' S.1=list(formula=~-1+sex+time+site)
#' p.1=list(formula=~-1+stratum:time)
#' r.1=list(formula=~-1+time+sex+site)
#' Psi.1=list(formula=~-1+stratum:tostratum)
#' # Run top model from paper but only for observable strata
#' mod_msld=crm(ps_dp,ddl,model.parameters=list(S=S.1,p=p.1,r=r.1,Psi=Psi.1),
#'               method="nlminb",hessian=TRUE,save.matrices=TRUE)
#' # Run top model from paper for all strata using simulated annealing (commented out)
#' # mod_msld=crm(ps_dp,ddl,model.parameters=list(S=S.1,p=p.1,r=r.1,Psi=Psi.1),
#' #                     method="SANN",itnmax=6e6,hessian=TRUE)
#'}
NULL


