// TMB Version: Mixed-effect Multi-State Cormack-Jolly-Seber + recovery model with unobservable states
// Jeff Laake; 9 Jan 2020 - added dmat/gamma report

#include <TMB.hpp>                              // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);                            // number of capture histories
  DATA_INTEGER(m);                            // number of capture occasions
  DATA_INTEGER(nS);                           // number of states excluding death states (newly dead and dead)
  DATA_IMATRIX(ch);                           // capture history matrix; uses numeric values for states
  DATA_IVECTOR(frst);                         // occasion first seen for each history
  DATA_IVECTOR(frstrec);                      // first record in ddl of S,r,p for each history
  DATA_IVECTOR(frstrecpsi);                   // first record in ddl of Psi foreach history
  DATA_VECTOR(freq);                          // frequency of each capture history
  DATA_MATRIX(tint);                          // time interval between occasions for each history-interval
  
  DATA_INTEGER(nrowphi);                      // number of rows in the simplified design matrix for Phi - survival
  // last column in each design matrix is either -1 (estimated) or a fixed value
  DATA_MATRIX(phidm);                         // design matrix for Phi
  DATA_VECTOR(phifix);                        // phi fixed values
  DATA_IVECTOR(phiindex);                     // phi indices
  DATA_INTEGER(phi_nre);                      // number of random effects for phi
  DATA_INTEGER(phi_krand);                    // number of columns in phi random effect DM
  DATA_MATRIX(phi_randDM);                    // phi random effect DM
  DATA_IVECTOR(phi_randDM_i);                 // phi random DM indices
  DATA_IMATRIX(phi_randIndex);                // phi random effect indices for DM
  DATA_IVECTOR(phi_randIndex_i);              // phi random effect index indices
  DATA_IVECTOR(phi_counts);                   // count of phi random effect indices by id
  DATA_IMATRIX(phi_idIndex);                  // phi random effect indices by id
  DATA_IVECTOR(phi_idIndex_i);                // phi random effect id index indices
  
  DATA_INTEGER(nrowr);                        // number of rows in the simplified design matrix for r - recovery
  // last column in each design matrix is either -1 (estimated) or a fixed value
  DATA_MATRIX(rdm);                           // design matrix for r
  DATA_VECTOR(rfix);                          // r fixed values
  DATA_IVECTOR(rindex);                       // r indices
  DATA_INTEGER(r_nre);                        // number of random effects for r
  DATA_INTEGER(r_krand);                      // number of columns in r random effect DM
  DATA_MATRIX(r_randDM);                      // r random effect DM
  DATA_IVECTOR(r_randDM_i);                   // r random DM indices
  DATA_IMATRIX(r_randIndex);                  // r random effect indices for DM
  DATA_IVECTOR(r_randIndex_i);                // r random effect index indices
  DATA_IVECTOR(r_counts);                     // count of r random effect indices by id
  DATA_IMATRIX(r_idIndex);                    // r random effect indices by id
  DATA_IVECTOR(r_idIndex_i);                  // r random effect id index indices
  
  DATA_INTEGER(nrowp);                        // number of rows in the simplified design matrix for p
  DATA_MATRIX(pdm);                           // design matrix for p
  DATA_VECTOR(pfix);                          // p fixed values
  DATA_IVECTOR(pindex);                       // p indices
  DATA_INTEGER(p_nre);                        // number of random effects for p
  DATA_INTEGER(p_krand);                      // number of columns in p random effect DM
  DATA_MATRIX(p_randDM);                      // p random effect DM
  DATA_IVECTOR(p_randDM_i);                   // p random DM indices
  DATA_IMATRIX(p_randIndex);                  // p random effect indices for DM; index into p_u
  DATA_IVECTOR(p_randIndex_i);                // p random effect index indices
  DATA_IVECTOR(p_counts);                     // count of p random effect indices by id
  DATA_IMATRIX(p_idIndex);                    // p random effect indices by id; index into u_phi to construct phi_u
  DATA_IVECTOR(p_idIndex_i);                  // p random effect id index indices
  
  DATA_INTEGER(nrowpsi);                      // number of rows in the simplified design matrix for psi
  DATA_MATRIX(psidm);                         // design matrix for psi
  DATA_VECTOR(psifix);                        // psi fixed values
  DATA_IVECTOR(psiindex);                     // psi indices
  DATA_INTEGER(psi_nre);                      // number of random effects for psi
  DATA_INTEGER(psi_krand);                    // number of columns in psi random effect DM
  DATA_MATRIX(psi_randDM);                    // psi random effect DM
  DATA_IVECTOR(psi_randDM_i);                 // psi random DM indices
  DATA_IMATRIX(psi_randIndex);                // psi random effect indices for DM; index into phi_u
  DATA_IVECTOR(psi_randIndex_i);              // psi random effect index indices
  DATA_IVECTOR(psi_counts);                   // count of psi random effect indices by id
  DATA_IMATRIX(psi_idIndex);                  // psi random effect indices by id; index into u_phi to construct phi_u
  DATA_IVECTOR(psi_idIndex_i);                // psi random effect id index indices
  
  DATA_INTEGER(getreals);                     // if 1, report reals and std errors
  
  PARAMETER_VECTOR(phibeta);                  // parameter vector for Phi
  PARAMETER_VECTOR(rbeta);                    // parameter vector for r
  PARAMETER_VECTOR(pbeta);                    // parameter vector for p
  PARAMETER_VECTOR(psibeta);                  // parameter vector for Psi
  PARAMETER_VECTOR(log_sigma_phi);
  PARAMETER_VECTOR(log_sigma_r);
  PARAMETER_VECTOR(log_sigma_p);
  PARAMETER_VECTOR(log_sigma_psi);
  PARAMETER_VECTOR(u_phi);
  PARAMETER_VECTOR(u_r);
  PARAMETER_VECTOR(u_p);
  PARAMETER_VECTOR(u_psi);
  
  Type g=0; 
  
  int nrows;                           // number of entries in design matrix m-1 values for each of nS states
  nrows=nS*(m-1);                
 
  int i,j,k,bindex,bindex2,k2,idx,i2;  // indices and counters
  int L; 
  vector<Type> uniquephi(nrowphi);     // all unique phi values    
  vector<Type> phi(nrows);             // temp vector for Phis for an individual
  vector<Type> uniquer(nrowr);         // all unique r values    
  vector<Type> r(nrows);               // temp vector for rs for an individual
  vector<Type> uniquep(nrowp);         // all unique p values    
  vector<Type> p(nrows);               // temp vector for ps for an individual
  vector<Type> uniquepsi(nrowpsi);     // temp vector for psis 
  Type psisum;                         // sum of psi for each state to normalize with
  
  phi.setZero();
  p.setZero();
  r.setZero();
  
  array<Type> psi(m-1,nS,nS);             // matrix for psis for each occasion 
  array<Type> gamma(m-1,2*nS+1,2*nS+1);   // transition probability matrices for individual i
  array<Type> dmat(m-1,nS+2,2*nS+1);      // observation probability matrices for individual i
  array<double> allgamma(n,m-1,2*nS+1,2*nS+1);   // transition probability matrices for all individuals
  array<double> alldmat(n,m-1,nS+2,2*nS+1);      // observation probability matrices  for all individuals
  Type u;                                 // sum of state probabilities
  vector<Type> pS(2*nS+1);              // update vector for prob of being in state j=1,2*nS + 1       
  vector<Type> S(2*nS+1);               // prob of being in state j=1,2*nS + 1 for each occasion
  vector<Type> v(2*nS+1);               // temporary update vector
  vector<Type> vec;                   // temporary vector
  Type Lglki=0;                       // log-likelihood accumulator
  Type mu;
  int nphicounts=n;                   // number of counts for phi random effects by id
  if(phi_nre==0)nphicounts=0;
  int nrcounts=n;                      // number of counts for r random effects by id
  if(r_nre==0)nrcounts=0;
  int npcounts=n;                     // number of counts for p random effects by id
  if(p_nre==0)npcounts=0;
  int npsicounts=n;                   // number of counts for psi random effects by id
  if(psi_nre==0)npsicounts=0;
  
  if(phi_krand>0)	                                     // likelihood contribution for n(0,1) re for phi
    for (int i=0;i<=phi_nre-1;i++)	   	     
      g-= dnorm(u_phi(i),Type(0),Type(1),true);
  
  if(r_krand>0)	                                        // likelihood contribution for n(0,1) re for r
    for (int i=0;i<=r_nre-1;i++)	   	     
      g-= dnorm(u_r(i),Type(0),Type(1),true);
  
  if(p_krand>0)	                                        // likelihood contribution for n(0,1) re for p
    for (int i=0;i<=p_nre-1;i++)
      g-= dnorm(u_p(i),Type(0),Type(1),true);
  
  if(psi_krand>0)	                                        // likelihood contribution for n(0,1) re for psi
    for (int i=0;i<=psi_nre-1;i++)
      g-= dnorm(u_psi(i),Type(0),Type(1),true);
  
  uniquephi=phidm*phibeta;                              // compute unique parameter sets on link scale
  uniquer=rdm*rbeta;    
  uniquep=pdm*pbeta;
  uniquepsi=psidm*psibeta;
  alldmat.setZero();
  allgamma.setZero();
  
  for(i=1;i<=2;i++)                             // loop over capture histories - one per capture history
  {
    vector<Type> p_u(p_idIndex.cols());        // define random effects vector for p, Phi,r and psi used                  
    vector<Type> phi_u(phi_idIndex.cols());    // just for this capture history copied from full vectors
    vector<Type> r_u(r_idIndex.cols());
    vector<Type> psi_u(psi_idIndex.cols());
    p_u.setZero();	   
    phi_u.setZero();	   
    r_u.setZero();	   
    if(nphicounts >0)                          // if any random effects for phi, copy values from u_phi to phi_u
    {
      if(phi_counts(i-1)==0)
        phi_u(0)=0;
      else
        for(j=0;j<=phi_counts(i-1)-1;j++)
          phi_u(j)=u_phi(phi_idIndex(phi_idIndex_i(i-1)-1,j)-1);
    } 
    
    if(nrcounts >0)                          // if any random effects for r, copy values from u_r to r_u
    {
      if(r_counts(i-1)==0)
        r_u(0)=0;
      else
        for(j=0;j<=r_counts(i-1)-1;j++)
          r_u(j)=u_r(r_idIndex(r_idIndex_i(i-1)-1,j)-1);
    } 
    
    if(npcounts >0)                           // if any random effects for p, copy values from u_p to p_u
    {
      if(p_counts(i-1)==0)
        p_u(0)=0;
      else
        for(j=0;j<=p_counts(i-1)-1;j++)
          p_u(j)=u_p(p_idIndex(p_idIndex_i(i-1)-1,j)-1);
    } 
    
    if(npsicounts >0)                           // if any random effects for psi, copy values from u_psi to psi_u
    {
      if(psi_counts(i-1)==0)
        psi_u(0)=0;
      else
        for(j=0;j<=psi_counts(i-1)-1;j++)
          psi_u(j)=u_psi(psi_idIndex(psi_idIndex_i(i-1)-1,j)-1);
    } 
    //  compute phi and p values for the individual   
    bindex=frstrec(i-1);                               // initialize indices into index values for the ith history
    bindex2=frstrecpsi(i-1);
    Rcout << "\n" << frstrec(i-1);
    Rcout << "\n" << frstrecpsi(i-1);
    for (j=frst(i-1);j<=m-1;j++)
    {                       
      for (k=1;k<=nS;k++)                        
      {   
        i2=bindex+(j-frst(i-1))*nS+k;
        idx=pindex(i2-1)-1;
        if(pfix(idx)< -0.5)
        {
          mu=0;
          if(npcounts>0)
            if(p_counts(i-1) > 0)	                        // random portion of mean if any
            {
              for(L=1;L<=p_krand;L++)
                if(p_randIndex(p_randIndex_i(i2-1)-1,L-1)>0)
                  mu+=p_randDM(p_randDM_i(i2-1)-1,L-1)*p_u(p_randIndex(p_randIndex_i(i2-1)-1,L-1)-1)*exp(log_sigma_p(L-1));
            }	           
            p((j-1)*nS+k-1)=1/(1+exp(-(uniquep(idx)+mu)));
        }
        else
          p((j-1)*nS+k-1)=pfix(idx);
        
        idx=phiindex(i2-1)-1;
        if(phifix(idx)< -0.5)
        {
          mu=0;
          if(nphicounts>0)
            if(phi_counts(i-1) > 0)	                        // random portion of mean if any
            {
              for(L=1;L<=phi_krand;L++)
                if(phi_randIndex(phi_randIndex_i(i2-1)-1,L-1)>0)
                  mu+=phi_randDM(phi_randDM_i(i2-1)-1,L-1)*phi_u(phi_randIndex(phi_randIndex_i(i2-1)-1,L-1)-1)*exp(log_sigma_phi(L-1));
            }	
            phi((j-1)*nS+k-1)=pow(1/(1+exp(-(uniquephi(idx)+mu))),tint(i-1,j-1)); 
        }
        else
          phi((j-1)*nS+k-1)=phifix(idx);     
        
       }
    }
    if(i<5)
    {
      Rcout << "\n\ni =" <<i;
      Rcout << "\n phi =" << phi ;
      Rcout << "\n p =" << p ;
    }    
  }
  g=0;
  return g;
}
