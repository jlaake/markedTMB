// TMB Version: Mixed-effect Multi-State Jolly-Seber model
// Jeff Laake; 13 Oct 2021

#include <TMB.hpp>                              // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);                            // number of capture histories
  DATA_INTEGER(m);                            // number of capture occasions
  DATA_INTEGER(nS);                           // number of states excluding death state
  DATA_IMATRIX(ch);                           // capture history matrix; uses numeric values for states
  DATA_IVECTOR(frst);                         // occasion first seen for each history
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
  
  DATA_INTEGER(nrowpent);                      // number of rows in the simplified design matrix for pent
  DATA_MATRIX(pentdm);                         // design matrix for pent
  DATA_VECTOR(pentfix);                        // pent fixed values
  DATA_IVECTOR(pentindex);                     // pent indices
  DATA_INTEGER(pent_nre);                      // number of random effects for pent
  DATA_INTEGER(pent_krand);                    // number of columns in pent random effect DM
  DATA_MATRIX(pent_randDM);                    // pent random effect DM
  DATA_IVECTOR(pent_randDM_i);                 // pent random DM indices
  DATA_IMATRIX(pent_randIndex);                // pent random effect indices for DM; index into phi_u
  DATA_IVECTOR(pent_randIndex_i);              // pent random effect index indices
  DATA_IVECTOR(pent_counts);                   // count of pent random effect indices by id
  DATA_IMATRIX(pent_idIndex);                  // pent random effect indices by id; index into u_phi to construct phi_u
  DATA_IVECTOR(pent_idIndex_i);                // pent random effect id index indices
  
  DATA_INTEGER(nrowpi);                      // number of rows in the simplified design matrix for pi
  DATA_MATRIX(pidm);                         // design matrix for pi
  DATA_VECTOR(pifix);                        // pi fixed values
  DATA_IVECTOR(piindex);                     // pi indices
  DATA_INTEGER(pi_nre);                      // number of random effects for pi
  DATA_INTEGER(pi_krand);                    // number of columns in pi random effect DM
  DATA_MATRIX(pi_randDM);                    // pi random effect DM
  DATA_IVECTOR(pi_randDM_i);                 // pi random DM indices
  DATA_IMATRIX(pi_randIndex);                // pi random effect indices for DM; index into phi_u
  DATA_IVECTOR(pi_randIndex_i);              // pi random effect index indices
  DATA_IVECTOR(pi_counts);                   // count of pi random effect indices by id
  DATA_IMATRIX(pi_idIndex);                  // pi random effect indices by id; index into u_phi to construct phi_u
  DATA_IVECTOR(pi_idIndex_i);                // pi random effect id index indices
  
  DATA_INTEGER(getreals);                     // if 1, report reals and std errors
  
  PARAMETER_VECTOR(phibeta);                  // parameter vector for Phi
  PARAMETER_VECTOR(pbeta);                    // parameter vector for p
  PARAMETER_VECTOR(psibeta);                  // parameter vector for Psi
  PARAMETER_VECTOR(pentbeta);                 // parameter vector for pent
  PARAMETER_VECTOR(pibeta);                   // parameter vector for pi
  PARAMETER_VECTOR(log_sigma_phi);
  PARAMETER_VECTOR(log_sigma_p);
  PARAMETER_VECTOR(log_sigma_psi);
  PARAMETER_VECTOR(log_sigma_pent);
  PARAMETER_VECTOR(log_sigma_pi);
  PARAMETER_VECTOR(u_phi);
  PARAMETER_VECTOR(u_p);
  PARAMETER_VECTOR(u_psi);
  PARAMETER_VECTOR(u_pent);
  PARAMETER_VECTOR(u_pi);
  
  Type g=0;
  
  int nrows;                           // number of entries in design matrix m-1 values for each of nS states
  nrows=nS*(m-1);
  int nT;                              // number of transitions excluding death
  nT=nS*nS*(m-1);
  
  int i,j,k,bindex,bindex2,bindex3,bindex4,k2,idx,i2;  // indices and counters
  int L;
  vector<Type> uniquephi(nrowphi);     // all unique phi values
  vector<Type> phi(nrows);             // temp vector for Phis for an individual
  vector<Type> uniquep(nrowp);         // all unique p values
  vector<Type> p(nrows);               // temp vector for ps for an individual
  vector<Type> uniquepsi(nrowpsi);     // temp vector for psis
  Type psisum;                         // sum of psi for each state to normalize with
  vector<Type> uniquepent(nrowpent);   // temp vector for pent
  Type pentsum;                        // sum of pent across time
  vector<Type> uniquepi(nrowpi);       // temp vector for pi
  Type pisum;                          // sum of pent across states within a time
  
  array<Type> psi(m-1,nS,nS);         // matrix for psis for each occasion
  vector<Type> pent(m-1);             // vector for pent by occasion
  vector<Type> cumsumpent(m-2);       // vector for cumulative sums of pent by occasion
  array<Type> pi(m-1,nS-1);           // matrix for pi by occasion / stratum
  array<Type> gamma(m-1,nS+1,nS+1);   // transition probability matrices for individual i
  array<Type> dmat(m-1,nS+1,nS+1);    // observation probability matrices for individual i
  array<double> allgamma(n,m-1,nS+1,nS+1);   // transition probability matrices for all individuals
  array<double> alldmat(n,m-1,nS+1,nS+1);      // observation probability matrices  for all individuals
  Type u;                             // sum of state probabilities
  vector<Type> pS(nS+1);              // update vector for prob of being in state j=1,nS + death
  vector<Type> S(nS+1);               // prob of being in state j=1,nS + death for each occasion
  vector<Type> v(nS+1);               // temporary update vector
  vector<Type> vec;                   // temporary vector
  Type Lglki=0;                       // log-likelihood accumulator
  Type mu;
  Type p0;
  int nphicounts=n;                   // number of counts for phi random effects by id
  if(phi_nre==0)nphicounts=0;
  int npcounts=n;                     // number of counts for p random effects by id
  if(p_nre==0)npcounts=0;
  int npsicounts=n;                   // number of counts for psi random effects by id
  if(psi_nre==0)npsicounts=0;
  int npentcounts=n;                   // number of counts for pent random effects by id
  if(pent_nre==0)npentcounts=0;
  int npicounts=n;                   // number of counts for pent random effects by id
  if(pi_nre==0)npicounts=0;
  
  if(phi_krand>0)	                                        // likelihood contribution for n(0,1) re for phi
    for (int i=0;i<=phi_nre-1;i++)
      g-= dnorm(u_phi(i),Type(0),Type(1),true);
  
  if(p_krand>0)	                                         // likelihood contribution for n(0,1) re for p
    for (int i=0;i<=p_nre-1;i++)
      g-= dnorm(u_p(i),Type(0),Type(1),true);
  
  if(psi_krand>0)	                                        // likelihood contribution for n(0,1) re for psi
    for (int i=0;i<=psi_nre-1;i++)
      g-= dnorm(u_psi(i),Type(0),Type(1),true);
  
  if(pent_krand>0)	                                      // likelihood contribution for n(0,1) re for pent
    for (int i=0;i<=pent_nre-1;i++)
      g-= dnorm(u_pent(i),Type(0),Type(1),true);
  
  if(pi_krand>0)	                                      // likelihood contribution for n(0,1) re for pi
    for (int i=0;i<=pi_nre-1;i++)
      g-= dnorm(u_pi(i),Type(0),Type(1),true);
  
  uniquephi=phidm*phibeta;                              // compute unique parameter sets on link scale
  uniquep=pdm*pbeta;
  uniquepsi=psidm*psibeta;
  uniquepent=pentdm*pentbeta;
  uniquepi=pidm*pibeta;
  alldmat.setZero();
  allgamma.setZero();
  
  // loop over capture histories
  for(i=1;i<=n;i++)                                    
  {
    vector<Type> p_u(p_idIndex.cols());        // define random effects vector for p, Phi and psi used
    vector<Type> phi_u(phi_idIndex.cols());    // just for this capture history copied from full vectors
    vector<Type> psi_u(psi_idIndex.cols());
    vector<Type> pent_u(pent_idIndex.cols());
    vector<Type> pi_u(pi_idIndex.cols());
    p_u.setZero();
    phi_u.setZero();
    if(nphicounts >0)                          // if any random effects for phi, copy values from u_phi to phi_u
    {
      if(phi_counts(i-1)==0)
        phi_u(0)=0;
      else
        for(j=0;j<=phi_counts(i-1)-1;j++)
          phi_u(j)=u_phi(phi_idIndex(phi_idIndex_i(i-1)-1,j)-1);
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
    
    if(npentcounts >0)                           // if any random effects for pent, copy values from u_pent to pent_u
    {
      if(pent_counts(i-1)==0)
        pent_u(0)=0;
      else
        for(j=0;j<=pent_counts(i-1)-1;j++)
          pent_u(j)=u_pent(pent_idIndex(pent_idIndex_i(i-1)-1,j)-1);
    }
    
    if(npicounts >0)                           // if any random effects for pi, copy values from u_pi to pi_u
    {
      if(pi_counts(i-1)==0)
        pi_u(0)=0;
      else
        for(j=0;j<=pi_counts(i-1)-1;j++)
          pi_u(j)=u_pi(pi_idIndex(pi_idIndex_i(i-1)-1,j)-1);
    }
    
    bindex=(i-1)*nrows;                               // initialize indices into index values for the ith history
    bindex2=(i-1)*nT;
    bindex3=(i-1)*(m-1);
    bindex4=(i-1)*(m-1)*(nS-1);
    pentsum=0;
    
    // loop over occasions
    for (j=1;j<=m-1;j++)
    {
      // compute pent values over time 
      i2=bindex3+j;
      idx=pentindex(i2-1)-1;
      if(pentfix(idx)< -0.5)
      {
        mu=0;
        if(npentcounts>0)
          if(pent_counts(i-1) > 0)	                        // random portion of mean if any
          {
            for(L=1;L<=pent_krand;L++)
              if(pent_randIndex(pent_randIndex_i(i2-1)-1,L-1)>0)
                mu+=pent_randDM(pent_randDM_i(i2-1)-1,L-1)*pent_u(pent_randIndex(pent_randIndex_i(i2-1)-1,L-1)-1)*exp(log_sigma_pent(L-1));
          }
          if((uniquepent(idx)+mu)< -700)
            pent(j-1)=exp(-700);
          else  
            if((uniquepent(idx)+mu)> 700)
              pent(j-1)=exp(700);
            else
              pent(j-1)=exp(uniquepent(idx)+mu);
      }
      else
        pent(j-1)=pentfix(idx);
      pentsum+=pent(j-1);
      
      // loop over strata
      pisum=0;
      for (k=1;k<=nS;k++)
      {
        // compute p values by state for the jth occasion
        i2=bindex+(j-1)*nS+k;
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
            if((uniquep(idx)+mu)< -25)
              p((j-1)*nS+k-1)=1/(1+exp(25));
            else  
              if((uniquep(idx)+mu)> 25)
                p((j-1)*nS+k-1)=1/(1+exp(-25));
              else
                p((j-1)*nS+k-1)=1/(1+exp(-(uniquep(idx)+mu)));
        }
        else
          p((j-1)*nS+k-1)=pfix(idx);
        
        // compute phi values by state for the jth occasion
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
        
        // compute pi values for each state (excluding N) for the jth occasion
        if(k>1)
        {
          i2=bindex4+(j-1)*(nS-1)+k-1;
          idx=piindex(i2-1)-1;
          if(pifix(idx)< -0.5)
          {
            mu=0;
            if(npicounts>0)
              if(pi_counts(i-1) > 0)	                        // random portion of mean if any
              {
                for(L=1;L<=pi_krand;L++)
                  if(pi_randIndex(pi_randIndex_i(i2-1)-1,L-1)>0)
                    mu+=pi_randDM(pi_randDM_i(i2-1)-1,L-1)*pi_u(pi_randIndex(pi_randIndex_i(i2-1)-1,L-1)-1)*exp(log_sigma_pi(L-1));
              }
              if((uniquepi(idx)+mu)< -700)
                pi(j-1,k-2)=exp(-700);
              else  
                if((uniquepi(idx)+mu)> 700)
                  pi(j-1,k-2)=exp(700);
                else
                  pi(j-1,k-2)=exp(uniquepi(idx)+mu);
          }
          else
            pi(j-1,k-2)=pifix(idx);
          pisum+=pi(j-1,k-2);
        }
        
        //  compute psi values for the individual
        psisum=0;
        for(k2=1;k2<=nS;k2++)
        {
          i2=bindex2+(k-1)*nS+k2;
          idx=psiindex(i2-1)-1;
          if(psifix(idx)< -0.5)
          {
            mu=0;
            if(npsicounts>0)
              if(psi_counts(i-1) > 0)
              {
                for(L=1;L<=psi_krand;L++)
                  if(psi_randIndex(psi_randIndex_i(i2-1)-1,L-1)>0)
                    mu+=psi_randDM(psi_randDM_i(i2-1)-1,L-1)*psi_u(psi_randIndex(psi_randIndex_i(i2-1)-1,L-1)-1)*exp(log_sigma_psi(L-1));
              }
              if((uniquepsi(idx)+mu)>700)
                psi(j-1,k-1,k2-1)=exp(700);
              else
                if((uniquepsi(idx)+mu)< -700)
                  psi(j-1,k-1,k2-1)=exp(-700);
                else
                  psi(j-1,k-1,k2-1)=exp(uniquepsi(idx)+mu);
          }
          else
            psi(j-1,k-1,k2-1)=psifix(idx);
          psisum+=psi(j-1,k-1,k2-1);
        }
        for(k2=1;k2<=nS;k2++)
          psi(j-1,k-1,k2-1)=psi(j-1,k-1,k2-1)/psisum;
        
      } // end of loop over states
      bindex2=bindex2+nS*nS;
      
      
      // normalize pi over strata
      for (k=2;k<=nS;k++)
        pi(j-1,k-2)=pi(j-1,k-2)/pisum;
    } // end of loop over occasions
    
    // normalize pent to sum to 1 by looping over occasions
    for (j=1;j<=m-1;j++)
      pent(j-1)=pent(j-1)/pentsum;
    // now adjust pent such that all are moved out of N; first compute cummulative sums
    for (j=1;j<=m-2;j++)
    {
      if(j>1)
        cumsumpent(j-1)=cumsumpent(j-2)+pent(j-1);
      else
        cumsumpent(0)=pent(0);
    }  
    for (j=2;j<=m-1;j++)
      pent(j-1)=pent(j-1)/(1-cumsumpent(j-2));
    
    if(getreals>0)                                               // if requested report phi,p,psi,pent and f0 values
    {
      ADREPORT(phi);
      ADREPORT(p);
      ADREPORT(psi);
      ADREPORT(pent);
    }
    
    //  compute transition matrices for each occasion
    gamma.setZero();                        // initialize all transitions to zero
    bindex=1;
    for(j=1;j<=m-1;j++)                        // loop over intervals
    {
      for(k=1;k<=nS;k++)                     // loop over states creating p and gamma values
      {
        for(k2=1;k2<=nS;k2++)
        {
          if(k==1)
          {
            if(k2>1)
              gamma(j-1,0,k2-1)=pent(j-1)*pi(j-1,k2-2);  // entry to population from state N
            else
              gamma(j-1,0,0)=1-pent(j-1); //stay in state N (not enter at this occasion)
          }
          else
            gamma(j-1,k-1,k2-1)=psi(j-1,k-1,k2-1)*phi(bindex-1);    // adjust psi for survival
          allgamma(i-1,j-1,k-1,k2-1)=asDouble(gamma(j-1,k-1,k2-1));
        }
        gamma(j-1,k-1,nS)=1-phi(bindex-1);              // add death state value for each state
        allgamma(i-1,j-1,k-1,nS)=asDouble(gamma(j-1,k-1,nS));
        bindex++;
      }
      gamma(j-1,nS,nS)=1;                             // death is an absorbing state
      allgamma(i-1,j-1,nS,nS)=1;
    }
    
    //  compute state dependent observation matrices for each occasion
    dmat.setZero();
    bindex=1;
    for(j=1;j<=m-1;j++)
    {
      for(k=1;k<=nS;k++)
      {
        dmat(j-1,k,k-1)=p(bindex-1);
        alldmat(i-1,j-1,k,k-1)=asDouble(dmat(j-1,k,k-1));
        dmat(j-1,0,k-1)=1-dmat(j-1,k,k-1);
        alldmat(i-1,j-1,0,k-1)=asDouble(dmat(j-1,0,k-1));
        bindex++;
      }
      dmat(j-1,0,nS)=1;
      alldmat(i-1,j-1,0,nS)=asDouble(dmat(j-1,0,nS));
    }
    //  HMM algorithm
    pS.setZero();                                      // initialize values to 0
    Lglki=0;
    S.setZero();
    S(ch(i-1,frst(i-1)-1)-1)=1;                        // set state prob to 1 for observed state at first observation
    for(j=frst(i-1)+1;j<=m;j++)                      // loop over possible occasions from first(i)+1 to m
    {
      for(k=1;k<=nS+1;k++)
      {
        pS(k-1)=0;
        for(k2=1;k2<=nS+1;k2++)
        {
          pS(k-1)+= S(k2-1)*gamma(j-2,k2-1,k-1);
        }
      }
      for(k=1;k<=nS+1;k++) 
      {
        v(k-1)=pS(k-1)*dmat(j-2,ch(i-1,j-1),k-1); // v is temp state vector alpha in Z&M
      }
      if(v.sum()==0) {
        Rcout << "\n Check Psi or p values set to 0";
        Rcout << "\n i = " << i << " ch = " << ch(i-1);
      }
      u=v.sum();                                       // sum across states
      S=v/u;                                          // update S;S is normalized alpha(j) (phi) in Z&M
      Lglki+=log(u);    	                            // accumulate log-likelihood value
      
    }
    if(freq(i-1)==0)
      p0=exp(Lglki);
    else
    {
      g-=freq(i-1)*Lglki;
      g+=freq(i-1)*log(1-p0);
    }
    
  }  // end of loop over individuals
  
  REPORT(alldmat);
  REPORT(allgamma);
  return g;
}
