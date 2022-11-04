#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "head.h"
#include "head_multi.h"




// SCHEMA BBC EN ENERGIE INTERNE  JCP 2009

/////////////////////  MAIN  ////////////////////
int funcBBC_JCP2009_multi(int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int EOS,
                          int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                          double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                          double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
 
  int ind_air, ind_etain; 
  ind_air=0; ind_etain=1;

    // variables
  double t=0;
  double dt,dt_new;
  int n=0,err=0;
  double minDMsRC; 
  char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=10; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=4; 
  // l'idee  : appliquer les cond aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double dmm, dmp;
  double TAUtemp;
  double dX;
  int zone;
  double dUdx;
  int n_iter;
  double* pVAR=malloc(2*sizeof(double));

  double dx=(b-a)/nx; // pas d'espace


  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2

  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage
  // PW
  double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
  double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
  // solution 
  double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* U=malloc((N+1)*sizeof(double)); //  u (frontières)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* E=malloc((N)*sizeof(double)); //  energie totale (centrées aux mailles)

  // solution intermedaire
  double* TAU2=malloc((N)*sizeof(double)); //  tau au temsp n+1/2 (centrées aux mailles)
  double* U2=malloc((N+1)*sizeof(double)); //  u au temsp n+1/2 (frontières)
  double* U4=malloc((N+1)*sizeof(double)); //  u au temsp n+1/4 (frontières)
  double* EPS2=malloc((N)*sizeof(double)); //  epsilon au temsp n+1/2 (centrées aux mailles)


  double* RHO=malloc((N)*sizeof(double));
  double* P=malloc((N)*sizeof(double));
  double* Q=malloc((N)*sizeof(double));   // artificial viscosity
  double* C=malloc((N)*sizeof(double));
  double* DMsRC=malloc((N)*sizeof(double));  // dX/c   centres aux mailles
  double* DMc=malloc((N)*sizeof(double));  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles
  double* DMd=malloc((N+1)*sizeof(double));  // Delta m = (Dm + Dm)/2   dé-centres aux mailles
  
  double* P2=malloc((N)*sizeof(double));
  double* C2=malloc((N)*sizeof(double));

  // Cinétique
  double* TAUn=malloc((N)*sizeof(double));
  double* Un=malloc((N+1)*sizeof(double));
  
  double* T=malloc((N)*sizeof(double));  // température 

  double** matFM=alloctabd(N,3);          // fraction massique
  double* lambda_out_unuse=malloc(3*sizeof(double));

  double* G=malloc((N)*sizeof(double)); // fraction massique
  double* S=malloc((N)*sizeof(double)); // fraction massique
  
  int* tabNiter=malloc((N)*sizeof(int)); // nb iterration par mailles
  for(int i=0; i<N; ++i){tabNiter[i]=0;}

  double* tabEPSILON=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-1e-15;}

  double* test_consis_thermo=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){test_consis_thermo[i]=-99;}


  int* indFLUIDE=malloc((N)*sizeof(int)); // indice du fluide sur chaque maille
  

  // Maillage
  for(int i=0; i<=N; i++){
    X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
  }

  err=init_fluides(tst_multi, indFLUIDE, N, ind_air, ind_etain);
  if(err){printf("Erreur fluide\n"); return 1;}
  

  // Initialisation 
  double x;

  double* Wpw=malloc(10*sizeof(double));
  double* Wav=malloc(10*sizeof(double));
  for(int i=0; i<=N-1; i++){
    //printf("  --i= %d\n", i);
    x=(X[i]+X[i+1])/2.0;

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, X[i],X[i+1], Wav);
    if(err){return err;}
    
    RHO0c[i]=Wav[0];    
    TAU[i]=Wav[1];
    EPS[i]=Wav[4];
    
    P[i]=Wpw[2];
    T[i]=Wpw[3];  
    matFM[i][0]=Wpw[4];
    matFM[i][1]=Wpw[5];
    matFM[i][2]=Wpw[6];
  }

  // variables défini aux frontières des mailles
  for(int i=0; i<=N; i++){
    x=X[i];

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, x-dx/2, x+dx/2, Wav);
    if(err){return err;}
    
    RHO0d[i]=Wav[0];
    U[i]=Wav[2];
  }
  free(Wpw);
  free(Wav);


  //Init des tableaux  C
  for(int i=0; i<=N-1; i++){
    if(indFLUIDE[i]==ind_etain){
      err= G_S_cin_3_phase(epsilon, P[i], TAU[i], EPS[i], matFM[i], &C[i], &G[i], &S[i]);  if(err){return err;}
    }
    else if(indFLUIDE[i]==ind_air){
      EGS(tst_air,TAU[i], EPS[i], &P[i], &C[i]);
    }
    DMc[i]=dx*RHO0c[i];
    DMsRC[i]=DMc[i]*(TAU[i]/C[i]);
    //printf("i=%d P= %lf, C=%lf, U=%lf, EPS=%lf, XsC=%lf \n",i,P[i],C[i],U[i],REPS[i]/RHO0c[i], XsC[i]);           
  }


  for(int i=ideb; i<=ifin+1; i++){
    DMd[i]=(DMc[i-1]+DMc[i])/2;          
  }

  // energie totale
  for(int i=ideb; i<=ifin; i++){ 
    dmm=(DMc[i-1]+DMc[i])*(U[i]*U[i])/DMc[i];
    dmp=(DMc[i+1]+DMc[i])*(U[i+1]*U[i+1])/DMc[i];
    E[i]=EPS[i] + (1./8)*(dmm+dmp);
  }


    // minimum de DMsRC  (calcul de dt)
  minDMsRC=DMsRC[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(minDMsRC>DMsRC[j]){
      minDMsRC=DMsRC[j];
    }
  }

  //printf("minDMsRC= %.8lf\n",minDMsRC );
  dt = dt_first*CFL * minDMsRC; // pas de temps
  printf("nx= %d et dt= %g\n", nx, dt);
  
  // energie totale du système
  int nbiter=10*ceil(Tf/dt);
  //printf("nbiter :%d\n",nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double));  // masse totale
  double* VOL=malloc(nbiter*sizeof(double));   // volume total
  

  sum1=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*RHO0c[i]*E[i];
    sum3+=dx*RHO0d[i]*U[i];
    sum4+=(X[i+1]-X[i])/TAU[i];
    sum5+=dx*RHO0c[i]*TAU[i];
  }
  sum3+=dx*RHO0d[ifin+1]*U[ifin+1];
  Etot[0]=sum1;
  IMPUL[0]=sum3;
  Mtot[0]=sum4;
  VOL[0]=sum5;
  


// debut de la boucle ============================================================================= 
  while((t<Tf)&&(n<Nmax)){
    // cond aux limites
    err=condlim(ind_cond, ideb,ifin, nbghosts, TAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U, iu);
    if(err){return err;}

    // 1 - Computation of the aritificial viscosity and first prediction of the velocity at time n+1/4
    for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
      Q[i]=qvis(q,Cq,Cl,&U[i], TAU[i], C[i], &DMc[i]);
    }

    for(int i=ideb; i<=ifin+1; i++){
      U4[i]=U[i] - (dt/4)*( (P[i]+Q[i]) - (P[i-1]+Q[i-1]) )/DMd[i];           
    }
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U4, iu);
    if(err){return err;}

    // 2 - Prediction at time n+1/2
    for(int i=ideb; i<=ifin; i++){
      TAU2[i]=TAU[i] + 0.5*(dt/DMc[i])*(U4[i+1]-U4[i]);           
    }
    err=condlim(ind_cond, ideb,ifin, nbghosts, TAU2, 0);
    if(err){return err;}

    for(int i=ideb; i<=ifin; i++){
      EPS2[i]=EPS[i] -(P[i]+Q[i])*(TAU2[i]-TAU[i]);        
    }
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS2, 0);
    if(err){return err;}
    

    for(int i=0; i<=N-1; i++){
      dX=X[i+1]-X[i];
      dUdx=(U[i+1]-U[i])/dX;
      if(indFLUIDE[i]==ind_etain){
        err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU2[i], EPS2[i], &P2[i], &C2[i], &T[i], pVAR, matFM[i], lambda_out_unuse, 
                 Ntau, dX, dt/2, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAU[i], C[i],
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                 SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){return err;}
        if(C2[i]==0){C2[i]=1e-8;}
      }
      else if(indFLUIDE[i]==ind_air){
        //EGS(tst_air, TAU2[i], EPS2[i], &P2[i], &C2[i]);
        EGS_cin_phase(tst_air, TAU2[i], EPS2[i], &P2[i], &C2[i], 
                      cin_phase, dX, dt/2, Ttau, 
                      &U[i], P[i], TAU[i], C[i]);
      }
    }

    
    for(int i=ideb; i<=ifin+1; i++){ ///////////
      U2[i]=U[i] - (dt/2)*( (P2[i]+Q[i]) - (P2[i-1]+Q[i-1]) )/DMd[i];
    }
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U2, iu);
    if(err){return err;}

    // Cinétique
    for(int i=0; i<=N-1; i++){
      Un[i]= U[i];
      TAUn[i]=TAU[i];
    }
    Un[N]= U[N];

    // 3 - Finalization at time n+1
    for(int i=0; i<=N; i++){
      X[i]+= dt*U2[i];
      U[i]=2*U2[i]-U[i];
    }

    for(int i=ideb; i<=ifin; i++){
      TAUtemp=TAU[i] + (dt/DMc[i])*(U2[i+1] - U2[i]);   // TAU au temps n+1
      EPS[i]-= (P2[i]+Q[i])*(TAUtemp-TAU[i]);
      TAU[i]=TAUtemp;
    }

    err=condlim(ind_cond, ideb,ifin, nbghosts, TAU2, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS2, 0);
    if(err){return err;}

    for(int i=0; i<=N-1; i++){
      dX=X[i+1]-X[i];
      dUdx=(Un[i+1]-Un[i])/dX;
      if(indFLUIDE[i]==ind_etain){
        err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], matFM[i], 
                 Ntau, dX, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                 SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){return err;}
        if(C[i]==0){C[i]=1e-8;}
      }
      else if(indFLUIDE[i]==ind_air){
        //EGS(tst_air,TAU[i], EPS[i], &P[i], &C[i]);
        EGS_cin_phase(tst_air, TAU[i], EPS[i], &P[i], &C[i], 
                      cin_phase, dX, dt, Ttau, 
                      &Un[i], P[i], TAUn[i], C[i]);
      }
    }


    // energie totale
    for(int i=ideb; i<=ifin; i++){ 
      dmm=(DMc[i-1]+DMc[i])*(U[i]*U[i])/DMc[i];
      dmp=(DMc[i+1]+DMc[i])*(U[i+1]*U[i+1])/DMc[i];
      E[i]=EPS[i] + (1./8)*(dmm+dmp);
    }



    ///////////////
    // Energie du système
    sum1=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=DMc[i]*E[i];
      sum3+=dx*RHO0d[i]*U[i];
      sum4+=(X[i+1]-X[i])/TAU[i];
      sum5+=dx*RHO0c[i]*TAU[i];
    }
    sum3+=dx*RHO0d[ifin+1]*U[ifin+1];
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    

    // Mise a jour de XsC
    for(int i=ideb; i<=ifin; i++){
      DMsRC[i]=DMc[i]*TAU[i]/C[i];
    }
     
    t=t+dt;  n=n+1;
    if(aff==1){ printf("n= %d, t= %g, dt= %g\n",n,t,dt); }

    // Calcul du pas de temps
    minDMsRC=DMsRC[ideb];
    for(int j=ideb+1; j<=ifin; j++){  if (minDMsRC>DMsRC[j]){   minDMsRC=DMsRC[j];  }  }
    dt_new = CFL * minDMsRC; // pas de temps
    dt=pastemps(Tf,t,dt_new,dt);

  }
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);

  // energie totale
  for(int i=ideb; i<=ifin; i++){ 
    Xc[i]=(X[i+1]+X[i])/2;
    RHO[i]=1./TAU[i];
    
    if(indFLUIDE[i]==ind_etain){ G_S_cin_3_phase(epsilon, P[i], TAU[i], EPS[i], matFM[i], &C[i], &G[i], &S[i]); }
    else if(indFLUIDE[i]==ind_air){ G_EGS(tst_air, TAU[i], EPS[i], &G[i]); }
  }

  // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, RHO, TAU, U, E, EPS, P, T, C, G, S, matFM);
  if(err){return err;}

  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}
  
  FILE* Res2;
  if((Res2 = fopen("nb_iter.txt", "w+")) != NULL){
    for(int i=ideb; i<=ifin; i++){  fprintf(Res2,"%d ",tabNiter[i]);  }
    fprintf(Res2, "0\n");
    for(int i=ideb; i<=ifin; i++){  fprintf(Res2,"%g ",tabEPSILON[i]);  }
    fprintf(Res2, "-1e-10\n");
    fclose(Res2);
  } 
  else{printf("erreur impression\n"); return 1;}


  free(X); free(Xc); free(indFLUIDE);
  free(RHO0d); free(RHO0c);
  free(DMc); free(DMd);
  
  free(TAU); free(U); free(EPS);

  free(P); free(Q); free(C); free(DMsRC);
  free(E); free(RHO);

  free(P2); free(C2);
  free(TAU2); free(U2); free(U4); free(EPS2);
  
  freetab(matFM); free(T); free(lambda_out_unuse);
  //  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
  
}
