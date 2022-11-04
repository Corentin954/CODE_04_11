#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "head.h"
#include "head_multi.h"




// SCHEMA : One-step 2nd Order Staggered-grid Scheme

/////////////////////  MAIN  ////////////////////
int funcKOVAo2_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS, double zeta,
                          int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                          double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                          double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
  
  printf(" zeta=%g\n",zeta );

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
  int zone;
  double dUdx;
  int n_iter;
  double minDXsC;
  double* pVAR=malloc(2*sizeof(double));
  double Ctriple=1e-12;


  // Calcul du flux
  double dXc,rc2,dXl,dXr,DxUl,DxUr,DxxU,DxU,rc2l,rc2r,DxRC2,DxTAU0;
  double DtU,DttU;
  double dX,DtP,dXcl,dXcr,DxPl,DxPr,DxxP,DxP,tau02,tau03,DxRHO0,DttP;


  double dx=(b-a)/nx; // pas d'espace

  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2

  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage

  double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
  double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
  double* RHO=malloc((N)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
  // solution 
  double* TAU=malloc((N)*sizeof(double)); //  tau   (centrées aux mailles)
  double* U=malloc((N+1)*sizeof(double)); //  u  (frontières)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)
  double* E=malloc((N)*sizeof(double)); //    enegie totale (frontieres)
  double* EKIN=malloc((N+1)*sizeof(double)); //  energie cinetique (frontières)

  double* Ustar=malloc((N+1)*sizeof(double));
  double* Pstar=malloc((N)*sizeof(double));

  double* P=malloc((N)*sizeof(double));
  double* C=malloc((N)*sizeof(double));
  double* T=malloc((N)*sizeof(double)); // Température 
  double* G=malloc((N)*sizeof(double)); // derive fondamentale
  double* S=malloc((N)*sizeof(double)); // entropie
  
  double* DK=malloc((N+1)*sizeof(double)); // kin fix

  // Cinétique
  double* TAUn=malloc((N)*sizeof(double));
  double* Un=malloc((N+1)*sizeof(double));
  
  double** matFM=alloctabd(N,3);          // fraction massique
  double* lambda_out_unuse=malloc(3*sizeof(double));

  int* tabNiter=malloc((N)*sizeof(int)); // nb iterration par mailles
  for(int i=0; i<N; ++i){tabNiter[i]=0;}
  double* tabEPSILON=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-1e-15;}
  double* test_consis_thermo=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){test_consis_thermo[i]=-99;}
  int* indFLUIDE=malloc((N)*sizeof(int)); // indice du fluide sur chaque maille
  
  // Tableau pour imprimer les derivées de P* et U*
  double* LdtU=malloc((N+1)*sizeof(double));
  double* LdttU=malloc((N+1)*sizeof(double));
  double* LdtP=malloc((N)*sizeof(double));
  double* LdttP=malloc((N)*sizeof(double));

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
    Xc[i]=x;

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, X[i],X[i+1], Wav);
    if(err){return err;}

    TAU[i]=Wav[1];
    EPS[i]=Wav[4];
    RHO0c[i]=Wav[0];    

    if(indFLUIDE[i]==ind_etain){
      err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], matFM[i], 
                     Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
                     nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                     SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){printf(" initialisation \n"); return err;}
      
      G[i]=pVAR[0];
      S[i]=pVAR[1];
      //err= G_S_cin_3_phase(epsilon, P[i], TAU[i], EPS[i], matFM[i], &C[i], &G[i], &S[i]);  if(err){return err;}
    }
    else if(indFLUIDE[i]==ind_air){
      EGS(tst_air, TAU[i], EPS[i], &P[i], &C[i]);
      G_EGS(tst_air, TAU[i], EPS[i], &G[i]);
      //printf("Init  i= %d | G=%g\n",i,G[i] );
    }

  }

  // variables défini aux frontières des mailles
  for(int i=0; i<=N; i++){
    x=X[i];

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, x-dx/2, x+dx/2, Wav);
    if(err){return err;}

    RHO0d[i]=Wav[0];
    U[i]=Wav[2];
    EKIN[i]=Wav[8]/RHO0d[i];
  }
  free(Wpw);
  free(Wav);

  
  // energie totale
  for(int i=ideb; i<=ifin; i++){ 
    E[i]=EPS[i] + (0.5)*(EKIN[i]+EKIN[i+1]);
  }


  // minimum de DXsC  (calcul de dt)
  minDXsC=dx/C[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if(minDXsC>dx/C[j]){
      minDXsC=dx/C[j];
    }
  }
  dt = dt_first*CFL * minDXsC; // pas de temps
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

    // Calcul du flux Lagrange
    //  U*
    for(int i=ideb; i<=ifin+1; i++){
      DxP=(P[i]-P[i-1])/dx;
      
      DtU=-(1./RHO0d[i])*DxP;

      if(sch==00){
        Ustar[i]=U[i] + zeta*(dt/2)*DtU ;// ordre 1;
      }
      else if(sch==01){
        rc2=0.5*(C[i]/TAU[i] + C[i-1]/TAU[i-1]);
        rc2=rc2*rc2;
        
        DxU=(U[i+1]-U[i-1])/(2*dx);
        DxxU=(U[i+1]+U[i-1]-2*U[i])/(dx*dx);
        
        rc2l=(C[i-1]/TAU[i-1])*(C[i-1]/TAU[i-1]);
        rc2r=(C[i]/TAU[i])*(C[i]/TAU[i]);
        DxRC2=(rc2r-rc2l)/dx;

        DxTAU0=(1./RHO0c[i] - 1./RHO0c[i-1])/dx;

        tau02=(1./(RHO0d[i]*RHO0d[i]));

        DttU=(1./RHO0d[i])*rc2*DxU*DxTAU0;
        DttU+=tau02*DxU*DxRC2;
        DttU+=tau02*rc2*DxxU;

        Ustar[i]=U[i] + (dt/2)*DtU + zeta*(dt*dt/6)*DttU; // ordre 2
      }

      // Pour print
      LdtU[i]=DtU;
      LdttU[i]=DttU;
    }
    //  P*
    for(int i=ideb-1; i<=ifin+1; i++){
      DxU=(U[i+1]-U[i])/dx;
      rc2=(C[i]/TAU[i])*(C[i]/TAU[i]);
      
      DtP=-(1./RHO0c[i])*(rc2)*DxU;
      
      if(sch==00){
        Pstar[i]=P[i] + zeta*(dt/2)*DtP; // ordre 1
      }
      else if(sch==01){
        DxP=(P[i+1]-P[i-1])/(2*dx);
        DxxP=(P[i+1]+P[i-1]-2*P[i])/(dx*dx);

        tau02=(1./RHO0c[i])*(1./RHO0c[i]);
        tau03=tau02/RHO0c[i];

        DxRHO0=(RHO0d[i+1]-RHO0d[i])/dx;

        DttP=tau02*rc2*DxxP;
        DttP+=-tau03*rc2*DxP*DxRHO0;
        DttP+=2*tau02*(1./TAU[i])*rc2*G[i]*(DxU*DxU);
        
        // printf(" i=%d | DttP= %g + %g + %g\n",i,tau02*rc2*DxxP,-tau03*rc2*DxP*DxRHO0, 2*tau02*(1./TAU[i])*rc2*G[i]*(DxU*DxU) );
        // printf("        (1) tau02=%g rc2=%g DxxP=%g \n",tau02,rc2,DxxP);
        // printf("            DxxP=(%g-%g)/(%g + %g)\n",DxPr,DxPl,dXcl,dXcr);
        // printf("        (2) tau03=%g rc2=%g DxP=%g DxRHO0=%g \n",tau03,rc2,DxP,DxRHO0);

        Pstar[i]=P[i] + (dt/2)*DtP + zeta*(dt*dt/6)*DttP; // ordre 2
      }
      // Pour print 
      LdtP[i]=DtP;
      LdttP[i]=DttP;
    }
    

    // Mise à jour de la solution
    for(int i=ideb; i<=ifin; i++){
      TAU[i]+=(dt/(RHO0c[i]*dx))*(Ustar[i+1]-Ustar[i]);
      EPS[i]-=(dt/(RHO0c[i]*dx))*Pstar[i]*(Ustar[i+1]-Ustar[i]);
      //printf("i= %d : TAU=%g, U=%g, E=%g \n",i,TAU[i],U[i],E[i] );
    }
    for(int i=ideb; i<=ifin+1; i++){
      U[i]-=(dt/(RHO0d[i]*dx))*(Pstar[i]-Pstar[i-1]);
      EKIN[i]-=(dt/(RHO0d[i]*dx))*Ustar[i]*(Pstar[i]-Pstar[i-1]);
    }
    
    
    // Mise à jour du maillage
    for(int i=ideb; i<=ifin+1; i++){
      X[i]+=U[i]*dt;
      //X[i]+=(Ustar[i]-U[i])*dt; ///RHO0d[i];
    }
    

    // Kinetic energy fix
    for (int i=ideb; i<=ifin+1; i++){
      DK[i]=RHO0d[i]*(EKIN[i]-0.5*U[i]*U[i]);
      EKIN[i]-=DK[i]/RHO0d[i];
      //EKIN[i]=0.5*U[i]*U[i];
    }
    for (int i=ideb; i<=ifin; i++){
      EPS[i]+=0.5*(DK[i]+DK[i+1])/RHO0c[i];
    }
    

    // Energie totale
    for(int i=ideb; i<=ifin; i++){ 
      E[i]=EPS[i] + (0.5)*(EKIN[i]*RHO0d[i]+EKIN[i+1]*RHO0d[i+1])/RHO0c[i];
    }


    // cond aux limites
    err=condlim(ind_cond, ideb,ifin, nbghosts, TAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, EPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, U, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, EKIN, 0);
    if(err){return err;}
    
    for(int j=0; j<nbghosts-1; j++){
      // à gauche
      X[ideb-1-j]=X[ideb]-(X[ideb+1+j]-X[ideb]);
      // à droite
      X[ifin+2+j]=X[ifin+1]+(X[ifin+1]-X[ifin-j]);
    }
    for(int i=0; i<=N-1; i++){
      Xc[i]=(X[i+1]+X[i])/2;
    }
    
    
    // Mise à jour des tableaux
    for(int i=0; i<=N-1; i++){ 
      dUdx=(U[i+1]-U[i-1])/(2*dx);

      if(indFLUIDE[i]==ind_etain){
        err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], matFM[i], 
                     Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
                     nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                     SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){printf(" cycle=%d (funcKOVAo2_multi)\n",n); return err;}
        
        G[i]=pVAR[0];
        S[i]=pVAR[1];
        //err= G_S_cin_3_phase(epsilon, P[i], TAU[i], EPS[i], matFM[i], &C[i], &G[i], &S[i]); // S n'est pas utilisée
        //if(err){return err;}

        //printf(" xA= %g, xB= %g, xC= %g\n",matFM[i][0],matFM[i][1],matFM[i][2] );
      }
      else if(indFLUIDE[i]==ind_air){
        EGS(tst_air, TAU[i], EPS[i], &P[i], &C[i]);
        G_EGS(tst_air, TAU[i], EPS[i], &G[i]);

        //EGS_cin_phase(tst_air,  TAU[i], EPS[i], &P[i], &C[i], cin_phase, dx, dt, Ttau, &Un[i], P[i], TAUn[i], C[i]);
        //G_BIZ_cin_phase( TAU[i], EPS[i], &G[i], cin_phase, dx, dt, Ttau, &Un[i], P[i], TAUn[i], C[i]);
      }
      
      if(C[i]==0){C[i]=Ctriple; }
      //printf("i=%d | TAU=%g, U=%g, E=%g, EPS=%g, C=%g, P=%g \n",i,TAU[i],U[i],E[i],EPS[i],C[i],P[i] );
    }
    
    /*
    for(int i=ideb; i<=ifin; i++){
     printf("i=%d | TAU=%g U=%g EPS=%g EKIN=%g P=%g   \n",i,TAU[i],U[i],EPS[i],EKIN[i],P[i]);
     printf("        U*=%g DtU=%g DttU=%g   zeta=%g\n",Ustar[i],LdtU[i],LdttU[i],zeta);
     printf("        P*=%g DtP=%g DttP=%g\n",Pstar[i],LdtP[i],LdttP[i]);
     printf("        Xf=%g U=%g dt=%g\n",X[i],U[i],dt);
    }
    printf("i=%d | U=%g EKIN=%g \n",ifin+1,U[ifin+1],EKIN[ifin+1] );
    printf("        U*=%g DtU=%g DttU=%g\n",Ustar[ifin+1],LdtU[ifin+1],LdttU[ifin+1]);
    printf("        P*=%g DtP=%g DttP=%g\n",Pstar[ifin+1],LdtP[ifin+1],LdttP[ifin+1]);
    printf("        Xf=%g U=%g dt=%g\n",X[ifin+1],U[ifin+1],dt);
    */


    // Energie du système
    sum1=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RHO0c[i]*E[i];
      sum3+=dx*RHO0d[i]*U[i];
      sum4+=(X[i+1]-X[i])/TAU[i];
      sum5+=dx*RHO0c[i]*TAU[i];
    }
    sum3+=dx*RHO0d[ifin+1]*U[ifin+1];
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;

    t=t+dt;  n=n+1;
    if(aff==1){ printf("n= %d, t= %g, dt= %g\n",n,t,dt); }


    // Calcul du pas de temps
    int indice=ideb;
    dX=X[ideb+1]-X[ideb];
    //dX=dx;
    minDXsC=dX/C[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      dX=X[j+1]-X[j];
      //dX=dx;
      if(minDXsC>dX/C[j]){  minDXsC=dX/C[j]; indice=j; }
    }
    //printf(" min(dX/C)=%g  indice=%d  C=%g  dx=%g dU=%g\n", minDXsC, indice, C[indice], X[indice+1]-X[indice], U[indice+1]-U[indice] );
    dt_new = CFL * minDXsC; // pas de temps
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

  FILE* Res;
  if((Res = fopen("nb_iter.txt", "w+")) != NULL){
    for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%d ",tabNiter[i]);  }
    fprintf(Res, "0\n");
    for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%g ",tabEPSILON[i]);  }
    fprintf(Res, "-1e-10\n");
    fclose(Res);
  } 
  else{printf("erreur impression\n"); return 1;}


  free(X); free(Xc); free(indFLUIDE);
  free(RHO0d); free(RHO0c);
  
  free(TAU); free(U); free(EPS); free(EKIN);
  free(Ustar);  free(Pstar);

  free(P); free(C); 
  free(E); free(RHO);

  free(G); free(S);
  
  freetab(matFM); free(T); free(lambda_out_unuse);
  //  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  // Cinétique
  free(TAUn); free(Un);

  free(tabNiter);
  free(tabEPSILON);
  free(test_consis_thermo);
  
  free(LdtU);
  free(LdttU);
  free(LdtP);
  free(LdttP);


  return 0;
  
}
