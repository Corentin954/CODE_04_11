#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_multi.h"




// SCHEMA DE TYPE RUNGE-KUTTA FORMULE EN ENERGIE INTERNE

/////////////////////  MAIN  ////////////////////
int funcRKint_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int kin_fix, int dpi, int so, int EOS,
                    int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                    double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                    double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
  
  double Gamma=7./5;
  
  int ind_air, ind_etain; 
  ind_air=0; ind_etain=1;

  if(kin_fix!=0 && kin_fix!=1){printf("Erreur kin_fix= %d (kinetic energy fix)\n",kin_fix ); return 1;}
  if(dpi!=0 && dpi!=1){printf("Erreur dpi= %d\n",dpi ); return 1;}

	// variables
	double t=0;
	double dt, dt_new;
	int n=0, err=0;
  double maxC; 
	char* W0;
  double sum1, sum2, sum3, sum4, sum5;
  int N;
  int nbghosts=2*(so+1); // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts2=(so+1); // nombre de mailles fantomes pour les fonctions qui dependent des variables ci-dessus
  // l'idee  : appliquer les ocnd aux lim sur les var. sur nbghosts mailles qui calculer la valeur des fcts sur les nbghosts2 mailles
  int ideb, ifin;
  double** Ck; double** Ckbar; double** dk; double** Qbar; double** rk;
  int sodPi; // spatial order for dPi
  double dX;
  int zone;
  double dUdx;
  double* pVAR=malloc(2*sizeof(double));

  // spatial order for dPi
  sodPi=so;


  // Coefficient
  Ck=alloctabd(5,5); Ckbar=alloctabd(5,5);
  dk=alloctabd(5,5); Qbar=alloctabd(5,5);
  rk=alloctabd(5,5);
  initCoef(Ck, Ckbar, dk, Qbar, rk);

  
  // Coefficient de Btucher
  int to;  // time order
  double** A; double* THETA; double* ALPHA;
  // RK1 et RK2
  if(sch==0||sch==1){   to=1;  }
  // SSPRK3(3,3) Spiteri-Ruuth
  else if (sch==2){   to=2;  }
  // RK4 Classic
  else if (sch==3){   to=3;  }
  // RK5 Cash-Karp 
  else if (sch==4){   to=5;  }
  // Dormand-Prince 
  else if (sch==5){   to=6;  }
  // HEUN  ordre 3
  else if (sch==6){   to=2;  }
  // RK3 RALSTON
  else if (sch==7){   to=2;  }
  // RK2 ou RK3 Bogaki-Shampine
  else if (sch==8 || sch==9){    to=3;  }
  // KUTTA odre 3      
  else if (sch==10){   to=2;  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){   to=1;  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){   to=1;  }
  // EULER forward
  else if (sch==13){   to=0;  }
  // SSP RK4 OPTIMUM
  else if (sch==14){   to=4;  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){   to=3;  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){   to=4;  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){   to=2;  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){   to=2;  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){   to=3;  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){   to=1;  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){   to=3;  }
  else{ printf("Erreur dans le choix de sch (RKint)\n"); return 1;}
  

  A=alloctabd(to, to);
  THETA=malloc((to+1)*sizeof(double));
  ALPHA=malloc(to*sizeof(double));
	initButcher(sch, A, THETA,ALPHA);
  /*
  for(int i=0; i<to; i++){
    for(int j=0; j<=i; j++){
      printf("A[%d][%d]= %lf ",i,j,A[i][j] );
    }
    printf("\n");
  }
  for (int i = 0; i < to+1; ++i){
    printf("theta[%d]= %lf ",i,THETA[i] );
  }
  printf("\n");
  */

  double dx=(b-a)/nx; // pas d'espace


  // Calcul de la taille des tableaux
  N=nx+2*nbghosts;
  ideb=nbghosts;
  ifin=nbghosts+nx-1;
  int N2=nbghosts+nx+nbghosts2; // cond au lim de Q et P sur i=nbghosts- nbghosts2, N2
  int ideb2=nbghosts-nbghosts2; 
  int ifin2=N2-1; 


  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double)); // maillage
  double* Xc=malloc((N)*sizeof(double)); // maillage
  // PW
  double* RHO0c=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
  double* RHO0d=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
	
  double* RHO0cpw=malloc((N)*sizeof(double));   //  rho0  (centrées aux mailles)
	double* RHO0dpw=malloc((N+1)*sizeof(double)); //  rho0  (dé-centrées | frontières)
   
	// solution AVERAGE
	double* RTAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* RU=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* REPS=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)
  double* REKIN=malloc((N+1)*sizeof(double));

  // solution POINT-WISE
  double* RTAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* RUpw=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* REPSpw=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)
  double* REKINpw=malloc((N+1)*sizeof(double));

  double* TAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* Upw=malloc((N+1)*sizeof(double)); //  rho0.u (frontières)
  double* EPSpw=malloc((N)*sizeof(double)); //  rho0.epsilon (centrées aux mailles)

  // AVERAGE (plot)  
  double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* EPS=malloc((N)*sizeof(double)); //  epsilon (centrées aux mailles)  
  double* U=malloc((N+1)*sizeof(double)); //  u  (frontières aux mailles)
  double* RHO=malloc((N)*sizeof(double));   //  p (centrées aux mailles)  
  double* E=malloc((N)*sizeof(double));   //  p (centrées aux mailles)  

  // Point-Wise
  double* Epw=malloc((N)*sizeof(double));
  double* RHOpw=malloc((N)*sizeof(double));
  double* Ppw=malloc((N)*sizeof(double));
  double* Qpw=malloc((N)*sizeof(double));   // artificial viscosity
  double* Cpw=malloc((N)*sizeof(double));
  
  double* DM0c=malloc((N)*sizeof(double));  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles

  double* dPI=malloc((N)*sizeof(double)); // delta(P+Q)
  double** dPIstar=alloctabd(2,N); // delta(P+Q)
  
  // tableaux temporaires
  // AVERAGE
  double** RTAUstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** RUstar=alloctabd(to,N+1); //  rho0.u (frontières)
  double** REPSstar=alloctabd(to,N); //  rho0.tau (centrées aux mailles)
  double** REKINstar=alloctabd(to,N+1);

  // Point-Wise
  double** Pstar=alloctabd(to,N);
  double** Qstar=alloctabd(to,N);  // artificial viscosity
  double** Cstar=alloctabd(to,N);
  double** Ustar=alloctabd(to,N+1); // vitesse donc frontières aux mailles
  double** TAUstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** EPSstar=alloctabd(to,N); //  rho0.tau  (centrées aux mailles)
  double** DMstar=alloctabd(to,N);  // Delta m = (X[i+1]-X[i])*(RHOC[i]/RTAU[i])   centres aux mailles
  double** Xstar=alloctabd(to,N+1); // maillage

  // Cinétique
  double* TAUpw_n=malloc((N)*sizeof(double));   
  double* Upw_n=malloc((N+1)*sizeof(double)); 

  // kinetic fix  
  double* DK=malloc((N+1)*sizeof(double));
  double* R2U2pw=malloc((N+1)*sizeof(double)); //(point-wise)
  double* RRpw=malloc((N)*sizeof(double)); //(point-wise)

  double* RU2pw=malloc((N+1)*sizeof(double)); //(point-wise)
  
  double* PVg=malloc((N)*sizeof(double)); // Compression adiabatique

  double* Tpw=malloc((N)*sizeof(double));  // température 
  double* FM=malloc((N)*sizeof(double)); // fraction massique
  
  double** matFM=alloctabd(N,3);          // fraction massique
  double* lambda_out_unuse=malloc(3*sizeof(double));
  
  double* Gpw=malloc((N)*sizeof(double)); // fraction massique
  double* Spw=malloc((N)*sizeof(double)); // fraction massique

  int* tabNiter=malloc((N)*sizeof(int)); // nb iterration par mailles
  for(int i=0; i<N; ++i){tabNiter[i]=0;}

  double* tabEPSILON=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-1e-15;}

  double* test_consis_thermo=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-99;}

  int* indFLUIDE=malloc((N)*sizeof(int)); // indice du fluide sur chaque maille
  
  // Maillage
  for(int i=0; i<=N; i++){
    X[i]=(a-nbghosts*dx)+dx*i;
    //printf("i= %d, X=%lf\n",i,X[i]);
  }

  err=init_fluides(tst_multi, indFLUIDE, N, ind_air, ind_etain);
  if(err){printf("Erreur fluide\n"); return 1;}

  // Initialisation POINT-WISE
  double x;
  // variables défini aux centres des mailles
  double* Wpw=malloc(10*sizeof(double));
  double* Wav=malloc(10*sizeof(double));
  for(int i=0; i<=N-1; i++){
    //printf("  --i= %d\n", i);
    x=(X[i]+X[i+1])/2.0;

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, X[i],X[i+1], Wav);
    if(err){return err;}
    
    RHO0c[i]=Wav[0];
    DM0c[i]=dx*RHO0c[i];
    
    RTAU[i]=1.0;
    REPS[i]=Wav[7];
    
    RTAUpw[i]=1.0;
    
    EPS[i]=Wav[4];
    
    Ppw[i]=Wpw[2];
    Tpw[i]=Wpw[3];  
    matFM[i][0]=Wpw[4];
    matFM[i][1]=Wpw[5];
    matFM[i][2]=Wpw[6];

    //printf("Init i=%d | x=%g lambda=[%g,%g,%g], eps=%g, P=%g\n",i,x,matFM[i][0],matFM[i][1],matFM[i][2],EPS[i],Ppw[i] );

    PVg[i]=Ppw[i]/pow(RHO0c[i],Gamma);
  }

  // variables défini aux frontières des mailles
  for(int i=0; i<=N; i++){
    x=X[i];

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, x-dx/2, x+dx/2, Wav);
    if(err){return err;}
    
    RHO0d[i]=Wav[0];
    
    RU[i]=Wav[5];
    U[i]=Wav[2];

    REKIN[i]=Wav[8];
  }
  free(Wpw);
  free(Wav);

  for(int i=0; i<=N-1; i++){
    RHO0cpw[i]=phi(so,&RHO0c[i],Ck);
    RHO0dpw[i]=phi(so,&RHO0d[i],Ck);
    TAUpw[i]=1./RHO0cpw[i];
    RUpw[i]=phi(so,&RU[i],Ck);
    Upw[i]=phi(so,&U[i],Ck);
    REPSpw[i]=phi(so,&REPS[i],Ck);
    EPSpw[i]=phi(so,&EPS[i],Ck);

    if(indFLUIDE[i]==ind_etain){
      err= G_S_cin_3_phase(epsilon, Ppw[i], TAUpw[i], EPSpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]);  if(err){return err;}
    }
    else if(indFLUIDE[i]==ind_air){
      EGS(tst_air, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i]);
      G_EGS(tst_air, TAUpw[i], EPSpw[i], &Gpw[i]);
    }
  }
  RUpw[N]=phi(so,&RU[N],Ck);
  Upw[N]=phi(so,&U[N],Ck);



  //for(int i=nbghosts-nbghosts2; i<=N2-1; i++){
  for(int i=ideb2; i<=ifin2; i++){
    Qpw[i]=qvis(q, Cq, Cl, &Upw[i], TAUpw[i], Cpw[i], &DM0c[i]);
  }
  
  /*
  for(int i=0; i<=N-1; i++){
    printf("Init i=%d x=%g  Point-Wise |  tau=%g, eps=%g, P=%g, u=%g, q=%g, c=%g\n",i,X[i],TAUpw[i],EPSpw[i],Ppw[i],Upw[i],Qpw[i],Cpw[i] );
    printf("                     Average    |  eps=%g, u=%g, REPS[i]= %g, RU[i]= %g\n",EPS[i],U[i],REPS[i],RU[i] );
  }
  */


  // Calcul du pas de temps
  maxC=Cpw[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    if (maxC<Cpw[j]){  maxC=Cpw[j];  }
  }
  dt= dt_first *CFL *dx/maxC; // pas de temps
  
  // energie totale du système
  int nbiter=5*ceil(Tf/dt);
  //printf("nbiter :%d\n",nbiter);  
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double)); // masse totale
  double* VOL=malloc(nbiter*sizeof(double)); // volume total
  
  // Energie du système
  //for(int i=nbghosts-nbghosts2; i<=N2; i++){
  for(int i=ideb2; i<=ifin2+1; i++){
    R2U2pw[i]=RUpw[i]*RUpw[i];  // (rho.u)^2
    RRpw[i]=(X[i+1]-X[i])*RHO0c[i]/RTAUpw[i];
  }
 
  sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*(REPS[i] + REKIN[i]);
    sum3+=dx*RU[i];
    sum4+=phi(so, &RRpw[i], Ckbar);
    sum5+=dx*RTAU[i];
  }
  sum1+=dx*REKIN[ifin+1];
  sum3+=dx*RU[ifin+1];
  Etot[0]=sum1;
  IMPUL[0]=sum3; 
  Mtot[0]=sum4;
  VOL[0]=sum5;

  
// debut de la boucle =============================================================================
	while((t<Tf)&&(n<Nmax)){

    // cond aux limites  symétrie    
    err=condlim(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, RU, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // attention frontières aux mailles
    if(err){return err;}
  
    // Cinétique
    for(int i=0; i<=N-1; i++){
      TAUpw_n[i]=TAUpw[i];
      Upw_n[i]=Upw[i];
    }
    Upw_n[N]=Upw[N];

    // dPI=delta(P+Q)
    if(dpi==1){
      //for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ ////
      for(int i=ideb2; i<=ifin2; i++){ ////
        dPI[i]=delta(sodPi,0,&Ppw[i], dk)+delta(sodPi,0,&Qpw[i], dk);
      }
    }
    
	  // Boucle RUNGE-KUTTA
    for(int k=0; k<to; k++){
      // ordre en temps
      for(int i=ideb; i<=ifin; i++){
        RTAUstar[k][i]=RTAU[i]+A[k][0]*(dt/dx)*(Upw[i+1]-Upw[i]);
        if(dpi==0){  RUstar[k][i]=RU[i]-A[k][0]*(dt/dx)*((Ppw[i]+Qpw[i])-(Ppw[i-1]+Qpw[i-1]));  }
        else if(dpi==1){  RUstar[k][i]=RU[i]-A[k][0]*(dt/dx)*phi(sodPi, &dPI[i], Ckbar);  }
        REPSstar[k][i]=REPS[i]-A[k][0]*(dt/dx)*YdZ(so,0,&Ppw[i],&Qpw[i],&Upw[i],dk,Ckbar);
        REKINstar[k][i]=REKIN[i]-A[k][0]*(dt/dx)*YdZ(so,1,&Ppw[i],&Qpw[i],&Upw[i],dk,Ckbar);
      }
      if(dpi==0){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*((Ppw[ifin+1]+Qpw[ifin+1])-(Ppw[ifin]+Qpw[ifin]));  }
      else if(dpi==1){  RUstar[k][ifin+1]=RU[ifin+1]-A[k][0]*(dt/dx)*phi(sodPi, &dPI[ifin+1], Ckbar);  }
      REKINstar[k][ifin+1]=REKIN[ifin+1]-A[k][0]*(dt/dx)*YdZ(so,1,&Ppw[ifin+1],&Qpw[ifin+1],&Upw[ifin+1],dk,Ckbar);

      //for(int i=nbghosts - nbghosts2; i<=N2; i++){
      for(int i=ideb2; i<=ifin2+1; i++){
        Xstar[k][i]= X[i] + A[k][0]*dt*Upw[i];
      }

      
      for(int j=1; j<=k; j++){ // calcul des sous pas de temps
        //
        for(int i=ideb; i<=ifin; i++){
    	    RTAUstar[k][i]+=A[k][j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
          if(dpi==0){   RUstar[k][i]-=A[k][j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
          else if(dpi==1){    RUstar[k][i]-=A[k][j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][i], Ckbar);   }
          REPSstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar);
          REKINstar[k][i]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar);

    	  }
        if(dpi==0){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
        else if(dpi==1){   RUstar[k][ifin+1]-=A[k][j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][ifin+1], Ckbar);  }
        REKINstar[k][ifin+1]-=A[k][j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1], &Qstar[j-1][ifin+1], &Ustar[j-1][ifin+1],dk,Ckbar);

        for(int i=ideb2; i<=ifin2+1; i++){
          Xstar[k][i]+=A[k][j]*dt*Ustar[j-1][i];
        }

      }

      // cond aux limites
      err=condlim(ind_cond, ideb,ifin, nbghosts, RTAUstar[k], 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin, nbghosts, REPSstar[k], 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, RUstar[k], iu);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKINstar[k], 0);
      if(err){return err;}


      // Energy kinetic fix
      // CRAS 2016
      if(kin_fix==11 || kin_fix==12){ 
        for(int i=ideb2; i<=ifin2+1; i++){ 
          RUpw[i]=phi(so, &RUstar[k][i],Ck);
        }
        if(kin_fix==11){
          for(int i=ideb2; i<=ifin2+1; i++){  ////  
            REKINpw[i]=phi(so, &REKINstar[k][i],Ck);
          }
          for(int i=ideb2; i<=ifin2+1; i++){ 
            DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0dpw[i]);
          }
          for (int i=ideb; i<=ifin; i++){
            REKINstar[k][i]-=phi(so, &DK[i],Ckbar);
            REPSstar[k][i]+=phiQK(so, &DK[i],Qbar); 
          }
          REKINstar[k][ifin+1]-=phi(so, &DK[ifin+1],Ckbar);
        }
        // JCP 2019 Dakin & al.
        else if(kin_fix==12){
          for(int i=ideb2; i<=ifin2+1; i++){  ////
            RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0dpw[i]);
          }
          for (int i=ideb; i<=ifin+1; i++){
            DK[i]=REKINstar[k][i]-phi(so, &RU2pw[i], Ckbar);
          } 
          for (int i=ideb; i<=ifin; i++){
            REKINstar[k][i]-= DK[i];
            REPSstar[k][i]+= (DK[i]+DK[i+1])/2.0; 
          }
          REKINstar[k][ifin+1]-= DK[ifin+1]; 
        }

        // On recalcule les variables dépendant de REPS
        // cond aux limites  symétrie  
        err=condlim(ind_cond, ideb, ifin, nbghosts, REPSstar[k], 0);
        if(err){return err;}
        err=condlim(ind_cond, ideb, ifin, nbghosts, REKINstar[k], 0);
        if(err){return err;} 
      }


      // On veut que U soit PW or PUstar est AV
      for(int i=ideb2; i<=ifin2; i++){  ////
        Ustar[k][i]=phi(so, &RUstar[k][i], Ck)/RHO0dpw[i];
        TAUstar[k][i]=phi(so, &RTAUstar[k][i], Ck)/RHO0cpw[i];
        EPSstar[k][i]=phi(so, &REPSstar[k][i], Ck)/RHO0cpw[i];
      }
      Ustar[k][ifin2+1]=phi(so, &RUstar[k][ifin2+1], Ck)/RHO0dpw[ifin2+1];

  	  for(int i=ideb2; i<=ifin2; i++){  ////
        dX=Xstar[k][i+1]-Xstar[k][i];
        dUdx=(Upw[i+1]-Upw[i])/dX;

        if(indFLUIDE[i]==ind_etain){
          err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAUstar[k][i], EPSstar[k][i], &Pstar[k][i], &Cstar[k][i], &Tpw[i], pVAR, matFM[i], lambda_out_unuse, 
                 Ntau, dX, A[k][k]*dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, Ppw[i], TAUpw[i], Cpw[i],
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                 SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
          if(err){return err;}
        }
        else if(indFLUIDE[i]==ind_air){
          //EGS(tst_air, TAUstar[k][i], EPSstar[k][i], &Pstar[k][i], &Cstar[k][i]); // On a besoin d'argument PW

          EGS_cin_phase(tst_air,  TAUstar[k][i], EPSstar[k][i], &Pstar[k][i], &Cstar[k][i],
                        cin_phase, dX, A[k][k]*dt, Ttau,
                        &Upw[i], Ppw[i], TAUpw[i], Cpw[i]);

          //if(tst==3){  Pstar[k][i]=PVg[i]*pow(1.0/TAUstar[k][i], Gamma);  }
        }
        Qstar[k][i]=qvis(q, Cq, Cl, &Ustar[k][i], TAUstar[k][i], Cstar[k][i], &DM0c[i]);
  	  }

      // dPIstar
      if(dpi==1){
        for(int i=ideb2; i<=ifin2; i++){ 
          dPIstar[k][i]=delta(sodPi,0,&Pstar[k][i], dk)+delta(sodPi,0,&Qstar[k][i], dk);
        }
      }

    } // fin de boucle sur k ///////////


    // Solution au temps n+1
    for(int i=ideb; i<=ifin; i++){
      RTAU[i]+=THETA[0]*(dt/dx)*(Upw[i+1]-Upw[i]);
      if(dpi==0){  RU[i]-=THETA[0]*(dt/dx)*((Ppw[i]+Qpw[i])-(Ppw[i-1]+Qpw[i-1]));  }
      else if(dpi==1){  RU[i]-=THETA[0]*(dt/dx)*phi(sodPi, &dPI[i], Ckbar);  }
      REPS[i]-=THETA[0]*(dt/dx)*YdZ(so,0,&Ppw[i],&Qpw[i],&Upw[i],dk,Ckbar);
      REKIN[i]-=THETA[0]*(dt/dx)*YdZ(so,1,&Ppw[i],&Qpw[i],&Upw[i],dk,Ckbar);    
    }
    if(dpi==0){  RU[ifin+1]-=THETA[0]*(dt/dx)*((Ppw[ifin+1]+Qpw[ifin+1])-(Ppw[ifin]+Qpw[ifin]));  }
    else if(dpi==1){  RU[ifin+1]-=THETA[0]*(dt/dx)*phi(sodPi, &dPI[ifin+1], Ckbar);  }
    REKIN[ifin+1]-=THETA[0]*(dt/dx)*YdZ(so,1,&Ppw[ifin+1],&Qpw[ifin+1],&Upw[ifin+1],dk,Ckbar);

    //for(int i=nbghosts-nbghosts2; i<=N2; i++){
    for(int i=ideb2; i<=ifin2+1; i++){
      X[i]+=dt*THETA[0]*Upw[i];
    }


    for (int j=1; j<=to; j++){ // calcul des sous pas de temps
      for(int i=ideb; i<=ifin; i++){
        RTAU[i]+=THETA[j]*(dt/dx)*(Ustar[j-1][i+1]-Ustar[j-1][i]);
        if(dpi==0){   RU[i]-=THETA[j]*(dt/dx)*((Pstar[j-1][i]+Qstar[j-1][i])-(Pstar[j-1][i-1]+Qstar[j-1][i-1]));   }
        else if(dpi==1){    RU[i]-=THETA[j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][i], Ckbar);   }
        REPS[i]-=THETA[j]*(dt/dx)*YdZ(so,0,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar); 
        REKIN[i]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][i],&Qstar[j-1][i],&Ustar[j-1][i],dk,Ckbar);  
      }
      if(dpi==0){   RU[ifin+1]-=THETA[j]*(dt/dx)*((Pstar[j-1][ifin+1]+Qstar[j-1][ifin+1])-(Pstar[j-1][ifin]+Qstar[j-1][ifin]));  }
      else if(dpi==1){   RU[ifin+1]-=THETA[j]*(dt/dx)*phi(sodPi, &dPIstar[j-1][ifin+1], Ckbar);  }
      REKIN[ifin+1]-=THETA[j]*(dt/dx)*YdZ(so,1,&Pstar[j-1][ifin+1],&Qstar[j-1][ifin+1],&Ustar[j-1][ifin+1],dk,Ckbar);  
      
      //for(int i=nbghosts-nbghosts2; i<=N2; i++){ // permet la gestion des cond aux lim du maillage
      for(int i=ideb2; i<=ifin2+1; i++){ // permet la gestion des cond aux lim du maillage
        X[i]+=dt*THETA[j]*Ustar[j-1][i];
      }
    }
  

    // cond aux limites  symétrie
    err=condlim(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, RU, iu);
    if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // ATENTION frontières aux mailles
    if(err){return err;}

    
    // Mise a jour des rho.Z  Point-Wise
    //for(int i=nbghosts-nbghosts2; i<=N2-1; i++){ //// 
    for(int i=ideb2; i<=ifin2; i++){ //// 
      RTAUpw[i]=phi(so, &RTAU[i], Ck);
      REPSpw[i]=phi(so, &REPS[i], Ck);
      RUpw[i]=phi(so, &RU[i], Ck);
    }
    RUpw[ifin2+1]=phi(so, &RU[ifin2+1], Ck);


    // Energy kinetic fix
    // CRAS 2016
    if(kin_fix==1 || kin_fix==2 || kin_fix==11 || kin_fix==12){ 
      if(kin_fix==1 || kin_fix==11){
        for(int i=ideb2; i<=ifin2+1; i++){ ////  
          REKINpw[i]=phi(so, &REKIN[i],Ck);
        }
        for (int i=ideb2; i<=ifin2+1; i++){
          DK[i]=REKINpw[i]-RUpw[i]*RUpw[i]/(2.0*RHO0dpw[i]);
        }
        for (int i=ideb; i<=ifin; i++){
          REKIN[i]-=phi(so, &DK[i],Ckbar);
          REPS[i]+=phiQK(so, &DK[i],Qbar); 
        }
        REKIN[ifin+1]-=phi(so,&DK[ifin+1],Ckbar);
      }
      // JCP 2019 Dakin & al.
      else if(kin_fix==2 || kin_fix==12){
        for(int i=ideb2; i<=ifin2+1; i++){ ////
          RU2pw[i]=RUpw[i]*RUpw[i]/(2.0*RHO0dpw[i]);
        }
        for (int i=ideb; i<=ifin+1; i++){
          DK[i]=REKIN[i]-phi(so, &RU2pw[i], Ckbar);
        } 
        for (int i=ideb; i<=ifin; i++){
          REKIN[i]-=DK[i];
          REPS[i]+=(DK[i]+DK[i+1])/2.0;
        }
        REKIN[ifin+1]-= DK[ifin+1]; 
      }

      // On recalcule les variables dépendant de REPS
      // cond aux limites  symétrie  
      err=condlim(ind_cond, ideb,ifin, nbghosts, REPS, 0);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, REKIN, 0); // ATENTION frontières aux mailles
      if(err){return err;} 
      
      // Mise a jour des rho.Z  Point-Wise
      for(int i=ideb2; i<=ifin2; i++){  
        REPSpw[i]=phi(so, &REPS[i], Ck);
      } 
    } //// fin kin fix

    for(int i=ideb2; i<=ifin2+1; i++){  Upw[i]=RUpw[i]/RHO0dpw[i];  }

    // Mise à jour des tableaux
    for(int i=ideb2; i<=ifin2; i++){ //// 
      dX=X[i+1]-X[i];
      dUdx=(Upw_n[i+1]-Upw_n[i])/dX;

      if(indFLUIDE[i]==ind_etain){
        err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, RTAUpw[i]/RHO0cpw[i], REPSpw[i]/RHO0cpw[i], &Ppw[i], &Cpw[i], &Tpw[i], pVAR, matFM[i], matFM[i], 
                              Ntau, dX, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, Ppw[i], TAUpw_n[i], Cpw[i],
                              nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                              SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
        if(err){return err;}
      }
      else if(indFLUIDE[i]==ind_air){
        //EGS(tst_air, RTAUpw[i]/RHO0cpw[i], REPSpw[i]/RHO0cpw[i], &Ppw[i], &Cpw[i]);
        
        EGS_cin_phase(tst_air, RTAUpw[i]/RHO0cpw[i], REPSpw[i]/RHO0cpw[i], &Ppw[i], &Cpw[i],
                        cin_phase, dX, dt, Ttau,
                        &Upw_n[i], Ppw[i], TAUpw_n[i], Cpw[i]);

        //if(tst==3){  Ppw[i]=PVg[i]*pow(RHO0cpw[i]/RTAUpw[i], Gamma);  }
      }
    }


    for(int i=ideb2; i<=ifin2; i++){  
      Qpw[i]=qvis(q,Cq,Cl,&Upw[i], RTAUpw[i]/RHO0c[i], Cpw[i],&DM0c[i]);  
    }

    ///////////////
    // Energie du système
    for(int i=ideb2; i<=ifin2+1; i++){
      R2U2pw[i]=RUpw[i]*RUpw[i];  // (rho.u)^2
      RRpw[i]=(X[i+1]-X[i])*RHO0cpw[i]/RTAUpw[i];
    }
    sum1=0; sum2=0; sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){  
      sum1+=dx*(REPS[i] + REKIN[i]);
      sum3+=dx*RU[i];
      sum4+=phi(so, &RRpw[i], Ckbar);
      sum5+=dx*RTAU[i];
    }
    sum1+=dx*REKIN[ifin+1];
    sum3+=dx*RU[ifin+1];
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;

	  t=t+dt;
	  n=n+1;
    if(aff==1){
	    printf("n= %d, t= %g, dt= %g\n",n,t,dt);
	  }

    // Calcul du pas de temps
    maxC=Cpw[ideb];
    for(int j=ideb+1; j<=ifin; j++){
      if (maxC<Cpw[j]){  maxC=Cpw[j];  }
    }
    dt_new = CFL *dx/maxC; // pas de temps
    dt=pastemps(Tf,t,dt_new,dt);

  }
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);
  // energie totale
  for(int i=ideb2; i<=ifin2; i++){ 
    Xc[i]=fxc(so, &X[i], rk);
    Epw[i]=REPSpw[i]/RHO0cpw[i]+(((RUpw[i]/RHO0dpw[i]+RUpw[i+1]/RHO0dpw[i+1])/2)*(RUpw[i]/RHO0dpw[i]+RUpw[i+1]/RHO0dpw[i+1])/2)/2; 
    RHOpw[i]=RHO0cpw[i]/RTAUpw[i];
    
    TAU[i]=phidiv(so, &RTAUpw[i], &RHO0cpw[i], Ckbar);
    EPS[i]=phidiv(so, &REPSpw[i], &RHO0cpw[i], Ckbar);
    U[i]=phidiv(so, &RUpw[i], &RHO0dpw[i], Ckbar);
  }
  U[ifin+1]=phidiv(so, &RUpw[ifin+1], &RHO0dpw[ifin+1], Ckbar);
  
  for(int i=ideb; i<=ifin; i++){
    if(indFLUIDE[i]==ind_etain){  G_S_cin_3_phase(epsilon, Ppw[i], RTAUpw[i]/RHO0cpw[i], REPSpw[i]/RHO0cpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]); }
    else if(indFLUIDE[i]==ind_air){ G_EGS(tst_air, RTAUpw[i]/RHO0cpw[i], REPSpw[i]/RHO0cpw[i], &Gpw[i]); }

    E[i]=phi(so, &Epw[i], Ckbar);
    RHO[i]=phi(so, &RHOpw[i], Ckbar);
  }

  // printf_sol
  err=print_sol(0, ideb, ifin, X, Xc, RHO, TAU, U, E, EPS, Ppw, Tpw, Cpw, Gpw, Spw, matFM);
  if(err){return err;}

  // printf_sol
  err=print_sol_RKint(0, ideb, ifin, X, Xc, RHO, TAU, U, E, EPS, Ppw, Tpw, Cpw, Gpw, Spw, RTAU, RU, REPS);
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

  free(dPI); freetab(dPIstar); 

  free(X);  free(Xc); free(E); 
  free(RHO0cpw); free(RHO0dpw);
  free(RHO0c); free(RHO0d); free(DM0c);   free(RHO); 

  free(RTAU); free(RU); free(REPS);
  free(TAU); free(U); free(EPS);

  free(RTAUpw); free(RUpw); free(REPSpw); 
  free(TAUpw); free(Upw); free(EPSpw); 

  freetab(RTAUstar); freetab(RUstar); freetab(REPSstar);
  freetab(TAUstar); freetab(Ustar); freetab(EPSstar); 
  freetab(Pstar); freetab(Qstar); freetab(Cstar); freetab(Xstar);

  free(Ppw); free(Qpw); free(Cpw); free(Gpw); free(Spw); free(Epw);

  free(Ck); free(Ckbar); free(dk); free(Qbar); free(rk);
  freetab(matFM); free(Tpw); free(lambda_out_unuse);

  // Energy kinetic fix
  free(R2U2pw); free(DK); free(RRpw);
  free(REKIN); free(REKINpw);
  freetab(REKINstar); 
  free(Etot); 
  free(RU2pw);
  free(IMPUL);
  free(Mtot);
  free(VOL);

  free(TAUpw_n);   
  free(Upw_n); 
  free(PVg); // Compression adiabatique
  free(FM); // fraction massique
  free(tabNiter); // nb iterration par mailles
  free(tabEPSILON); // valeur du critère
  free(indFLUIDE); // indice du fluide sur chaque maille
  

  return 0;
  
}
