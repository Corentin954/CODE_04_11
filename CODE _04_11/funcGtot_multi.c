#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "head_multi.h"
//#include "head.h"


/*
  Schéma lagrangien de type Godunov 
*/

int sgn(double v);
int choix_sch(int sch, int* psolveur, int* pmethd, int* plimiteur);


/////////////////////  MAIN  ////////////////////
int funcGtot_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS,
	                 int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
	
  double Ctriple=1e-10; // Valeur de C sur le point triple (sinon C=0)
  
  int ind_air=0, ind_etain=1; 
	
	// variables
  double t, dt, dt_new;
  int n=0;
  int err=0;
  double rhoc, rhocL, rhocR, rhoc2L, rhoc2R;
  char* W0;
  double mu;
  double phiplusU, phimoinsU;
  double phiplusP, phimoinsP;
  double phi_lim;
  double r_plus, r_moins;
  double ms, mi;
  double sum1, sum2, sum3, sum4, sum5;
  double dX;
  int zone;
  double dUdx;
  double delta_eps;
  double* pVAR=malloc(2*sizeof(double));
  
  double t_start,t_end,t_sum=0;

  int nb_iter_point_fix=2;  // Extension faiblement non-linéaire ordre 2 du solveur de Riemann acoustique lagrangien
  
  // Extention d'ordre 2 et 3 du solveur de Riemann
  int nmax_iter;
  double pR,pL,uR,uL,rhoR,rhoL,cR,cL,gR,gL;
  double det, inc; 
  int iter;  
  double fR, fL, dfR, dfL, dU, dP;
  double a0R,a0L,a1R,a1L,a2R,a2L,a3R,a3L,Delta;
  int sum_sys_deg=0;
  double k2R,k2L,k3R,k3L;
  double a0,a1,a2,a3;

  // Déclaration pour le schéma Dukowicz
  double WL, PL, AL, SSL, RHOL;
  double WR, PR, AR, SSR, RHOR;
  double WMIN, WMAX, PLMIN, PRMIN;
  double BL, BR;
  double Aa, Bb, Cc, Dd;
  double DD, W12, P12;
  double* tabDD=malloc(4*sizeof(double));
  double* tabW12=malloc(4*sizeof(double));

  // Solveur de Riemann exacte
  double gamma_l,gamma_r;
  double pzero_l, pzero_r;
  double gp1s2_l, gp1s2_r, gm1s2_l, gm1s2_r, gp1s2g_l, gp1s2g_r, gm1s2g_l, gm1s2g_r;
  double rc_l, rc_r, r_l, r_r, p_l, p_r, u_l, u_r, c_l, c_r ;
  double a_l, a_r, f_l, f_r, fp_l, fp_r;
  double pi_l, pi_r;
  int nb_newton=5;

  // Solveur Gallice
  double Pjump, Ujump;
  double zepsL, zepsR;

  double ztau1L, ztau1R;
  double ztau2L, ztau2R;
  double ztauL, ztauR;
  double r, DeltaL, DeltaR;

  double aL, aR, bL, bR;
  double dx=(b-a)/nx; // pas d'espace

  // Calcul de la taille des tableaux
  int nbghostsPW=4; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int nbghosts=2*4; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int N=nx+2*nbghosts;
  int ideb=nbghosts;
  int ifin=nbghosts+nx-1;
  int idebPW=nbghosts - nbghostsPW;
  int ifinPW=ifin+nbghostsPW;

  // Coefficient  Reconstruction spatiale
  double** Ck=alloctabd(5,5); double** Ckbar=alloctabd(5,5);
  double** dk=alloctabd(5,5); double** Qbar=alloctabd(5,5);
  double** rk=alloctabd(5,5);
  initCoef(Ck, Ckbar, dk, Qbar, rk);


  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double));   // maillage
  double* Xc=malloc((N)*sizeof(double));    // maillage
  double* RHO0=malloc((N)*sizeof(double)); // masse donc centrées aux mailles

  // Point-wise
  double* TAU=malloc((N)*sizeof(double)); //  tau  (centrées aux mailles)
  double* U=malloc((N)*sizeof(double));   //  u (frontières)
  double* E=malloc((N)*sizeof(double));   //  energie totale (centrées aux mailles)

  //RK1 Matsuno + solveur acoustique ordre 1
  double* TAU_pr=malloc((N)*sizeof(double)); 
  double* U_pr=malloc((N)*sizeof(double));   
  double* E_pr=malloc((N)*sizeof(double));   
  double* lambda_unuse=malloc((3)*sizeof(double));   

  double* RHO=malloc((N)*sizeof(double));
  double* P=malloc((N)*sizeof(double));
  double* EPS=malloc((N)*sizeof(double));
  double* C=malloc((N)*sizeof(double));
  double* RC2=malloc((N)*sizeof(double));
  double* G=malloc((N)*sizeof(double));
  double* S=malloc((N)*sizeof(double));
  double* PU=malloc((N)*sizeof(double));

  //cintétique 
  double* TAUn=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* Un=malloc((N)*sizeof(double));   //  rho0.u (frontières)

  // Flux
  double* Ustar=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar=malloc((N+1)*sizeof(double));
  double* PUstar=malloc((N+1)*sizeof(double));

  double* Ustar_01=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar_01=malloc((N+1)*sizeof(double));

  // Diagnostic
  double* RE=malloc((N)*sizeof(double));
  double* RU=malloc((N)*sizeof(double));
  double* RR=malloc((N)*sizeof(double));
  double* RT=malloc((N)*sizeof(double));
  
  double* T=malloc((N)*sizeof(double));  // température 
  double** matFM=alloctabd(N,3);         // fraction massique

  int* tabNiter=malloc((N)*sizeof(int)); // nb iterration par mailles
  for(int i=0; i<N; ++i){tabNiter[i]=0;}

  double* tabEPSILON=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-1e-15;}

  double* test_consis_thermo=malloc((N)*sizeof(double)); // valeur du critère
  for(int i=0; i<N; ++i){tabEPSILON[i]=-99;}

  int* indFLUIDE=malloc((N)*sizeof(int)); // indice du fluide sur chaque maille
  
  double* Zmax=malloc((N)*sizeof(double)); // calcul du pas de temps gallice

  // Maillage
  for(int i=0; i<=N; i++){
	X[i]=(a-nbghosts*dx)+dx*i;
  //printf("i= %d, X=%lf\n",i,X[i]);
  }

  err=init_fluides(tst_multi, indFLUIDE, N, ind_air, ind_etain);
  if(err){printf("Erreur fluide\n"); return 1;}

  // Initialisation 
  double x;
  // variables défini aux centres des mailles
  double* Wpw=malloc(10*sizeof(double));
  double* Wav=malloc(10*sizeof(double));
  //#pragma omp parallel for
  for(int i=0; i<=N-1; i++){
    x=(X[i]+X[i+1])/2.0;

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
    err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, X[i], X[i+1], Wav);
    //if(err){return err;}
    
    RHO0[i]=Wav[0];
    TAU[i]=Wav[1];
    U[i]=Wav[2];
    E[i]=Wav[3];
    EPS[i]=Wav[4];
    
    P[i]=Wpw[2];
    T[i]=Wpw[3];  
    matFM[i][0]=Wpw[4];
    matFM[i][1]=Wpw[5];
    matFM[i][2]=Wpw[6];

    if(indFLUIDE[i]==ind_etain){
      err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], matFM[i], 
                     Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
                     nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                     SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      //if(err){printf(" initialisation \n"); return err;}
      
      G[i]=pVAR[0];
      S[i]=pVAR[1];
      //err= G_S_cin_3_phase(epsilon, P[i], TAU[i], EPS[i], matFM[i], &C[i], &G[i], &S[i]);  if(err){return err;}
    }
    else if(indFLUIDE[i]==ind_air){
      EGS(tst_air, TAU[i], EPS[i], &P[i], &C[i]);
      G_EGS(tst_air, TAU[i], EPS[i], &G[i]);
      //printf("Init  i= %d | G=%g\n",i,G[i] );
    }

    if(C[i]==0){C[i]=Ctriple;  }
    //printf("i=%d : ind_fluide= %d, TAU= %g, U= %g, E= %g, EPS=%g, P= %g, RHO= %g, RC2=%g, C=%g, G=%g\n",i,indFLUIDE[i],TAU[i],U[i],E[i],EPS[i],P[i],RHO0[i],RC2[i],C[i],G[i] );
  }
  free(Wpw);
  free(Wav);

  /*
  for(int i=0; i<=N-1; i++){
    printf("i=%d | TAU=%g, U=%g, E=%g, C=%g, P=%g G=%g, lambda=[%g %g %g]\n",i,TAU[i],U[i],E[i],C[i],P[i],G[i],matFM[i][0],matFM[i][1],matFM[i][2] );
  }
  */

  //for(int i=0; i<=N-1; i++){ G[i]=2.0;}

  // Calcul du 1er pas de temps
  dt=(X[ideb+1]-X[ideb])/C[ideb];
  for(int j=ideb+1; j<=ifin; j++){
    dt = fmin(dt, (X[j+1]-X[j])/C[j] );
  }
  dt = dt_first* CFL * dt; // pas de temps
  printf("nx= %d et dt= %g\n", nx, dt);

  // energie totale du système
  int nbiter=10*ceil(Tf/dt); 
  if(nbiter<=0){printf("Erreur nbiter= %d\n",nbiter ); return 1;}
  //printf("nbiter= %d\n", nbiter);
  double* Etot=malloc(nbiter*sizeof(double));
  double* IMPUL=malloc(nbiter*sizeof(double)); // impulsion
  double* Mtot=malloc(nbiter*sizeof(double));  // masse totale
  double* VOL=malloc(nbiter*sizeof(double));   // volume total

  sum1=0;  sum3=0.; sum4=0.; sum5=0.;
  #pragma omp parallel for
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*RHO0[i]*E[i];
    sum3+=dx*RHO0[i]*U[i];
    sum4+=(X[i+1]-X[i])/TAU[i];
    sum5+=dx*RHO0[i]*TAU[i];
  }
  Etot[0]=sum1;
  IMPUL[0]=sum3; 
  Mtot[0]=sum4;
  VOL[0]=sum5;

  int solveur; // solevur d'ordre 1
  int methd;  // méthode de montée à l'ordre 2  
  // 0 si la méthode est d'ordre 1
  int limiteur; // limiteur utilisé dans la montée en ordre
  err=choix_sch(sch, &solveur, &methd, &limiteur);



  // debut de la boucle ============================================================================= 
  while((t<Tf)&&(n<Nmax)){
    // Calcul du flux 

    if(solveur==0){ // Despres
	    #pragma omp parallel for 
      for(int i=ideb; i<=ifin+1; i++){
        rhoc=(C[i-1]/TAU[i-1]+ C[i]/TAU[i])/2; 
        Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
        Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("  i=%d, Ustar=%g, Pstar=%g, rhoc=%g, TAU= %g , C= %g\n",i, Ustar[i],Pstar[i],rhoc,TAU[i],C[i]);
      }
    }
    else if(solveur==1){ // Jaouen 2001
	    #pragma omp parallel for 
      for(int i=ideb; i<=ifin+1; i++){
        rhoc2L=C[i-1]*C[i-1]/TAU[i-1]; 
        rhoc2R=C[i]*C[i]/TAU[i];
        ms=fmax(rhoc2R,rhocL);
        mi=fmax(TAU[i],TAU[i-1]);
        rhoc=sqrt(ms/mi);
        Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
        Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
      }
    }
    else if(solveur==2){ // Solveur Acoustique aux faces ordre 1 
      #pragma omp parallel for 
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU[i-1]; 
        rhocR=C[i]/TAU[i];
        Ustar[i]=(rhocL*U[i-1]+rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR) ;
        Pstar[i]=(rhocR*P[i-1]+rhocL*P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
      }
    }
    /*
    else if(solveur==3){ // Extension faiblement non-linéaire à l'ordre 2
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU[i-1];
        rhocR=C[i]/TAU[i];
        Ustar[i]=(rhocL*U[i-1]+rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR);
        Pstar[i]=(rhocR*P[i-1]+rhocL+P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]);
        PUstar[i]=Pstar[i]*Ustar[i];
      }
      for(int k=0; k<nb_iter_point_fix; k++){
        
        for(int i=ideb; i<=ifin+1; i++){
          //rhocL=(1./TAU[i-1])*(C[i-1] - 0.5*G[i-1]*(Ustar[i]-U[i-1]));
          //rhocR=(1./TAU[i])*(C[i] + 0.5*G[i]*(Ustar[i]-U[i]));
          
          rhocL = rhocL + 1.00*((1./TAU[i-1])*(C[i-1] - 0.5*G[i-1]*(Ustar[i]-U[i-1]))-rhocL);
          rhocR = rhocR + 1.00*((1./TAU[i])*(C[i] + 0.5*G[i]*(Ustar[i]-U[i]))-rhocR);
          
          Ustar[i]=(rhocL*U[i-1] + rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR);
          Pstar[i]=(rhocR*P[i-1] + rhocL*P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]);
          PUstar[i]=Pstar[i]*Ustar[i];
          
          //printf("i=%d | Ustar=%g Pstar=%g rhocL=%g rhocR=%g G=%g C=%g \n",i,Ustar[i],Pstar[i],rhocL,rhocR,G[i],C[i] );
        }
      }
    }
    */
    else if(solveur==3){ // Extension faiblement non-linéaire à l'ordre 2
	    #pragma omp parallel for 
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU[i-1];
        rhocR=C[i]/TAU[i];
        Ustar[i]=(rhocL*U[i-1]+rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR);
        Pstar[i]=(rhocR*P[i-1]+rhocL+P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]);
        PUstar[i]=Pstar[i]*Ustar[i];
        
        pR=P[i]; pL=P[i-1];
        uR=U[i]; uL=U[i-1];
        rhoR=1./TAU[i]; rhoL=1./TAU[i-1];
        cR=C[i]; cL=C[i-1];
        gR=G[i]; gL=G[i-1];

        // Test de l'existance et uncite de la solution
        a0R=pR - rhoR*cR*uR + (1./2)*rhoR*gR*uR*uR;
        a0L=pL + rhoL*cL*uL + (1./2)*rhoL*gL*uL*uL;

        a1R= rhoR*cR - rhoR*gR*uR;
        a1L= -rhoL*cL - rhoL*gL*uL;

        a2R= (1./2)*rhoR*gR;
        a2L= (1./2)*rhoL*gL;

        Delta= (a1R - a1L)*(a1R - a1L) - 4*(a2R - a2L)*(a0R - a0L);
        //printf("i=%d | delta=%g",i,Delta );
        if(Delta<0){
          sum_sys_deg+=1;
          printf("     nombre de mailles ou le système n'est pas résoluble : %d\n",sum_sys_deg );
        }

        // Méthode de Newton-Raphson
        iter=0;
        inc=1;
        nmax_iter=30;
        while(inc>1e-8 && iter<nmax_iter){
          fR= - (1./2)*rhoR*gR*(Ustar[i]-uR);
          fR= (fR - rhoR*cR)*(Ustar[i]-uR); 
          fR+= Pstar[i] - pR;
          
          fL= - (1./2)*rhoL*gL*(Ustar[i]-uL);
          fL= (fL + rhoL*cL)*(Ustar[i]-uL);
          fL+= Pstar[i] - pL;
          
          dfR=- rhoR*gR*(Ustar[i]-uR);
          dfR+= -rhoR*cR;

          dfL= - rhoL*gL*(Ustar[i]-uL);
          dfL+= rhoL*cL;

          det=dfL - dfR;
          //if(fabs(det)<1e-10){printf("Erreur Jacobienne non-inversible : Extension faiblement non-linéaire à l'ordre 3 (funcGtot_multi) det=%g\n",det); return 1;}

          dU= (-1./det)*(fL - fR);
          dP= (-1./det)*(fR*dfL - fL*dfR);

          Ustar[i] += dU;
          Pstar[i] += dP;
          PUstar[i]=Pstar[i]*Ustar[i];
          
          if(Ustar[i]==0){ 
            if(Pstar[i]==0){ inc=dP+dU; }
            else{ inc=fabs(dP/Pstar[i]); } 
          }
          else{
           if(Pstar[i]==0){ inc=fabs(dU/Ustar[i]); }
           else{ inc=fabs(dU/Ustar[i] + dP/Pstar[i]); }
          }
          iter++;

          //printf("    det=%g fR=%g fL=%g dfR=%g dfL=%g\n",det,fR,fL,dfR,dfL);
        }
        //if(iter==nmax_iter){printf("Erreur Newton convergent (%d iterations) n=%d (funcGtot_multi)\n",nmax_iter,n); return 1;}
        //printf("i=%d | Newton-Raphson itération= %d  dP/P=%g dU/U=%g inc=%g \n",i,iter,dP/Pstar[i],dU/Ustar[i],inc);
      }
    }
    else if(solveur==4){ // Extension faiblement non-linéaire à l'ordre 3
  	  #pragma omp parallel for
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU[i-1];
        rhocR=C[i]/TAU[i];
        Ustar[i]=(rhocL*U[i-1]+rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR);
        Pstar[i]=(rhocR*P[i-1]+rhocL+P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]);
        PUstar[i]=Pstar[i]*Ustar[i];
        
        pR=P[i]; pL=P[i-1];
        uR=U[i]; uL=U[i-1];
        rhoR=1./TAU[i]; rhoL=1./TAU[i-1];
        cR=C[i]; cL=C[i-1];
        gR=G[i]; gL=G[i-1];

        // Test de l'existance et uncite de la solution
        k2R=(1./2)*rhoR*gR;
        k2L=(1./2)*rhoL*gL;

        k3R= (1./(6*cR))*(gR+2)*(2*gR+5);
        k3L=-(1./(6*cL))*(gL+2)*(2*gL+5);

        a0R=pR - rhoR*cR*uR + k2R*(uR*uR) + k3R*(uR*uR*uR);
        a0L=pL + rhoL*cL*uL + k2L*(uL*uL) + k3L*(uL*uL*uL);
        a0=a0R-a0L;

        a1R=   rhoR*cR*uR - 2*k2R*uR + 3*k3R*(uR*uR);
        a1L= - rhoL*cL*uL - 2*k2L*uL + 3*k3L*(uL*uL);
        //a1R= rhoR*cR - rhoR*gR*uR;
        //a1L= -rhoL*cL - rhoL*gL*uL;
        a1=a1R-a1L;
        
        a2R= k2R*uR - 3*k3R*uR;
        a2L= k2L*uL - 3*k3L*uL;
        //a2R= (1./2)*rhoR*gR;
        //a2L= (1./2)*rhoL*gL;
        a2=a2R-a2L;

        a3R= k3R;
        a3L= k3L;
        a3=a3R-a3L;

        Delta = (a2*a2)*(a1*a1) + 18*a3*a2*a1*a0 - 27*(a3*a3)*(a0*a0) - 4*a3*(a1*a1*a1) - 4*(a2*a2*a2)*a0;
        //printf("n=%d i=%d | delta=%g",n,i,Delta );
        if(Delta>0){
          sum_sys_deg+=1;
          printf("     nombre de mailles ou le système n'a pas une unique solution : %d\n",sum_sys_deg );
        }


        // Méthode de Newton-Raphson
        iter=0;
        inc=1;
        nmax_iter=500;
        while(inc>1e-5 && iter<nmax_iter){
          fR= -(1./(6*cR))*(gR+2)*(2*gR+5)*(Ustar[i]-uR);
          fR= (fR - 0.5*rhoR*gR)*(Ustar[i]-uR);
          fR= (fR - rhoR*cR)*(Ustar[i]-uR) + Pstar[i] - pR;
          
          fL= (1./(6*cL))*(gL+2)*(2*gL+5)*(Ustar[i]-uL);
          fL= (fL - 0.5*rhoL*gL)*(Ustar[i]-uL);
          fL= (fL + rhoL*cL)*(Ustar[i]-uL) + Pstar[i] - pL;
          
          dfR= -(1./(2*cR))*(gR+2)*(2*gR+5)*(Ustar[i]-uR);
          dfR= (dfR - rhoR*gR)*(Ustar[i]-uR);
          dfR+= -rhoR*cR;

          dfL= (1./(2*cL))*(gL+2)*(2*gL+5)*(Ustar[i]-uL);
          dfL= (dfL - rhoL*gL)*(Ustar[i]-uL);
          dfL+= rhoL*cL;

          det=dfL - dfR;
          //if(fabs(det)<1e-10){printf("Erreur Jacobienne non-inversible : Extension faiblement non-linéaire à l'ordre 3 (funcGtot_multi) det=%g\n",det); return 1;}

          dU= (-1./det)*(fL - fR);
          dP= (-1./det)*(fR*dfL - fL*dfR);

          Ustar[i] += 0.2*dU;
          Pstar[i] += 0.2*dP;
          PUstar[i]=Pstar[i]*Ustar[i];
          
          if(Ustar[i]==0){ 
            if(Pstar[i]==0){ inc=dP+dU; }
            else{ inc=fabs(dP/Pstar[i]); } 
          }
          else{
           if(Pstar[i]==0){ inc=fabs(dU/Ustar[i]); }
           else{ inc=fabs(dU/Ustar[i] + dP/Pstar[i]); }
          }
          iter++;

          //printf("    det=%g fR=%g fL=%g dfR=%g dfL=%g\n",det,fR,fL,dfR,dfL);
          //printf("    det=%g inc=%g fR=%g fL=%g\n",det,inc,fR, fL);
        }
        //if(iter==nmax_iter){printf("Erreur Newton convergent (%d iterations) n=%d (funcGtot_multi)\n",nmax_iter,n); return 1;}
        //printf("i=%d | Newton-Raphson itération= %d  dP/P=%g dU/U=%g inc=%g \n",i,iter,dP/Pstar[i],dU/Ustar[i],inc);
        
        // test de solution par le newton
        /*
        printf("i=%d\n",i );
        printf("  p*-PR(u*)=%g\n", (Pstar[i]-(((a3R*Ustar[i]+a2R)*Ustar[i] + a1R)*Ustar[i] + a0R))/Pstar[i] );
        printf("  p*-PL(u*)=%g\n", (Pstar[i]-(((a3L*Ustar[i]+a2L)*Ustar[i] + a1L)*Ustar[i] + a0L))/Pstar[i] );
        */
      }
    }
    else if(solveur==5){ // Solveur de Dukovicz
      //    Dukowicz's formulae BUT with "as" and "As" assigned to the EOS-provided sound speed and *half* fundamental derivative
      //    John K. Dukovicz, "A General, Non-Iterative Riemann Solver for Godunov's Method", Journal of Computational Physics 61, 119-137 (1985)
      //    Dukovicz's formula (24) involves *expensive* EOS calls per cell per cycle with Rankine-Hugoniot *infinite strengh* shock computations
      //    Here, the two-shock approximation is retained but in the *weak shock limit* with sound speed & fundamental derivative computed from the EOS
       //    For a perfect gas, Dukowicz's formula (24) is very close to the fundamental derivative, different from *half* the fundamental derivative here
    
       //      LEfT AS EXERCISE: extension/modification of the Dukowicz's formulae for NEGATIVE values of fundamental derivative (e.g. Bizarrium EOS)
      #pragma omp parallel for
      for(int i=ideb; i<=ifin+1; i++){

        WL = U[i-1];
        PL = P[i-1];
        AL = G[i-1]*0.5;    // Half the fundamental derivative -> our choice to satisfy the WEAK SHOCK limit
        SSL = C[i-1];
        RHOL = 1./TAU[i-1];

        WR = U[i];
        PR = P[i];
        AR = G[i]*0.5;
        SSR = C[i];
        RHOR = 1./TAU[i];

        if(sch==7){      //  Dukowicz's STRONG SHOCK limit formula (24) in the case of a perfect gas EOS
          AL = 2*AL;     //  See also Juan Cheng and Chi-Wang Shu, Journal of Computational Physics 227, 1567-1596 (2007)
          AR = 2*AR;     //  A factor 2 (two) has to be applied to the coefficients...
        }

        WMIN = WR - 0.5*SSR/AR;
        WMAX = WL + 0.5*SSL/AL;

        PLMIN = PL - 0.25*RHOL*SSL*SSL/AL;
        PRMIN = PR - 0.25*RHOR*SSR*SSR/AR;

        BL = RHOL*AL;
        BR = RHOR*AR;

        Aa = (BR-BL)*(PRMIN-PLMIN);
        Bb = BR*WMIN*WMIN - BL*WMAX*WMAX;
        Cc = BR*WMIN-BL*WMAX;
        Dd = BR*BL*(WMIN-WMAX)*(WMIN-WMAX);

        DD = sqrt(fmax(0, Dd-Aa));
        tabW12[0] = (Bb+PRMIN-PLMIN)/(Cc-fabs(DD)*sgn(WMAX-WMIN));
        DD = sqrt(fmax(0, Dd+Aa));
        tabW12[1] = (Bb-PRMIN+PLMIN)/(Cc-fabs(DD)*sgn(WMAX-WMIN));
        
        Aa = (BL+BR)*(PLMIN-PRMIN);
        Bb = BL*WMAX+BR*WMIN;
        Cc = 1./(BL+BR);
        DD = sqrt(fmax(0, Aa-Dd));
        tabW12[2] = (Bb+DD)*Cc;
        DD = sqrt(fmax(0, -Aa-Dd));
        tabW12[3] = (Bb-DD)*Cc;

        //printf("i=%d WMIN=%g WMAX=%g | ",i,WMIN,WMAX );
        //printf("W12=%g  W12=%g  W12=%g  W12=%g | ",tabW12[0],tabW12[1],tabW12[2],tabW12[3] );

        //  Case A: W12-WMIN > 0 and W12-WMAX < 0
        if (tabW12[0]-WMIN >= 0. && tabW12[0]-WMAX <= 0.){
          W12=tabW12[0];
          //printf(" case A\n");
        }
        //  Case B: W12-WMIN < 0 and W12-WMAX > 0
        else if (tabW12[1]-WMIN <=0. && tabW12[1]-WMAX >=0.){
          W12=tabW12[1];
          //printf(" case B\n" );
        }
        //  Case C: W12-WMIN > 0 and W12-WMAX > 0
        else if (tabW12[2]-WMIN >= 0. && tabW12[2]-WMAX >=0.){
          W12=tabW12[2];
          //printf(" case C \n");
        }
        //  Case D: W12-WMIN < 0 and W12-WMAX < 0
        else{ //} if (tabW12[3]-WMIN <= 0. && tabW12[3]-WMAX <=0.){
          W12=tabW12[3];
          //printf(" case D\n");
        }
        //else{printf(" erreur Dukowicz case  i= %d\n",i );}

        P12 = 0.5*(PLMIN+PRMIN+BR*fabs(W12-WMIN)*(W12-WMIN)-BL*fabs(W12-WMAX)*(W12-WMAX));

        Ustar[i] = W12 ;
        Pstar[i] = P12 ;
        PUstar[i]=W12*P12;

      }
    }
    else if(solveur==6){ // Solveur de Riemann exacte
      #pragma omp parallel for
      for(int i=ideb; i<=ifin+1; i++){
        //printf(" maille %d\n",i );

        // Left and right Stiffened Gas coefficients computed from the EOS-provided fundamental derivative and sound velocity
        gamma_l = 2.0*G[i-1] -1;
        gamma_r = 2.0*G[i] -1;

        pzero_l = (1./gamma_l)*(1./TAU[i-1])*C[i-1]*C[i-1] - P[i-1];
        pzero_r = (1./gamma_r)*(1./TAU[i])*C[i]*C[i] - P[i];

        // Newton's iterative procedure from Godounov et coll. pp. 114-117 (Mir, 1979)
        gp1s2_l = (gamma_l+1)/2;
        gp1s2_r = (gamma_r+1)/2;

        gm1s2_l = (gamma_l-1)/2;
        gm1s2_r = (gamma_r-1)/2;

        gp1s2g_l = (gamma_l+1)/(2*gamma_l);
        gp1s2g_r = (gamma_r+1)/(2*gamma_r);

        gm1s2g_l = (gamma_l-1)/(2*gamma_l);
        gm1s2g_r = (gamma_r-1)/(2*gamma_r);

        rc_l = C[i-1]/TAU[i-1];
        rc_r = C[i]/TAU[i];

        r_l = 1./TAU[i-1];  r_r = 1./TAU[i];
        p_l = P[i-1]; p_r = P[i];
        u_l = U[i-1]; u_r = U[i];
        c_l = C[i-1]; c_r = C[i];

        //# Initial guess using the *acoustic* formulae
        
        Ustar[i] = ( rc_l*U[i-1] + rc_r*U[i] + (P[i-1] - P[i]) ) / (rc_l + rc_r);
        Pstar[i] = ( rc_r*P[i-1] + rc_l*P[i] + rc_l*rc_r*(U[i-1] - U[i]) ) / (rc_l + rc_r);

        for(int k=0; k<nb_newton; k++){
          //printf("  iteration Newton= %d/%d |  ",k,nb_newton );
          pi_l = (Pstar[i]+pzero_l)/(p_l+pzero_l);
          pi_r = (Pstar[i]+pzero_r)/(p_r+pzero_r);

          if(sch==8){ //riemann == "two-shock_Godunov" #  Shock formulae only for *cheaper* approximate Riemann solver
            a_l = sqrt(r_l*(gp1s2_l*(Pstar[i]+pzero_l) + gm1s2_l*(p_l+pzero_l)));
            a_r = sqrt(r_r*(gp1s2_r*(Pstar[i]+pzero_r) + gm1s2_r*(p_r+pzero_r)));
            f_l = (Pstar[i]-p_l)/(rc_l*sqrt(gp1s2g_l*pi_l+gm1s2g_l));
            f_r = (Pstar[i]-p_r)/(rc_r*sqrt(gp1s2g_r*pi_r+gm1s2g_r));
            fp_l = ((gamma_l+1)*pi_l+(3*gamma_l-1))/(4*gamma_l*rc_l*sqrt(pow(gp1s2g_l*pi_l+gm1s2g_l,3)));
            fp_r = ((gamma_r+1)*pi_r+(3*gamma_r-1))/(4*gamma_r*rc_r*sqrt(pow(gp1s2g_r*pi_r+gm1s2g_r,3)));
          }
          else if(sch==9){ // riemann == "two-rarefaction_Godunov" #  Rarefaction formulae only

            if(fabs(Pstar[i]-p_l) <= 1.e-11*p_l){ //  # Acoustic formulae to prevent NaN between cells at equilibrium
              a_l = rc_l;
              f_l = (Pstar[i]-p_l)/rc_l;
              fp_l = 1./rc_l;
            }
            else{
              if(fabs(pi_l-1.0)<=1e-13){a_l=0;}
              else{a_l = gm1s2g_l*rc_l*(1-pi_l)/(1-pow(pi_l,gm1s2g_l));}
              f_l = 2./(gamma_l-1)*c_l*(pow(pi_l,gm1s2g_l) - 1);
              fp_l = 1./(gamma_l*(p_l+pzero_l))*c_l*pow(pi_l,gm1s2g_l);
            }

            if( fabs(Pstar[i]-p_r) <= 1.e-11*p_r){ //  # Acoustic formulae to prevent NaN between cells at equilibrium
              a_r = rc_r;
              f_r = (Pstar[i]-p_r)/rc_r;
              fp_r = 1./rc_r;
            }
            else{
              if(fabs(pi_r-1.0)<=1e-13){a_r=0;}
              else{a_r = gm1s2g_r*rc_r*(1-pi_r)/(1-pow(pi_r,gm1s2g_r));}
              f_r = 2./(gamma_r-1)*c_r*(pow(pi_r,gm1s2g_r) - 1);
              fp_r = 1./(gamma_r*(p_r+pzero_r))*c_r*pow(pi_r,gm1s2g_r);
            }
          }
          else if(sch==10){ // riemann == "exact_Godunov" # Standard formulae for the EXACT Riemann solver

            if( (Pstar[i]-p_l) >= - 1.e-11*fabs(p_l)){ //  # Shock formulae To prevent NaN between cells at equilibrium
              a_l = sqrt(r_l*(gp1s2_l*(Pstar[i]+pzero_l) + gm1s2_l*(p_l+pzero_l)));
              f_l = (Pstar[i]-p_l)/(rc_l*sqrt(gp1s2g_l*pi_l+gm1s2g_l));
              fp_l = ((gamma_l+1)*pi_l+(3*gamma_l-1))/(4*gamma_l*rc_l*sqrt(pow(gp1s2g_l*pi_l+gm1s2g_l,3)));
            }
            else{
              if(fabs(pi_l-1.0)<=1e-13){a_l=0;}
              else{a_l = gm1s2g_l*rc_l*(1-pi_l)/(1-pow(pi_l,gm1s2g_l));}
              f_l = 2./(gamma_l-1)*c_l*(pow(pi_l,gm1s2g_l) - 1);
              fp_l = 1./(gamma_l*(p_l+pzero_l))*c_l*pow(pi_l,gm1s2g_l);
            }

            if( (Pstar[i]-p_r) >= -1.e-11*fabs(p_r)){
              a_r = sqrt(r_r*(gp1s2_r*(Pstar[i]+pzero_r) + gm1s2_r*(p_r+pzero_r)));
              f_r = (Pstar[i]-p_r)/(rc_r*sqrt(gp1s2g_r*pi_r+gm1s2g_r));
              fp_r = ((gamma_r+1)*pi_r+(3*gamma_r-1))/(4*gamma_r*rc_r*sqrt(pow(gp1s2g_r*pi_r+gm1s2g_r,3)));
            
              //printf("\n           a_r = sqrt(r_r*(gp1s2_r*(Pstar[i]+pzero_r) + gm1s2_r*(p_r+pzero_r)))\n");
              //printf("          %g = sqrt(%g*(%g+%g)) = sqrt(%g*%g)\n",a_r, r_r,gp1s2_r*(Pstar[i]+pzero_r), gm1s2_r*(p_r+pzero_r),r_r,gp1s2_r*(Pstar[i]+pzero_r)+ gm1s2_r*(p_r+pzero_r));

            }
            else{
              if(fabs(pi_r-1.0)<=1e-13){a_r=0;}// printf("         a_r=0 car pi_r=%.14g",pi_r); }
              else{a_r = gm1s2g_r*rc_r*(1-pi_r)/(1-pow(pi_r,gm1s2g_r));}
              f_r = 2./(gamma_r-1)*c_r*(pow(pi_r,gm1s2g_r) - 1);
              fp_r = 1./(gamma_r*(p_r+pzero_r))*c_r*pow(pi_r,gm1s2g_r);

              //printf("\n          a_r = gm1s2g_r*rc_r*(1-pi_r)/(1-pow(pi_r,gm1s2g_r)) \n");
              //printf("          %g = %g*%g*(1-%g)/(1-pow(%g,%g)) \n",a_r, gm1s2g_r,rc_r,pi_r,pi_r,gm1s2g_r);
            }
            
            /*
            printf("\n          Newton=%d | testl=(p* - p_l)/p_l= %.14g, testr=(p* - p_r)/p_r= %.14g \n",k,(Pstar[i] - p_l)/p_l, (Pstar[i] - p_r)/p_r );
            printf("          a_l=%g f_l=%g fp_l=%g \n",a_l, f_l, fp_l);
            printf("          a_r=%g f_r=%g fp_r=%g \n",a_r, f_r, fp_r);
            
            printf("\n          pi_l=(p*+pzero_l)/(p_l+pzero_l)\n" );
            printf("          %.14g = (%g+%g)/(%g+%g)\n",pi_l,Pstar[i],pzero_l,p_l,pzero_l);

            printf("\n          pi_r=(p*+pzero_r)/(p_r+pzero_r)\n");
            printf("          %.14g = (%g+%g)/(%g+%g)\n",pi_r,Pstar[i],pzero_r,p_r,pzero_r);


            printf("\n          gamma_l=%g pzero_l=%g G_l=%g \n",gamma_l, pzero_l,G[i-1]);
            printf("          gamma_r=%g pzero_r=%g G_r=%g \n",gamma_r, pzero_r, G[i]);
            
            printf("\n          gp1s2_l=%g gm1s2_l=%g gp1s2g_l=%g gm1s2g_l=%g  \n",gp1s2_l, gm1s2_l, gp1s2g_l, gm1s2g_l);
            printf("          gp1s2_r=%g gm1s2_r=%g gp1s2g_r=%g gm1s2g_r=%g  \n",gp1s2_r, gm1s2_r, gp1s2g_r, gm1s2g_r);
            */
          }

          Pstar[i] -= (f_l+f_r-(u_l-u_r))/(fp_l+fp_r);
          Ustar[i] = (a_l*u_l+a_r*u_r+p_l-p_r)/(a_l+a_r);
          PUstar[i] = Pstar[i]*Ustar[i];            

          //printf("       p*= %.14g, u*= %.14g\n",Pstar[i],Ustar[i] );
          //printf("             ----------------------------\n");
          
        }
        //printf("i=%d | p*= %.14g, u*= %.14g\n",i,Pstar[i],Ustar[i] );
      }
    }
    else if(solveur==7){ // Solveur de Gallice 
      /* Sur le cas-test de LeBlanc : acoustic et Gallice 2 mais pas Gallice 1  (Robustesse)
         Sur les cas-test de Sod et Bizarruim aucune diff entre les 3
      */
      #pragma omp parallel for
      for(int i=ideb; i<=ifin+1; i++){
        Pjump=P[i]-P[i-1]; 
        Ujump=U[i]-U[i-1];

        // cas avec convexité des isentropes
        
        zepsL= C[i-1]/TAU[i-1];
        zepsR= C[i]/TAU[i];
        
        // cas sans convexité des isentropes
        //   conforme au papier
        /*
        zepsL= P[i-1]/sqrt(2*EPS[i-1]);
        zepsR= P[i]/sqrt(2*EPS[i]);
        */
        // mon calcul
        /*
        zepsL= P[i-1]*P[i-1]/(2*EPS[i-1]);
        zepsR= P[i]*P[i]/(2*EPS[i]);
        */

        if(sch==11){ // Méthode 1
          ztau1L=sqrt( fmax(Pjump,0.)/TAU[i-1] );
          ztau1R=sqrt( fmax(-Pjump,0.)/TAU[i] );
          
          ztau2L=-Ujump/TAU[i-1];
          ztau2R=-Ujump/TAU[i];

          ztauL=fmax(ztau1L,ztau2L);
          ztauR=fmax(ztau1R,ztau2R);
        }
        else if(sch==12){ // Méthode 2
          r=TAU[i-1]*C[i]/(TAU[i]*C[i-1]);  //(rhocR/rhocL)

          DeltaL=r*r*Ujump*Ujump + 4*(1+r)*TAU[i-1]*Pjump;
          DeltaR=Ujump*Ujump - 4*r*(1+r)*TAU[i]*Pjump;

          if(DeltaL<=0 && DeltaR<=0){
            ztauL=zepsL; 
            ztauR=zepsR;
          }
          else if(DeltaL<=0 && DeltaR>0){
            ztauR=( -Ujump+sqrt(DeltaR) )/( 2*(1+r)*TAU[i] );
            ztauL=ztauR/r;
          }
          else if(DeltaL>0 && DeltaR<=0){
            ztauL=( -r*Ujump+sqrt(DeltaL) )/( 2*(1+r)*TAU[i-1] );
            ztauR=r*ztauL; 
          }
          else if(DeltaL>0 && DeltaR>0){
            ztauL=( -r*Ujump+sqrt(DeltaR) )/( 2*(1+r)*TAU[i-1] );
            ztauR=( -Ujump+sqrt(DeltaR) )/( 2*(1+r)*TAU[i] );
          }
        }
        else if(sch==13){ // Méthode 3
          aL=Ujump/TAU[i-1];  aR=Ujump/TAU[i];
          bL=Pjump/TAU[i-1];  bR=Pjump/TAU[i];
          
          // aL,aR de même signe et idem pour bL,bR  (car TAU>0)
          // on utilise donc uniquement aL et bL pour les tests
          if(aL==0 || bL==0){
            ztauL=zepsL; 
            ztauR=zepsR;
          }
          else if(aL>0 && bL>0){  // ++
            ztauL=sqrt(bL);
            ztauR=bL/aL;
          }
          else if(aL>0 && bL<0){  // +-
            ztauL=-bL/aL;
            ztauR=sqrt(-bR);
          }
          else if(aL<0 && bL>0){  // -+
            ztauL=fmax(sqrt(bL),-aL);
            ztauR=-aR;
          }
          else if(aL<0 && bL<0){  // --
            ztauL=-aL;
            ztauR=sqrt(bR);
          }
          //else{ printf(" erreur choix condition sur tau (Gallice 3)\n  aL=%g bL=%g  (aR=%g  bR=%g) \n", aL,bL,aR,bR); return 1; }
        }

        rhocL=fmax(zepsL, ztauL); 
        rhocR=fmax(zepsR, ztauR); 

        // Calcul du pas de temps
        Zmax[i]=fmax(rhocL,rhocR);
        
        Ustar[i]=(rhocL*U[i-1]+rhocR*U[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR) ;
        Pstar[i]=(rhocR*P[i-1]+rhocL*P[i])/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d | Ustar=%g Pstar=%g z=max(zeps, ztau)  zL=%g=max(%g, %g)  zR=%g=max(%g, %g) \n",i, Ustar[i],Pstar[i],rhocL,zepsL,ztauL,rhocR,zepsR,ztauR);
      }
    }
    else{printf("Erreur dans le choix du solveur =%d (funcGtot_multi)\n",solveur); return 1;}

    if(methd==0){}
    else if(methd==1){ // GAD d'ordre 2
      delta_eps=1e-30;
      #pragma omp parallel for
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU[i-1];
        rhocR=C[i]/TAU[i];
        
        mu= (rhocL+rhocR)*dt/(dx*(RHO0[i-1]+RHO0[i]));
      
        if(limiteur==0){  // ordre 2 sans limiteur
          phimoinsU=1.0;  phimoinsP=1.0;
          phiplusU=1.0;   phiplusP=1.0;
        }
        else if(limiteur==1){ // ordre 2 avec limiteur minmod
          phimoinsU=MinMod((Ustar[i+1]-U[i]+delta_eps)/(Ustar[i]-U[i-1]+delta_eps));
          phimoinsP=MinMod((Pstar[i+1]-P[i]+delta_eps)/(Pstar[i]-P[i-1]+delta_eps));

          phiplusU=MinMod((U[i-1]-Ustar[i-1]+delta_eps)/(U[i]-Ustar[i]+delta_eps));
          phiplusP=MinMod((P[i-1]-Pstar[i-1]+delta_eps)/(P[i]-Pstar[i]+delta_eps));
        }
        else if(limiteur==2){ // ordre 2 avec limiteur superbee
          phimoinsU=Superbee((Ustar[i+1]-U[i]+delta_eps)/(Ustar[i]-U[i-1]+delta_eps));
          phimoinsP=Superbee((Pstar[i+1]-P[i]+delta_eps)/(Pstar[i]-P[i-1]+delta_eps));

          phiplusU=Superbee((U[i-1]-Ustar[i-1]+delta_eps)/(U[i]-Ustar[i]+delta_eps));
          phiplusP=Superbee((P[i-1]-Pstar[i-1]+delta_eps)/(P[i]-Pstar[i]+delta_eps));
        }
        
        Ustar[i]+=(0.5)*(1.0-mu)*( phiplusU*(U[i]-Ustar[i]) - phimoinsU*(Ustar[i]-U[i-1]) );
        Pstar[i]+=(0.5)*(1.0-mu)*( phiplusP*(P[i]-Pstar[i]) - phimoinsP*(Pstar[i]-P[i-1]) );
        PUstar[i]=Pstar[i]*Ustar[i];
      }
    }
    else{printf("Erreur dans le choix de la méthode d'ordre 2   methd=%d (funcGtot_multi)\n",methd); return 1;}


    if(sch==3){  // RK1 MAtsuno + solveur acoustique
      // Mise à jour de la sol U
      for(int i=ideb; i<=ifin; i++){
        TAU_pr[i]=TAU[i] + (dt/(dx*RHO0[i]))*(Ustar[i+1]-Ustar[i]);
        U_pr[i]=U[i] - (dt/(dx*RHO0[i]))*(Pstar[i+1]-Pstar[i]);
        E_pr[i]=E[i] - (dt/(dx*RHO0[i]))*(PUstar[i+1]-PUstar[i]);
        
        //printf("i= %d : TAU=%g, U=%g, E=%g \n",i,TAU[i],U[i],E[i] );
      }

      // cond aux limites sur les valeurs moyennes
      err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, TAU_pr, 0);
      if(err){return err;}
      err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, E_pr, 0);
      if(err){return err;}
      err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, U_pr, iu);
      if(err){return err;}
      err=condlim(ind_cond, ideb,ifin+1, nbghosts, Ustar, iu);
      if(err){return err;}

      // Mise à jour des tableaux
      //for(int i=0; i<=N; i++){
      for(int i=0; i<=N-1; i++){
        dUdx=(Ustar[i+1]-Ustar[i])/dx;
        EPS[i]= E_pr[i]-U_pr[i]*U_pr[i]/2;  // energie interne

        if(indFLUIDE[i]==ind_etain){
          err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU_pr[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], lambda_unuse, 
                       Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
          if(err){printf(" cycle=%d (RK1 Matsuno)\n",n); return err;}
        }
        else if(indFLUIDE[i]==ind_air){
          EGS_cin_phase(tst_air,  TAU_pr[i], EPS[i], &P[i], &C[i],
                        cin_phase, dx, dt, Ttau, 
                        &Un[i], P[i], TAUn[i], C[i]);
        }
        if(C[i]==0){ C[i]=Ctriple; }
      }
      
      #pragma omp parallel for 
      for(int i=ideb; i<=ifin+1; i++){
        rhocL=C[i-1]/TAU_pr[i-1]; 
        rhocR=C[i]/TAU_pr[i];
        Ustar[i]=(rhocL*U_pr[i-1]+rhocR*U_pr[i])/(rhocL+rhocR) + (P[i-1]-P[i])/(rhocL+rhocR) ;
        Pstar[i]=(P[i]*rhocL+P[i-1]*rhocR)/(rhocL+rhocR) + (rhocL*rhocR/(rhocL+rhocR))*(U_pr[i-1]-U_pr[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
      }
    }


    // Copie de TAU et U pour la cintétique (Utile pour RK1 Matsuno)
    #pragma omp parallel for 
    for(int i=0; i<=N-1; i++){  TAUn[i]=TAU[i];  Un[i]=U[i];  }


	  // Mise à jour de la sol U
    #pragma omp parallel for 
	  for(int i=ideb; i<=ifin; i++){
	    TAU[i]+=(dt/(dx*RHO0[i]))*(Ustar[i+1]-Ustar[i]);
	    U[i]-=(dt/(dx*RHO0[i]))*(Pstar[i+1]-Pstar[i]);
	    E[i]-=(dt/(dx*RHO0[i]))*(PUstar[i+1]-PUstar[i]);

      //printf("i= %d : TAU=%g, U=%g, E=%g \n",i,TAU[i],U[i],E[i] );
	  }
    
	  // Mise à jour du maillage
    #pragma omp parallel for 
	  for(int i=ideb; i<=ifin+1; i++){
	  	X[i]+=Ustar[i]*dt;
	  }

	// cond aux limites sur les valeurs moyennes
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, TAU, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, E, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, U, iu);
	  if(err){return err;}
    err=condlim(ind_cond, ideb,ifin+1, nbghosts, Ustar, iu);
    if(err){return err;}


    // Mise à jour des tableaux
    //for(int i=0; i<=N; i++){ 
    t_start = omp_get_wtime();
    #pragma omp parallel for 
    for(int i=0; i<=N-1; i++){ 
      //dUdx=(Un[i+1]-Un[i-1])/2*dx;
      dUdx=(Ustar[i+1]-Ustar[i])/dx;
      EPS[i]= E[i]-U[i]*U[i]/2;  // energie interne 

	    if(indFLUIDE[i]==ind_etain){
	      err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAU[i], EPS[i], &P[i], &C[i], &T[i], pVAR, matFM[i], matFM[i], 
                     Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, P[i], TAUn[i], C[i],
	                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
			               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
			  //if(err){printf(" cycle=%d \n",n); return err;}
			  //if(err){ return 1;}
			  
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
      RC2[i]=C[i]*C[i]/(TAU[i]*TAU[i]);

      //printf("i=%d | TAU=%g, U=%g, E=%g, EPS=%g, C=%g, P=%g \n",i,TAU[i],U[i],E[i],EPS[i],C[i],P[i] );
    }
    t_end = omp_get_wtime();
    t_sum += t_end-t_start;
    
    // Calcul des valeurs ponctuelles à partir des valeurs moyennes (pour la prochaine itération)
    sum1=0;  sum3=0.; sum4=0.; sum5=0.;
    #pragma omp parallel for 
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RHO0[i]*E[i];
      sum3+=dx*RHO0[i]*U[i];
      sum4+=(X[i+1]-X[i])/TAU[i];
      sum5+=dx*RHO0[i]*TAU[i];
    }
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3; 
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    
	  t=t+dt;
	  n=n+1;
	  if(aff==1){  printf("n= %d, t= %g, dt= %g\n",n,t,dt); }

    // Calcul du pas de temps
    if(sch==11||sch==12){  // Gallice
      dt_new=(X[ideb+1]-X[ideb])/(TAU[ideb]*Zmax[ideb]);
      for(int j=ideb+1; j<=ifin; j++){
        dt_new = fmin(dt_new, (X[j+1]-X[j])/(TAU[j]*Zmax[j]) );
      }
    }
    else{  // Autres schemas
      dt_new=(X[ideb+1]-X[ideb])/C[ideb];
      for(int j=ideb+1; j<=ifin; j++){
        dt_new = fmin(dt_new, (X[j+1]-X[j])/C[j] );
      }
    }

    dt_new = CFL*dt_new; // pas de temps
    dt=pastemps(Tf,t,dt_new,dt);

  }
  //==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);
  if(sch==4){ printf("Nombre de mailles ou le système n'a pas été pas résoluble : %d\n",sum_sys_deg ); }
  printf("Temps EOS : %g\n",t_sum);

  // Valeurs moyennes pour tracé
  for(int i=ideb; i<=ifin; i++){
    Xc[i]=0.5*(X[i+1]+X[i]);
    RHO[i]=1./TAU[i];
  }


  // printf_sol
  err=print_sol(1, ideb, ifin, X, Xc, RHO, TAU, U, E, EPS, P, T, C, G, S, matFM);
  if(err){return err;}

  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}

  /*
  FILE* Res2;
  if((Res2 = fopen("nb_iter.txt", "w+")) != NULL){
    for(int i=ideb; i<=ifin; i++){  fprintf(Res2,"%d ",tabNiter[i]);  }
    fprintf(Res2, "0\n");
    for(int i=ideb; i<=ifin; i++){  fprintf(Res2,"%g ",tabEPSILON[i]);  }
    fprintf(Res2, "-1e-10\n");
    fclose(Res2);
  } 
  else{printf("erreur impression\n"); return 1;}
  */
  
  FILE* Res2;
  if((Res2 = fopen("nb_iter.txt", "w+")) != NULL){
    for(int i=ideb; i<=ifin; i++){  fprintf(Res2,"%d %g \n",tabNiter[i],tabEPSILON[i]);  }
    fclose(Res2);
  } 
  else{printf("erreur impression\n"); return 1;}

  free(U); free(TAU); free(E); 
  free(U_pr); free(TAU_pr); free(E_pr);   free(lambda_unuse);
  free(X); free(Xc); free(RHO0);
  free(Ustar); free(Pstar);  free(PUstar);  
  free(Ustar_01);   free(Pstar_01);
  free(P);  free(C); free(G); free(T); free(S); 
  free(EPS); free(RHO); 

  free(RT); free(RU); free(RR); free(RE);
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  free(RC2); free(PU);
  
  //cintétique 
  free(TAUn);  free(Un);

  free(matFM);    free(tabNiter);  free(tabEPSILON);  free(indFLUIDE);

  return 0;
}




int sgn(double v) {
  if (v < 0) return -1;
  if (v > 0) return 1;
  return 0;
}



// Fonction pour le choix du schéma
/*  Permet le choix de 
      le solveur ordre 1
        Voir main
      la méthode d'ordre 2
        methd=0   méthode ordre 1 (uniquement le solveur)
        methd=1   GAD
      le limiteur
        limiteur=0  pas de limiteur
        limiteur=1  MinMid
        limiteur=2  SuperBee
*/
int choix_sch(int sch, int* psolveur, int* pmethd, int* plimiteur){
  int solveur=0;
  int methd=0;
  int limiteur=0;

  if(sch==0 || sch==1 || sch==2){ solveur=sch; }             //Després, Jaouen, Solveur acoustique ordre 1 à deux états
  else if(sch==3){  solveur=2;  }                  // RK1 Matsuno + Solveur acoustique ordre 1 à deux états
  else if(sch==4){  solveur=3;  }                            //Extension faiblement non-linéaire à l'ordre 2
  else if(sch==5){ solveur=4;  }                             //Extension faiblement non-linéaire à l'ordre 3  
  else if(sch==6 || sch==7){  solveur=5;  }                  // Solveur de Dukovicz
  else if(sch==8 || sch==9 || sch==10){ solveur=6; }         // Solveur de Riemann exacte
  else if(sch==11||sch==12||sch==13){  solveur=7; }                   // Solveur de Gallice 1 et 2

  // Méthode GAD + solveur acoustique
  else if(sch==30 || sch==31 || sch==32){
    solveur=2; // solveur acoustique  
    methd=1; // GAD
    if(sch==30){limiteur=0;}
    else if(sch==31){limiteur=1;}
    else if(sch==32){limiteur=2;}
  }
  // Méthode GAD + solveur de Riemann exacte
  else if(sch==33 || sch==34 || sch==35){
    solveur=6; // solveur de Riemann exacte
    methd=1; // GAD
    if(sch==33){limiteur=0;}
    else if(sch==34){limiteur=1;}
    else if(sch==35){limiteur=2;}
  }

  *psolveur=solveur;
  *pmethd=methd;
  *plimiteur=limiteur;
}
