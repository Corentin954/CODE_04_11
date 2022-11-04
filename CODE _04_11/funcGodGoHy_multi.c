#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_multi.h"
//#include "head.h"


/*
  Schéma lagrangien GoHy d'ordre 2 et 3 en temps et arbitrairement élevée en espace
*/


/////////////////////  MAIN  ////////////////////
int funcGodGoHy_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int so, int EOS,
	                 int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
	
  double Ctriple=1e-10; // Valeur de C sur le point triple (sinon C=0)


  int ind_air, ind_etain; 
  ind_air=0; ind_etain=1;
	
	// variables
	double t, dt, dt_new;
	int n=0;
	int err=0;
	double rhoc, rhoc1;
	double maxC;
	char* W0;
	double mu;
	double phiplusU, phimoinsU;
  double phiplusP, phimoinsP;
	double phi_lim;
	double rhoc2, rhoc21;
	double r_plus, r_moins;
	double ms, mi;
	double sum1, sum2, sum3, sum4, sum5;
	double dX;
  int zone;
  double dUdx;
  double* pVAR=malloc(2*sizeof(double));

	double dx=(b-a)/nx; // pas d'espace

	// Calcul de la taille des tableaux
	int nbghostsPW=so; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
	int nbghosts=2*so; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
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
  double* RHO0pw=malloc((N)*sizeof(double)); // masse donc centrées aux mailles
	double* RHO0=malloc((N)*sizeof(double)); // masse donc centrées aux mailles
	//double** U=alloctabd(nx+2,3);           // solution  (tau, u, e)
 
  // Cell Average	
	double* RTAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* RU=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* RE=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* TAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* U=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* E=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)
  
  double* EPS=malloc((N)*sizeof(double));
  double* RHO=malloc((N)*sizeof(double));
  double* Pav=malloc((N)*sizeof(double)); 

  // Point-wise
  double* RTAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* RUpw=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* REpw=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* TAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* Upw=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* Epw=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* RHOpw=malloc((N)*sizeof(double));
  double* Ppw=malloc((N)*sizeof(double));
  double* EPSpw=malloc((N)*sizeof(double));
  double* Cpw=malloc((N)*sizeof(double));
  double* RC2pw=malloc((N)*sizeof(double));
  double* Gpw=malloc((N)*sizeof(double));
  double* Spw=malloc((N)*sizeof(double));
  double* PUpw=malloc((N)*sizeof(double));

  // Cinétique 
  double* TAUpw_n=malloc((N)*sizeof(double)); 
  double* Upw_n=malloc((N)*sizeof(double));   

  // FLux
  double* Ustar=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar=malloc((N+1)*sizeof(double));
  double* PUstar=malloc((N+1)*sizeof(double));

  double* Ustar_01=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar_01=malloc((N+1)*sizeof(double));

  // Diagnostic
  double* RRpw=malloc((N)*sizeof(double));
  double* RTpw=malloc((N)*sizeof(double));
  
  double* Tpw=malloc((N)*sizeof(double));  // température 
  double** matFM=alloctabd(N,3);          // fraction massique

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
  // variables défini aux centres des mailles

  double* Wpw=malloc(10*sizeof(double));
	double* Wav=malloc(10*sizeof(double));
	for(int i=0; i<=N-1; i++){
		//printf("  --i= %d\n", i);
	  x=(X[i]+X[i+1])/2.0;

    err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, Wpw);
	  err=u0_bar_multi(tst, indFLUIDE[i], ind_air, ind_etain, tst_air, a, b, X[i],X[i+1], Wav);
	  if(err){return err;}
	  
	  RHO0[i]=Wav[0];
    
    RTAU[i]=1.0;
    RU[i]=Wav[5];
    RE[i]=Wav[6];
    
    RTAUpw[i]=1.0;
	  
    EPS[i]=Wav[4];
    E[i]=Wav[3];
    
    Ppw[i]=Wpw[2];
	  Tpw[i]=Wpw[3];  
    matFM[i][0]=Wpw[4];
    matFM[i][1]=Wpw[5];
    matFM[i][2]=Wpw[6];
	}
  free(Wpw);
  free(Wav);


  for(int i=0; i<=N-1; i++){
    RHO0pw[i]=phi(so,&RHO0[i],Ck);
    RUpw[i]=phi(so,&RU[i],Ck);
    Upw[i]=RUpw[i]/RHO0pw[i];
    TAUpw[i]=1./RHO0pw[i];
    EPSpw[i]=phi(so,&EPS[i],Ck);
    Epw[i]=phi(so,&E[i],Ck);

	  if(indFLUIDE[i]==ind_etain){
	    err= G_S_cin_3_phase(epsilon, Ppw[i], TAUpw[i], EPSpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]);  if(err){return err;}
		}
	  else if(indFLUIDE[i]==ind_air){
      EGS(tst_air, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i]);
	    G_EGS(tst_air, TAUpw[i], EPSpw[i], &Gpw[i]);
    }

    if(Cpw[i]==0){Cpw[i]=Ctriple;  }
    RC2pw[i]=Cpw[i]*Cpw[i]/(TAUpw[i]*TAUpw[i]);

    PUpw[i]=Ppw[i]*Upw[i]; 
    //printf("i=%d : ind_fluide= %d, TAU= %g, U= %g, E= %g, EPS=%g, P= %g, RHO= %g, RC2=%g, C=%g, G=%g\n",i,indFLUIDE[i],TAUpw[i],Upw[i],Epw[i],EPSpw[i],Ppw[i],RHO0pw[i],RC2pw[i],Cpw[i],Gpw[i] );
	}

  /*
  for(int i=0; i<=N-1; i++){
    printf("i=%d | RTAU=%e, RU=%e, RE=%e, P=%e \n",i,RTAU[i],RU[i],RE[i],Ppw[i] );
  }
  */


  // minimum de XsC  (calcul de dt)
	maxC=Cpw[ideb];
	for (int j=ideb+1; j<=ifin; j++){
	  if (maxC<Cpw[j]){  maxC=Cpw[j];	}
  }
	dt = dt_first* CFL * dx / maxC; // pas de temps
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
  for(int i=idebPW; i<=ifinPW; i++){
    RRpw[i]=(X[i+1]-X[i])/TAUpw[i];
    RTpw[i]=RHO0pw[i]*TAUpw[i];
  }
  for(int i=ideb; i<=ifin; i++){
    sum1+=dx*RE[i];
    sum3+=dx*RU[i];
    sum4+=phi(so, &RRpw[i], Ckbar);
    sum5+=dx*phi(so, &RTpw[i], Ckbar);
  }
  Etot[0]=sum1;
  IMPUL[0]=sum3; 
  Mtot[0]=sum4;
  VOL[0]=sum5;



// debut de la boucle =============================================================================  
  while((t<Tf)&&(n<Nmax)){
    // Calcul du flux 
    if ( (sch == 13) || (sch == 14) || (sch == 15) ){ // GoHy d'ordre 2
      for(int i=ideb; i<=ifin+1; i++){
        // 2nd ordre
        Ustar[i]= f_inter(so, &Upw[i], rk);
        //Ustar[i]-= (dt/dx)*(Ppw[i] - Ppw[i-1])/(RHO0pw[i]+RHO0pw[i-1]);
        Ustar[i]+= (dt/2)*fdtu(&RHO0pw[i], &Ppw[i], dx);;
        
        Pstar[i]= f_inter(so, &Ppw[i], rk);
        //Pstar[i]-= (dt/2*dx)*(RC2pw[i]+RC2pw[i-1])*(Upw[i]-Upw[i-1])/(RHO0pw[i]+RHO0pw[i-1]);
        Pstar[i]+= (dt/2)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);

        //PUstar[i]=Pstar[i]*Ustar[i];
        PUstar[i]=f_inter(so, &PUpw[i], rk);
        PUstar[i]+= (dt/2)*f_inter(so, &Ppw[i], rk)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        PUstar[i]+= (dt/2)*f_inter(so, &Upw[i], rk)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);
        
        if( sch == 14 || sch== 15){ // ordre 2 avec limiteur 
          rhoc=Cpw[i-1]/TAUpw[i-1];
          rhoc1=Cpw[i]/TAUpw[i];
          // 1er ordre
          Ustar_01[i]=(rhoc*Upw[i-1]+rhoc1*Upw[i])/(rhoc+rhoc1) + (Ppw[i-1]-Ppw[i])/(rhoc+rhoc1) ;
          Pstar_01[i]=(Ppw[i]*rhoc+Ppw[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(Upw[i-1]-Upw[i]);

          r_moins=(Ustar[i]-Ustar[i-1]+1e-10)/(Upw[i]-Upw[i-1]+1e-10);
          //r_moins*=(RHO0pw[i]+RHO0pw[i-1])/(RHO0pw[i-1]);
          
          r_plus=(Ustar[i+1]-Ustar[i]+1e-10)/(Upw[i]-Upw[i-1]+1e-10);
          //r_plus*=(RHO0pw[i]+RHO0pw[i-1])/(RHO0pw[i]);

          if(sch==14){phi_lim=MinMod(fmin(r_plus, r_moins));}
          else if(sch==15){phi_lim=vanLeer(r_plus, r_moins);}
          else {printf("Erreur choix sch=%d\n",sch ); return 1;}
      
          Ustar[i]=(1.-phi_lim)*Ustar_01[i] + phi_lim*Ustar[i];
          Pstar[i]=(1.-phi_lim)*Pstar_01[i] + phi_lim*Pstar[i];
          PUstar[i]=(1.-phi_lim)*Pstar_01[i]*Ustar_01[i] + phi_lim*PUstar[i];

          //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
        }
      }
    }
    else if ( (sch == 20) || (sch == 21) || (sch == 22) ){ // GoHy ordre 3
      for(int i=ideb; i<=ifin+1; i++){
        Ustar[i]= f_inter(so, &Upw[i], rk);
        Ustar[i]+= (dt/2)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        Ustar[i]+= (dt*dt/6)*fdttu(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);
        
        Pstar[i]= f_inter(so, &Ppw[i], rk);
        Pstar[i]+= (dt/2)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);
        Pstar[i]+= (dt*dt/6)*fdttp(&RHO0pw[i], &RC2pw[i], &Gpw[i], &Ppw[i], &Upw[i], dx);
        
        /*
        PUstar[i]= f_inter(so, &PUpw[i], rk);
        PUstar[i]+= (dt*dt/3)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        PUstar[i]+= (dt/2)*((Ppw[i-1]+Ppw[i])/2)*(fdtu(&RHO0pw[i], &Ppw[i], dx) + (dt/3)*fdttu(&RHO0pw[i], &RC2pw[i], &Upw[i], dx));
        PUstar[i]+= (dt/2)*((Upw[i-1]+Upw[i])/2)*(fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx) + (dt/3)*fdttp(&RHO0pw[i], &RC2pw[i], &Gpw[i], &Ppw[i], &Upw[i], dx));
        */

        PUstar[i]= f_inter(so, &PUpw[i], rk);
        PUstar[i]+= (dt*dt/3)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        PUstar[i]+= (dt/2)*f_inter(so, &Upw[i], rk)*(fdtu(&RHO0pw[i], &Ppw[i], dx) + (dt/3)*fdttu(&RHO0pw[i], &RC2pw[i], &Upw[i], dx));
        PUstar[i]+= (dt/2)*f_inter(so, &Ppw[i], rk)*(fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx) + (dt/3)*fdttp(&RHO0pw[i], &RC2pw[i], &Gpw[i], &Ppw[i], &Upw[i], dx));
      
        if( sch == 21 || sch== 22){ // ordre 3 avec limiteur 
          rhoc=Cpw[i-1]/TAUpw[i-1];
          rhoc1=Cpw[i]/TAUpw[i];
          // 1er ordre
          Ustar_01[i]=(rhoc*Upw[i-1]+rhoc1*Upw[i])/(rhoc+rhoc1) + (Ppw[i-1]-Ppw[i])/(rhoc+rhoc1) ;
          Pstar_01[i]=(Ppw[i]*rhoc+Ppw[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(Upw[i-1]-Upw[i]);
          
          r_moins=(Ustar[i]-Ustar[i-1]+1e-10)/(Upw[i]-Upw[i-1]+1e-10);
          r_moins*=(RHO0pw[i]+RHO0pw[i-1])/(RHO0pw[i-1]);
          
          r_plus=(Ustar[i+1]-Ustar[i]+1e-10)/(Upw[i]-Upw[i-1]+1e-10);
          r_plus*=(RHO0pw[i]+RHO0pw[i-1])/(RHO0pw[i]);
          
          if(sch==21){phi_lim=MinMod(fmin(r_plus, r_moins));}
          else if(sch==22){phi_lim=vanLeer(r_plus, r_moins);}
          else {printf("Erreur choix sch=%d\n",sch ); return 1;}
          
          Ustar[i]=(1.- phi_lim)*Ustar_01[i] + phi_lim*Ustar[i];
          Pstar[i]=(1.- phi_lim)*Pstar_01[i] + phi_lim*Pstar[i];
          PUstar[i]=(1.- phi_lim)*Pstar_01[i]*Ustar_01[i] + phi_lim*PUstar[i];
          
          //printf("i=%d : U*=%g, P*=%g, phi=%g, r_plus=%g, r_moins=%g \n",i,Ustar[i],Pstar[i],phi,r_plus,r_moins );
        }
      }
    }
    else{printf("Erreur dans le choix de Ustar et Pstar (Gtot) sch=%d\n",sch); return 1;}


    // Calcul du pas de temps
    maxC=Cpw[ideb];
	  for(int j=ideb+1; j<=ifin; j++){
	    if(maxC<Cpw[j]){  maxC=Cpw[j];	}
	  }
	  dt_new = CFL * dx / maxC; // pas de temps	  
	  dt=pastemps(Tf,t,dt_new,dt);
    //printf(" maxC=%g\n",maxC );

	  // Mise à jour de la sol U
	  for(int i=ideb; i<=ifin; i++){
	    RTAU[i]+=(dt/dx)*(Ustar[i+1]-Ustar[i]);
	    RU[i]-=(dt/dx)*(Pstar[i+1]-Pstar[i]);
	    RE[i]-=(dt/dx)*(PUstar[i+1]-PUstar[i]);
	    
      //printf("i= %d : TAU=%g, U=%g, E=%g \n",i,TAU[i],U[i],E[i] );
	  }
	  
    
	  // Mise à jour du maillage
	  for(int i=ideb; i<=ifin+1; i++){
	  	X[i]+=Ustar[i]*dt;
	  }

	  // cond aux limites sur les valeurs moyennes
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, RTAU, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, RE, 0);
    if(err){return err;}
    err=condlim_Gtot(ind_cond, ideb,ifin, nbghosts, RU, iu);
	  if(err){return err;}

    // Cinétique
    for(int i=idebPW; i<=ifinPW; i++){
      TAUpw_n[i]=TAUpw[i];
      Upw_n[i]=Upw[i];
    }
    
    // Calcul des valeurs ponctuelles à partir des valeurs moyennes (pour la prochaine itération)
    for(int i=idebPW; i<=ifinPW; i++){
      RTAUpw[i]=phi(so, &RTAU[i], Ck);
      RUpw[i]=phi(so, &RU[i], Ck);
      REpw[i]=phi(so, &RE[i], Ck);

      TAUpw[i]=RTAUpw[i]/RHO0pw[i];
      Upw[i]=RUpw[i]/RHO0pw[i];
      Epw[i]=REpw[i]/RHO0pw[i];
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
    }


    // Mise à jour des tableaux
    for(int i=idebPW; i<=ifinPW; i++){ 
      dUdx=(Upw_n[i+1]-Upw_n[i-1])/2*dx;
      EPSpw[i]= Epw[i]-Upw[i]*Upw[i]/2;  // energie interne 

	    if(indFLUIDE[i]==ind_etain){
			  err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, epsilon, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i], &Tpw[i], pVAR, matFM[i], matFM[i], 
                     Ntau, dx, dt, Ttau, &zone, &tabNiter[i], &tabEPSILON[i], &test_consis_thermo[i], dUdx, Ppw[i], TAUpw_n[i], Cpw[i],
	                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
			               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
			  if(err){return err;}
			  
        Gpw[i]=pVAR[0];
        Spw[i]=pVAR[1];
        //err= G_S_cin_3_phase(epsilon, Ppw[i], TAUpw[i], EPSpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]); // S n'est pas utilisée
        //if(err){return err;}

			  //printf(" xA= %g, xB= %g, xC= %g\n",matFM[i][0],matFM[i][1],matFM[i][2] );
			}
		  else if(indFLUIDE[i]==ind_air){
        //EGS(tst_air, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i]);
		    //G_EGS(tst_air,TAUpw[i], EPSpw[i], &Gpw[i]);

        EGS_cin_phase(tst_air, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i],
                      cin_phase, dx, dt, Ttau, 
                      &Upw_n[i], Ppw[i], TAUpw_n[i], Cpw[i]);

        G_BIZ_cin_phase(TAUpw[i], EPSpw[i], &Gpw[i],
                        cin_phase, dx, dt, Ttau, 
                        &Upw_n[i], Ppw[i], TAUpw_n[i], Cpw[i]);
	    }
      
      if(Cpw[i]==0){Cpw[i]=Ctriple; }
      RC2pw[i]=Cpw[i]*Cpw[i]/(TAUpw[i]*TAUpw[i]);

      PUpw[i]=Ppw[i]*Upw[i];
      
      //printf("i=%d | TAUpw=%g, Upw=%g, Epw=%g, EPSpw=%g, P=%g, C=%g, RC2=%g, G=%g \n",i,TAUpw[i],Upw[i],Epw[i],EPSpw[i],Ppw[i],Cpw[i],RC2pw[i],Gpw[i] );
    }

    // Calcul des valeurs ponctuelles à partir des valeurs moyennes (pour la prochaine itération)
    sum1=0;  sum3=0.; sum4=0.; sum5=0.;
    for(int i=idebPW; i<=ifinPW; i++){
      RRpw[i]=(X[i+1]-X[i])/TAUpw[i];
      RTpw[i]=RHO0pw[i]*TAUpw[i];
    }
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RE[i];
      sum3+=dx*RUpw[i];
      sum4+=phi(so, &RRpw[i], Ckbar);
      sum5+=dx*phi(so, &RTpw[i], Ckbar);
    }
    Etot[n+1]=sum1;
    IMPUL[n+1]=sum3;  //printf("n=%d, impul=%g\n",n+1,sum3 );
    Mtot[n+1]=sum4;
    VOL[n+1]=sum5;
    
      
	  t=t+dt;
	  n=n+1;
	  if(aff==1){
	    printf("n= %d, t= %g, dt= %g\n",n,t,dt);
	  }
	  /*
    for(int i=ideb; i<=ifin; i++){
	    printf("i=%d : ind_fluide= %d, TAU= %g, U= %g, E= %g, P= %g, RHO= %g, C= %g, Ustar= %g, Pstar= %g\n",i,indFLUIDE[i],TAU[i],U[i],EPS[i],P[i],RHO0c[i],C[i],Ustar[i],Pstar[i] );
    }
    */
	}
//==========================================================================

  printf("Nombre d'itérations : %d\n",n);
  printf("Temps final : %g\n",t);

  for(int i=0; i<=N-1; i++){ 
    RHOpw[i]=1./TAUpw[i];
  }

  // Valeurs moyennes pour tracé
  for(int i=ideb; i<=ifin; i++){ 
    Xc[i]=fxc(so, &X[i], rk);
    EPS[i]=phi(so, &EPSpw[i], Ckbar);
    RHO[i]=phi(so, &RHOpw[i], Ckbar);
    TAU[i]=phi(so, &TAUpw[i], Ckbar);
    U[i]=phi(so, &Upw[i], Ckbar);
    E[i]=phi(so, &Epw[i], Ckbar);
    //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
  }

  // printf_sol
  err=print_sol(1, ideb, ifin, X, Xc, RHO, TAU, U, E, EPS, Ppw, Tpw, Cpw, Gpw, Spw, matFM);
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

  
  free(RU); free(RTAU); free(RE); 
  free(RUpw); free(RTAUpw); free(REpw); 
  free(Upw); free(TAUpw); free(Epw);  free(RHOpw); 
  free(X); free(Xc);
  free(Ustar); free(Pstar);  free(PUstar); 
  free(Ppw);  free(Cpw); free(Gpw); free(Tpw); free(Spw); 
  free(EPSpw); free(PUpw);

  free(EPS); free(U); free(TAU); free(E); free(RHO); 

  free(RRpw); free(RTpw);
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
}
