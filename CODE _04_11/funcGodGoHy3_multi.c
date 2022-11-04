#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_multi.h"


/*
  Schéma lagrangien de Godunov  GoHy d'ordre 3 
*/



// fonction pour passer d'une valeur point-wise à une valeur average sur la maille à l'ordre 3
double fpw_to_av(double* Ypw){
  return (11./12)*Ypw[0] + (1./24)*(Ypw[1] + Ypw[-1]);
}

// fonction pour passer d'une valeur average à une valeur point-wise sur la maille à l'ordre 3
double fav_to_pw(double* Yav){
  return (13./12)*Yav[0] - (1./24)*(Yav[1] + Yav[-1]);
}

// Reconstruction des positions des centres des mailles à l'ordre 3
double f_grid(double* Xd){
  return (9./16)*(Xd[1]+Xd[0]) + (-1./16)*(Xd[2]+Xd[-1]);
}


/////////////////////  MAIN  ////////////////////
int funcGodGoHy3_multi(int sch, int tst, double Tf, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS,
	                 int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){
	
  double Ctriple=1e-10; // Valeur de C sur le point triple (sinon C=0)

  int ordre_scheme=floor(sch/10);

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
	double phi;
	double rhoc2, rhoc21;
	double ms, mi;
	double sum1, sum2, sum3, sum4, sum5;
	double dX;
  int zone;

	double dx=(b-a)/nx; // pas d'espace

	// Calcul de la taille des tableaux
	int nbghostsPW=4; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
	int nbghosts=8; // nombre de mailles fantomes pour les variables du schéma (\rho, \tau, u, \epsilon, ekin)
  int N=nx+2*nbghosts;
  int ideb=nbghosts;
  int ifin=nbghosts+nx-1;
  int idebPW=ideb - nbghostsPW;
  int ifinPW=ifin+nbghostsPW;


  // Allocation des tableaux
  double* X=malloc((N+1)*sizeof(double));   // maillage
  double* Xc=malloc((N)*sizeof(double));    // maillage
	double* RHO0pw=malloc((N)*sizeof(double)); // masse donc centrées aux mailles
	//double** U=alloctabd(nx+2,3);           // solution  (tau, u, e)
	
  // Cell average
	double* RTAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
	double* RU=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* RE=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* TAU=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* U=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* E=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)


  // Point-wise
  double* RTAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* RUpw=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* REpw=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* TAUpw=malloc((N)*sizeof(double)); //  rho0.tau  (centrées aux mailles)
  double* Upw=malloc((N)*sizeof(double));   //  rho0.u (frontières)
  double* Epw=malloc((N)*sizeof(double));   //  rho0.epsilon (centrées aux mailles)

  double* Ppw=malloc((N)*sizeof(double));
  double* EPSpw=malloc((N)*sizeof(double));
  double* Cpw=malloc((N)*sizeof(double));
  double* RC2pw=malloc((N)*sizeof(double));
  double* Gpw=malloc((N)*sizeof(double));
  double* Spw=malloc((N)*sizeof(double));

  // Cell Average
  double* EPSav=malloc((N)*sizeof(double));

  double* Ustar=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar=malloc((N+1)*sizeof(double));
  double* PUstar=malloc((N+1)*sizeof(double));

  double* Ustar_01=malloc((N+1)*sizeof(double));  // aux frontières
  double* Pstar_01=malloc((N+1)*sizeof(double));
  
  double* Tpw=malloc((N)*sizeof(double));  // température 
  double** matFM=alloctabd(N,3);          // fraction massique

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

	double* W=malloc(10*sizeof(double));
	for(int i=0; i<=N-1; i++){
		//printf("  --i= %d\n", i);
	  x=(X[i]+X[i+1])/2.0;
    
	  err=u0_multi(tst, indFLUIDE[i], ind_air, ind_etain, a, b, x, W);
	  if(err){return err;}
	  
	  RHO0pw[i]=W[0];
    TAUpw[i]=1./W[0]; 
	  RTAUpw[i]=1.0; 
    Upw[i]=W[1];
	  RUpw[i]=W[0]*W[1];
	  Tpw[i]=W[3];  // if tst<100, T= ?
	  Ppw[i]=W[2];
    matFM[i][0]=W[4];
    matFM[i][1]=W[5];
    matFM[i][2]=W[6];
	  
	  if(indFLUIDE[i]==ind_etain){
	    err= epsilon_etain(1./W[0], W[2], Tpw[i], &EPSpw[i], matFM[i], nb_points, SA, SB, SC, SAB, SAC, SBC );
	    if(err){return err;}
	    //printf(" E= %g, V= %g, P= %g\n",EPS[i],1./W[0],W[2] );

	    err= G_S_cin_3_phase(epsilon, Ppw[i], TAUpw[i], EPSpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]); // S n'est pas utilisée
      if(err){return err;}

		}
	  else if(indFLUIDE[i]==ind_air){
      EPSpw[i]= epsilonEOS(tst_air, 1./W[0], W[2]);  // energie interne	  
      //printf("E= %g | V= %g, P= %g\n",EPS[i], TAU[i], P[i] );
      EGS(tst_air,TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i]);
	    G_EGS(tst_air,TAUpw[i], EPSpw[i], &Gpw[i]);

      //printf("P= %g | V= %g, E= %g\n",P[i], EPS[i], EPS[i] );
    }

    if(Cpw[i]==0){Cpw[i]=Ctriple;  }
    RC2pw[i]=Cpw[i]*Cpw[i]/(TAUpw[i]*TAUpw[i]);

    Epw[i]=EPSpw[i] + Upw[i]*Upw[i]/2;
    REpw[i]=RHO0pw[i]*Epw[i];
	  
    //printf("i=%d : ind_fluide= %d, TAU= %g, U= %g, E= %g, EPS=%g, P= %g, RHO= %g, RC2=%g, C=%g, G=%g\n",i,indFLUIDE[i],TAUpw[i],Upw[i],Epw[i],EPSpw[i],Ppw[i],RHO0pw[i],RC2pw[i],Cpw[i],Gpw[i] );
	  //printf("                       M= %g, C= %g, XsC= %g, Etot= %g\n",i,indFLUIDE[i],M[i],C[i],XsC[i], E[i] );
	}
	free(W);
   
  // Initialisation des valeurs moyennes sur la mailles
  if(ordre_scheme<=2){
    for(int i=ideb; i<=ifin; i++){
      RTAU[i]=RTAUpw[i];
      RU[i]=RUpw[i];
      RE[i]=REpw[i];
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
    }
  }
  else(ordre_scheme>=3){
    for(int i=ideb; i<=ifin; i++){
      RTAU[i]=fpw_to_av(&RTAUpw[i]);
      RU[i]=fpw_to_av(&RUpw[i]);
      RE[i]=fpw_to_av(&REpw[i]);
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
    }
  }

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
  for(int i=ideb; i<=ifin; i++){
  	//printf("i= %d\n",i );
    sum1+=dx*RE[i];
    sum3+=dx*RU[i];
    sum4+=(X[i+1]-X[i])*TAUpw[i];
    //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
    sum5+=dx*RHO0pw[i]*TAUpw[i];
  }
  Etot[0]=sum1;
  IMPUL[0]=sum3; 
  if (tst==0){
    IMPUL[0]-=0.9*t; 
  }
  Mtot[0]=sum4;  
  VOL[0]=sum5;


// debut de la boucle ============================================================================= 
  while((t<Tf)&&(n<Nmax)){
    // Calcul du flux 
    if(sch == 0){ // Despres
      for(int i=ideb; i<=ifin+1; i++){
        rhoc=(C[i-1]/TAU[i-1]+ C[i]/TAU[i])/2; 
        Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
        Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("  i=%d, Ustar=%g, Pstar=%g, rhoc=%g, TAU= %g , C= %g\n",i, Ustar[i],Pstar[i],rhoc,TAU[i],C[i]);
      }
    } 
    else if (sch == 1){ // Jaouen 2001
      for(int i=ideb; i<=ifin+1; i++){
        rhoc2=C[i-1]*C[i-1]/TAU[i-1]; 
        rhoc21=C[i]*C[i]/TAU[i];
        ms=fmax(rhoc21,rhoc2);
        mi=fmax(TAU[i],TAU[i-1]);
        rhoc=sqrt(ms/mi);
        Ustar[i]=(U[i-1]+U[i])/2 + (P[i-1]-P[i])/(2*rhoc) ;
        Pstar[i]=(P[i]+P[i-1])/2 + (rhoc/2)*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
      }
    }
    else if (sch == 2){ // Solveur Acoustique aux faces ordre 1
      for(int i=ideb; i<=ifin+1; i++){
        rhoc=C[i-1]/TAU[i-1]; 
        rhoc1=C[i]/TAU[i];
        Ustar[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
        Pstar[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]) ;
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
      }
    }
    else if ( (sch == 10) || (sch == 11) || (sch == 12) ){ // GAD d'ordre 2
      for(int i=ideb; i<=ifin+1; i++){
        rhoc=C[i-1]/TAU[i-1]; 
        rhoc1=C[i]/TAU[i];
        // 1er ordre
        Ustar[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
        Pstar[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]) ;
        // 2nd ordre
        mu= 0.5*(rhoc1+rhoc)*dt/(0.5*dx*(RHO0pw[i-1]+RHO0pw[i]));

        if (sch == 10 ){  // ordre 2 sans limiteur
          phimoinsU=1.0;  phimoinsP=1.0;
          phiplusU=1.0;   phiplusP=1.0;
        }
        else if( sch == 11 ){ // ordre 2 avec limiteur minmod
          phimoinsU=MinMod((Ustar[i+1]-U[i])/(Ustar[i]-U[i-1]));
          phimoinsP=MinMod((Pstar[i+1]-P[i])/(Pstar[i]-P[i-1]));

          phiplusU=MinMod((U[i-1]-Ustar[i-1])/(U[i]-Ustar[i]));
          phiplusP=MinMod((P[i-1]-Pstar[i-1])/(P[i]-Pstar[i]));
        }
        else if( sch == 12 ){ // ordre 2 avec limiteur superbee
          // printf("superbee\n");
          phimoinsU=Superbee((Ustar[i+1]-U[i])/(Ustar[i]-U[i-1]));
          phimoinsP=Superbee((Pstar[i+1]-P[i])/(Pstar[i]-P[i-1]));

          phiplusU=Superbee((U[i-1]-Ustar[i-1])/(U[i]-Ustar[i]));
          phiplusP=Superbee((P[i-1]-Pstar[i-1])/(P[i]-Pstar[i]));
        }
        Ustar[i]+=(0.5)*(1.0-mu)*( phiplusU*(U[i]-Ustar[i]) - phimoinsU*(Ustar[i]-U[i-1]) );
        Pstar[i]+=(0.5)*(1.0-mu)*( phiplusP*(P[i]-Pstar[i]) - phimoinsP*(Pstar[i]-P[i-1]) );
        PUstar[i]=Pstar[i]*Ustar[i];
        //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
      }
    }
    else if ( (sch == 13) || (sch == 14) || (sch == 15) ){ // GoHy d'ordre 2
      for(int i=ideb; i<=ifin+1; i++){
        // 2nd ordre
        Ustar[i]= (U[i-1] + U[i])/2;
        Ustar[i]-= (dt/dx)*(P[i] - P[i-1])/(RHO0[i]+RHO0[i-1]);
        
        Pstar[i]= (P[i-1] + P[i])/2;
        Pstar[i]-= (dt/2*dx)*(RC2[i]+RC2[i])*(U[i]-U[i-1])/(RHO0[i]+RHO0[i-1]);
      

        else if( sch == 13 || sch= 14){ // ordre 2 avec limiteur 
          rhoc=C[i-1]/TAU[i-1]; 
          rhoc1=C[i]/TAU[i];
          // 1er ordre
          Ustar_01[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
          Pstar_01[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]);
       
          r_moins=(Ustar[i]-Ustar[i-1])/(Upw[i]-Upw[i-1]);
          r_moins*=(RHO0pw[i]-RHO0pw[i-1])/(RHO0pw[i-1]);
          
          r_plus=(Ustar[i+1]-Ustar[i])/(Upw[i]-Upw[i-1]);
          r_plus*=(RHO0pw[i]-RHO0pw[i-1])/(RHO0pw[i]);

          if(sch==13){phi=MinMod(fmax(r_plus,r_moins));}
          else if(sch==14){phi=vanLeer(r_plus, r_moins);}
          else {printf("Erreur choix sch=%d\n",sch ); return 1;}
      
          Ustar[i]=(1-phi)*Ustar_01 + phi*Ustar[i];
          Pstar[i]=(1-phi)*Pstar_01 + phi*Pstar[i];
          PUstar[i]=Pstar[i]*Ustar[i];

          //printf("i=%d, Ustar=%f, Pstar=%f, rhoc=%f, rhoc1=%f P[i-1]=%f, P[i]=%f\n",i, Ustar[i],Pstar[i],rhoc,rhoc1,P[i-1], P[i] );
        }
      }
    }
    else if ( (sch == 20) || (sch == 21) || (sch == 22) ){ // GoHy ordre 3
      for(int i=ideb; i<=ifin+1; i++){
        Ustar[i]= (9./16)*(Upw[i-1] + Upw[i]) - (1./16)*(Upw[i-2] + Upw[i+1]);
        Ustar[i]+= (dt/2)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        Ustar[i]+= (dt*dt/6)*fdttu(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);
        
        Pstar[i]= (9./16)*(Ppw[i-1] + Ppw[i]) - (1./16)*(Ppw[i-2] + Ppw[i+1]);
        Pstar[i]+= (dt/2)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx);
        Pstar[i]+= (dt*dt/6)*fdttp(&RHO0pw[i], &RC2pw[i], &Gpw[i], &Ppw[i], &Upw[i], dx);
        
        PUstar[i]= (9./16)*(Ppw[i-1]*Upw[i-1] + Ppw[i]*Upw[i]) - (1./16)*(Ppw[i-2]*Upw[i-2] + Ppw[i+1]*Upw[i+1]);
        PUstar[i]+= (dt*dt/3)*fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx)*fdtu(&RHO0pw[i], &Ppw[i], dx);
        PUstar[i]+= (dt/2)*((Ppw[i-1]+Ppw[i])/2)*(fdtu(&RHO0pw[i], &Ppw[i], dx) + (dt/2)*fdttu(&RHO0pw[i], &RC2pw[i], &Upw[i], dx));
        PUstar[i]+= (dt/2)*((Upw[i-1]+Upw[i])/2)*(fdtp(&RHO0pw[i], &RC2pw[i], &Upw[i], dx) + (dt/2)*fdttp(&RHO0pw[i], &RC2pw[i], &Gpw[i], &Ppw[i], &Upw[i], dx));
      
        else if( sch == 21 || sch= 22){ // ordre 3 avec limiteur 
          rhoc=C[i-1]/TAU[i-1]; 
          rhoc1=C[i]/TAU[i];
          // 1er ordre
          Ustar_01[i]=(rhoc*U[i-1]+rhoc1*U[i])/(rhoc+rhoc1) + (P[i-1]-P[i])/(rhoc+rhoc1) ;
          Pstar_01[i]=(P[i]*rhoc+P[i-1]*rhoc1)/(rhoc+rhoc1) + (rhoc*rhoc1/(rhoc+rhoc1))*(U[i-1]-U[i]);
       
          r_moins=(Ustar[i]-Ustar[i-1])/(Upw[i]-Upw[i-1]);
          r_moins*=(RHO0pw[i]-RHO0pw[i-1])/(RHO0pw[i-1]);
          
          r_plus=(Ustar[i+1]-Ustar[i])/(Upw[i]-Upw[i-1]);
          r_plus*=(RHO0pw[i]-RHO0pw[i-1])/(RHO0pw[i]);

          if(sch==13){phi=MinMod(fmax(r_plus,r_moins));}
          else if(sch==14){phi=vanLeer(r_plus, r_moins);}
          else {printf("Erreur choix sch=%d\n",sch ); return 1;}
      
          Ustar[i]=(1-phi)*Ustar_01 + phi*Ustar[i];
          Pstar[i]=(1-phi)*Pstar_01 + phi*Pstar[i];
          PUstar[i]=(1-phi)*Pstar_01*Ustar_01[i] + phi*PUstar[i];

        }
      //printf("i=%d : U*=%g, P*=%g, PU*=%g \n",i,Ustar[i],Pstar[i],PUstar[i] );
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
	    TAU[i]+=(dt/dx*RHO0pw[i])*(Ustar[i+1]-Ustar[i]);
	    U[i]-=(dt/dx*RHO0pw[i])*(Pstar[i+1]-Pstar[i]);
	    E[i]-=(dt/dx*RHO0pw[i])*(PUstar[i+1]-PUstar[i]);
	    
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
	  }
	  
    
	  // Mise à jour du maillage
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
    
    // Calcul des valeurs ponctuelles à partir des valeurs moyennes (pour la prochaine itération)
    if(ordre_scheme<=2){
      for(int i=idebPW; i<=ifinPW; i++){
        TAUpw[i]=TAU[i];
        Upw[i]=U[i];
        Epw[i]=E[i];
      }
    else if(ordre_scheme=>3){
      for(int i=idebPW; i<=ifinPW; i++){
        TAUpw[i]=fav_to_pw(&TAU[i]);
        Upw[i]=fav_to_pw(&U[i]);
        Epw[i]=fav_to_pw(&E[i]);
      }
    }

    // Mise à jour des tableaux
    for(int i=idebPW; i<=ifinPW; i++){ 

      EPSpw[i]= Epw[i]-Upw[i]*Upw[i]/2;  // energie interne 

	    if(indFLUIDE[i]==ind_etain){
		    if(err){return err;}
			  err=fPTC_cin_phase(0, EOS, cin_phase, sch_cin_phase, TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i], &Tpw[i], matFM[i], matFM[i], Ntau, dx, dt, Ttau, &zone, &Upw[i], Ppw[i], C[i],
	                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
			               SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
			  if(err){return err;}
			  
        err= G_S_cin_3_phase(epsilon, Ppw[i], TAUpw[i], EPSpw[i], matFM[i], &Cpw[i], &Gpw[i], &Spw[i]); // S n'est pas utilisée
        if(err){return err;}

			  //printf(" xA= %g, xB= %g, xC= %g\n",matFM[i][0],matFM[i][1],matFM[i][2] );
			}
		  else if(indFLUIDE[i]==ind_air){
        EGS(tst_air,TAUpw[i], EPSpw[i], &Ppw[i], &Cpw[i]);
		    G_EGS(tst_air,TAUpw[i], EPSpw[i], &Gpw[i]);
	    }
      
      if(Cpw[i]==0){Cpw[i]=Ctriple; }
      RC2pw[i]=Cpw[i]*Cpw[i]/(TAUpw[i]*TAUpw[i]);
      
      //printf("i=%d | TAUpw=%g, Upw=%g, Epw=%g, EPSpw=%g, P=%g, C=%g, RC2=%g, G=%g \n",i,TAUpw[i],Upw[i],Epw[i],EPSpw[i],Ppw[i],Cpw[i],RC2pw[i],Gpw[i] );
    }


    sum1=0;  sum3=0.; sum4=0.; sum5=0.;
    for(int i=ideb; i<=ifin; i++){
      sum1+=dx*RHO0pw[i]*Epw[i];
      sum3+=dx*Upw[i];
	    sum4+=(X[i+1]-X[i])*TAUpw[i];
	    //sum4+=dx*phi(so,&RHO0c[i], Ckbar);
	    sum5+=dx*RHO0pw[i]*TAUpw[i];
	  }
	  Etot[n+1]=sum1;
	  IMPUL[n+1]=sum3; 
	  if (tst==0){
	    IMPUL[n+1]-=0.9*t; 
	  }
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

  // Valeurs moyennes pour tracé
  if(ordre_scheme<=2){
    for(int i=ideb; i<=ifin; i++){ 
      Xc[i]=(X[i+1]+X[i])/2;
      Pav[i]=Ppw[i];
      EPSav[i]=EPSpw[i];
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
    }
  else if(ordre_scheme=>3){
    for(int i=ideb; i<=ifin; i++){ 
      Xc[i]=f_grid(&X[i]);
      Pav[i]=fpw_to_av(&Ppw[i]);
      EPSav[i]=fpw_to_av(&EPSpw[i]);
      //printf("i= %d : RTAU=%g, RU=%g, RE=%g \n",i,RTAU[i],RU[i],RE[i] );
    }
  }

  // printf_sol
  err=print_sol(1, ideb, ifin, X, Xc, TAU, U, E, Pav, EPSav, Cpw, Gpw, Spw, matFM);
  if(err){return err;}


  err=print_nrj(n, Etot, IMPUL, Mtot, VOL);
  if(err){return err;}

  
  free(RU); free(RTAU); free(RE); 
  free(RUpw); free(RTAUpw); free(REpw); 
  free(Upw); free(TAUpw); free(Epw); 
  free(X); free(Xc);
  free(Ustar); free(Pstar);  free(PUstar); 
  free(Ppw);  free(Cpw); free(Gpw); free(Tpw); free(Spw); 
  free(EPSpw);

  free(Uav); free(TAUav); free(Eav); 
  free(Pav);  free(Cav); free(Gav); free(Sav); 
  free(EPSav);

  
  free(Etot); 
  free(IMPUL);
  free(Mtot);
  free(VOL);

  return 0;
}
