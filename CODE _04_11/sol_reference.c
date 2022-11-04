#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "head_multi.h"
#include "EOS_etain.h"
#include "sol_head.h"


//  Cree le fichier texte contenant le solution de reference
int main(int argc, char const *argv[]){
  
  int cas_test=atoi(argv[1]);

  int n,res; double err;
  double *X,*TAU,*RHO,*U,*EPS,*E,*P,*C,*G,*S,*T;
  double *lambdaBETA, *lambdaGAMMA, *lambdaLIQ;
  double Gamma,x,a,b,Tf,dx,eps;
  double epsilon=1e-12;

 
  // Calcul de solution exacte sur le meme nombre de mailles
  if(cas_test==0||cas_test==1){ // cas test de Sod ou LeBlanc
  	double p,u;
  	double mu2, Gamma;
  	double rhoL, uL, pL, rhoR, uR, pR, xdis;
  	double u12,p12;
    
    int ndetente=1e1;
    n=(ndetente+6);
    X=malloc(n*sizeof(double));
    TAU=malloc(n*sizeof(double));
    RHO=malloc(n*sizeof(double));
    U=malloc(n*sizeof(double));
    EPS=malloc(n*sizeof(double));
    E=malloc(n*sizeof(double));
    P=malloc(n*sizeof(double));
    C=malloc(n*sizeof(double));
    G=malloc(n*sizeof(double));
    T=malloc(n*sizeof(double));
    S=malloc(n*sizeof(double));

		/*   
		% Solution analytique pour une 1-onde de detente puis une 3-onde de choc 
		% x : points de calculs
		% xdis : position de la discontinuite
		% T : temps de la sol finale
		% wL, wR : donnee initiale
		function [W]=sol_ana_PU(x,xdis,T,wL,wR,p12,u12,gamma)
		*/

    // Sod
	  if(cas_test==0){
	    Gamma=1.4;
	    rhoL=1.0; rhoR=0.125;
	    uL=0; uR=0;
	    pL=1.0; pR=0.1;
	    a=0; b=1;
	    xdis=0.5*(b+a);
	    Tf=0.2;
	  }  
	  // LeBlanc
	  else if(cas_test==1){
	    Gamma=5./3;
	    rhoL=1; rhoR=1e-3;
	    uL=0; uR=0;
	    pL=(Gamma-1)/10; pR=(Gamma-1)*1e-12;
	    a=0; b=9;
	    xdis=(b+a)/3;
	    Tf=6;
	  }
	  else{printf("Erreur cas_test inconu =%d\n",cas_test ); return 1;}

    //mu2 = (Gamma-1.0)/(Gamma+1.0);
	  double cL=sqrt(Gamma*pL/rhoL);
	  double cR=sqrt(Gamma*pR/rhoR);

	  err=cal_PU(cas_test, &p12, &u12);
	  //printf("u*=%g et p*=%g\n",u12,p12 );

		//   LIVRE
		// vitesse de l'onde de choc   
		double aR=sqrt( rhoR*((Gamma+1)*p12/2 + (Gamma-1)*pR/2) ); // LIVRE
		double sigma = uR + aR/rhoR;
		//printf("aR=%g  sigma=%g \n",aR,sigma );

		// Etat 1 et 2 provenant du livre p119
		double c1=cL+(Gamma-1)*(uL-u12)/2;
		double rho1=Gamma*p12/(c1*c1);
		double rho2=rhoR*aR/(aR+rhoR*(uR-u12));
    //printf("rho1=%g  rho2=%g\n",rho1,rho2 );

		// Calcul des frontieres pour epsilon
		double eps_det_a=uL-cL,   x_det_a=xdis + eps_det_a*Tf;
		double eps_det_b=u12-c1,  x_det_b=xdis + eps_det_b*Tf;
		double eps_disc_cont=u12, x_disc_cont=xdis + eps_disc_cont*Tf;
		double eps_choc=sigma,    x_choc=xdis + eps_choc*Tf;
		//printf("frontières = %g  %g  %g  %g\n",eps_det_a,eps_det_b,eps_disc_cont,eps_choc);

    // etat gauche
    X[0]=a;
  	U[0]=uL;
	  C[0]=cL;
	  G[0]=0.5*(Gamma+1);
	  P[0]=pL;
	  RHO[0]=rhoL;
	  TAU[0]=1./rhoL;
	  EPS[0]=(1./(Gamma-1))*pL/rhoL;
	  E[0]=EPS[0] + 0.5*uL*uL;
    
    // détente
    dx=(x_det_b-x_det_a)/(ndetente-1);
    for(int i=1; i<1+ndetente; i++){
    	x=x_det_a+dx*(i-1);
    	//printf("i=%d x=%g\n",i,x );
    	eps=(x-xdis)/Tf;
    	X[i]=x;
    	U[i]=(Gamma-1)/(Gamma+1)*uL + 2/(Gamma+1)*(cL+eps);
		  C[i]=(Gamma-1)/(Gamma+1)*(uL-eps) + 2/(Gamma+1)*cL;
		  P[i]=pow(pL*(C[i]/cL),(2*Gamma/(Gamma-1)));
		  RHO[i]=Gamma*P[i]/(C[i]*C[i]);
		  TAU[i]=1./RHO[i];
	    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
	    E[i]=EPS[i] + 0.5*U[i]*U[i];
	    G[i]=0.5*(Gamma+1);
    }

    // discontinuité de contact
    int i=ndetente+1;
    X[i]=x_disc_cont - epsilon;
    U[i]=u12;
	  C[i]=c1;
	  P[i]=p12;
	  RHO[i]=rho1;
	  TAU[i]=1./RHO[i];
    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
    E[i]=EPS[i] + 0.5*U[i]*U[i];
	  G[i]=0.5*(Gamma+1);

    // choc
    i=ndetente+2;
    X[i]=x_disc_cont + epsilon;
    U[i]=u12;
	  P[i]=p12;
	  RHO[i]=rho2;
	  TAU[i]=1./RHO[i];
    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
    E[i]=EPS[i] + 0.5*U[i]*U[i];
    C[i]=sqrt(Gamma*P[i]*TAU[i]);
	  G[i]=0.5*(Gamma+1);

    i=ndetente+3;
    X[i]=x_choc - epsilon;
    U[i]=u12;
	  P[i]=p12;
	  RHO[i]=rho2;
	  TAU[i]=1./RHO[i];
    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
    E[i]=EPS[i] + 0.5*U[i]*U[i];
    C[i]=sqrt(Gamma*P[i]*TAU[i]);
	  G[i]=0.5*(Gamma+1);

    // etat droit
    i=ndetente+4;
    X[i]=x_choc + epsilon;
    U[i]=uR;
	  P[i]=pR;
	  RHO[i]=rhoR;
	  TAU[i]=1./RHO[i];
    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
    E[i]=EPS[i] + 0.5*U[i]*U[i];
    C[i]=sqrt(Gamma*P[i]*TAU[i]);
	  G[i]=0.5*(Gamma+1);

    i=ndetente+5;
    X[i]=b;
    U[i]=uR;
	  P[i]=pR;
	  RHO[i]=rhoR;
	  TAU[i]=1./RHO[i];
    EPS[i]=(1./(Gamma-1))*P[i]*TAU[i];
    E[i]=EPS[i] + 0.5*U[i]*U[i];
    C[i]=sqrt(Gamma*P[i]*TAU[i]);
	  G[i]=0.5*(Gamma+1);
  }
  else if(cas_test==3){ //Onde acoustique
    n=1e3;
    X=malloc(n*sizeof(double));
    TAU=malloc(n*sizeof(double));
    RHO=malloc(n*sizeof(double));
    U=malloc(n*sizeof(double));
    EPS=malloc(n*sizeof(double));
    E=malloc(n*sizeof(double));
    P=malloc(n*sizeof(double));
    C=malloc(n*sizeof(double));
    G=malloc(n*sizeof(double));
    T=malloc(n*sizeof(double));
    S=malloc(n*sizeof(double));

    Gamma=1.4;
	  double fac_eps=1e-8;
	  double kappa=2*M_PI;
	  double rho0=1.0; 
	  double p0=5./7;
	  
	  a=0.0; b=1.0; Tf=1.0;
	  dx=(b-a)/(n-1);
	  // Maillage
	  for(int i=0; i<n; i++){
	    x=a+dx*i;
	    X[i]=x;
	    RHO[i]= rho0 + fac_eps*sin(kappa*(x-Tf)); 
	    TAU[i]= 1./RHO[i];
	    U[i]= fac_eps*sin(kappa*(x-Tf));
	    P[i]= p0 + fac_eps*sin(kappa*(x-Tf));
	    EPS[i]= P[i]*TAU[i]/(Gamma-1);
	    E[i]= EPS[i] + 0.5*U[i]*U[i];
	    C[i]=sqrt(Gamma*P[i]*TAU[i]);
	    G[i]=0.5*(Gamma+1);
	  }
  }
  else if(cas_test>=101 && cas_test<=112){ // Choc étain
  	//printf(" Etain cas_test=%d (sol_reference)\n",cas_test );
    a=-1e-3;   b=1e-3; 
    //printf("Tf=%g\n",Tf );
    Tf=0.2e-6;
    //printf("Tf=%g\n", Tf);
    // IDEE : lancer le calcul de la vrai solution grace au Newton de L'hugoniot
    int nb_choc, nhalf;
    if(cas_test==101 || cas_test==106 || cas_test==107 || cas_test==108 || cas_test==109 || cas_test==110 || cas_test==111 || cas_test==112){
    	nb_choc=1;
    }
    else{ nb_choc=2; }
    
    n=2*(3*nb_choc+1);
    X=malloc(n*sizeof(double));
    TAU=malloc(n*sizeof(double));
    RHO=malloc(n*sizeof(double));
    U=malloc(n*sizeof(double));
    EPS=malloc(n*sizeof(double));
    E=malloc(n*sizeof(double));
    P=malloc(n*sizeof(double));
    C=malloc(n*sizeof(double));
    G=malloc(n*sizeof(double));
    T=malloc(n*sizeof(double));
    S=malloc(n*sizeof(double));
    lambdaBETA=malloc(n*sizeof(double));
    lambdaGAMMA=malloc(n*sizeof(double));
    lambdaLIQ=malloc(n*sizeof(double));

    double* lambda0=malloc(3*sizeof(double));
    double* lambda1=malloc(3*sizeof(double));
    double* lambdaf=malloc(3*sizeof(double));
    
    double* Etat0=malloc(13*sizeof(double));
    double* Etat1=malloc(13*sizeof(double));
    double* EtatFinal=malloc(13*sizeof(double));

		double Pf,uf,rhof,Vf,Ef,Tempf,Sf,Gibbsf,Gf,Df,Cf;

		double rho0,V0,E0,P0,T0,S0,Gibbs0,G0,C0,u0;
		double rho1,V1,E1,P1,T1,D1,u1,S1,Gibbs1,G1,C1;
		double choc1,choc2;
    
    res=etat(cas_test, Etat0, Etat1, EtatFinal);

	  // Remplisage des etats 0 et 1
	  /*
	  rho0=Etat0[0]; V0=1./rho0;
	  E0=Etat0[1];
	  P0=Etat0[2];
	  T0=Etat0[3];
	  S0=Etat0[4];     
	  Gibbs0=Etat0[5];
	  G0=Etat0[6];
	  C0=Etat0[7];
	  lambda0[0]=Etat0[9];
	  lambda0[1]=Etat0[10];
	  lambda0[2]=Etat0[11];
	  */

	  rho0=Etat0[0];
    V0=Etat0[1];
    //Etat0[2]=E0 + 0.5*u0*u0;
    E0=Etat0[3];
    P0=Etat0[4];
    T0=Etat0[5];
    C0=Etat0[6];
    G0=Etat0[7];
    S0=Etat0[8];
  	lambda0[0]=Etat0[9];
  	lambda0[0]=Etat0[10];
  	lambda0[0]=Etat0[11];
  	u0=Etat0[11];

    /*
	  rho1=Etat1[0]; V1=1./rho1;
	  E1=Etat1[1];
	  P1=Etat1[2];
	  T1=Etat1[3];
	  S1=Etat1[4];
	  Gibbs1=Etat1[5];
	  G1=Etat1[6];
	  C1=Etat1[7];
	  u1=Etat1[8];
	  lambda1[0]=Etat1[9];
	  lambda1[1]=Etat1[10];
	  lambda1[2]=Etat1[11];
	  D1=Etat1[12];
	  */

	  rho1=Etat1[0];
	  V1=Etat1[1];
	  //Etat1[2]=E1 + 0.5*u1*u1;
	  E1=Etat1[3];
	  P1=Etat1[4];
	  T1=Etat1[5];
	  C1=Etat1[6];
	  G1=Etat1[7];
	  S1=Etat1[8];
	  lambda1[0]=Etat1[9];
	  lambda1[1]=Etat1[10];
	  lambda1[2]=Etat1[11];
	  u1=Etat1[12];
	  D1=Etat1[13];

    /*
	  rhof=EtatFinal[0]; Vf=1./rhof;
	  Ef=EtatFinal[1];
	  Pf=EtatFinal[2];
	  Tempf=EtatFinal[3];
	  Sf=EtatFinal[4];
	  Gibbsf=EtatFinal[5];
	  Gf=EtatFinal[6];
	  Cf=EtatFinal[7];
	  uf=EtatFinal[8];
	  lambdaf[0]=EtatFinal[9];
	  lambdaf[1]=EtatFinal[10];
	  lambdaf[2]=EtatFinal[11];
	  Df=EtatFinal[12];
	  */

	  rhof=EtatFinal[0];
	  Vf=EtatFinal[1];
	  //Etat1[2]=E1 + 0.5*u1*u1;
	  Ef=EtatFinal[3];
	  Pf=EtatFinal[4];
	  Tempf=EtatFinal[5];
	  Cf=EtatFinal[6];
	  Gf=EtatFinal[7];
	  Sf=EtatFinal[8];
	  lambdaf[0]=EtatFinal[9];
	  lambdaf[1]=EtatFinal[10];
	  lambdaf[2]=EtatFinal[11];
	  uf=EtatFinal[12];
	  Df=EtatFinal[13];


 		
		// simple choc
		if(cas_test==101 || cas_test==106 || cas_test==107 || cas_test==108 || cas_test==109 || cas_test==110 || cas_test==111 || cas_test==112){
      
		  choc1=V0*sqrt((Pf-P0)/(V0-Vf))*Tf - uf*Tf;

		  //x=[a,-choc1-dx,-choc1+dx,choc1-dx,choc1+dx,b];
		  X[0]=a;  X[1]=-choc1-epsilon; X[2]=-choc1; X[3]=-choc1+epsilon; 
		  X[4]=choc1-epsilon; X[5]=choc1; X[6]=choc1+epsilon; X[7]=b;

	    TAU[0]=V0;   TAU[1]=V0;   TAU[2]=0.5*(Vf+V0);     TAU[3]=Vf;   
	    RHO[0]=rho0; RHO[1]=rho0; RHO[2]=0.5*(rhof+rho0); RHO[3]=rhof; 
	    U[0]=uf;     U[1]=uf;     U[2]=0.5*(uf);          U[3]=0;     
	    EPS[0]=E0;   EPS[1]=E0;   EPS[2]=0.5*(Ef+E0);     EPS[3]=Ef;   
	    P[0]=P0;     P[1]=P0;     P[2]=0.5*(Pf+P0);       P[3]=Pf;     
	    C[0]=C0; 		 C[1]=C0;     C[2]=0.5*(Cf+C0);			  C[3]=Cf;     
	    G[0]=G0; 		 G[1]=G0; 		G[2]=0.5*(Gf+G0); 			G[3]=Gf;     
	    S[0]=S0; 		 S[1]=S0; 		S[2]=0.5*(Sf+S0); 			S[3]=Sf;     
	    T[0]=T0; 		 T[1]=T0; 		T[2]=0.5*(Tempf-T0); 		T[3]=Tempf;     
	    lambdaBETA[0]=lambda0[0];  lambdaBETA[1]=lambda0[0];	lambdaBETA[2]=0.5*(lambdaf[0]+lambda0[0]);  lambdaBETA[3]=lambdaf[0];	
	    lambdaGAMMA[0]=lambda0[1]; lambdaGAMMA[1]=lambda0[1];	lambdaGAMMA[2]=0.5*(lambdaf[1]+lambda0[1]); lambdaGAMMA[3]=lambdaf[1];
	    lambdaLIQ[0]=lambda0[1];   lambdaLIQ[1]=lambda0[2];  	lambdaLIQ[2]=0.5*(lambdaf[2]+lambda0[2]);   lambdaLIQ[2]=lambdaf[2];
		}
		else {
		  choc1=V0*sqrt((P1-P0)/(V0-V1))*Tf - uf*Tf;
		  choc2=V1*sqrt((Pf-P1)/(V1-Vf))*Tf  + u1*Tf - uf*Tf;
		  //printf("choc1=%g choc2=%g\n",choc1,choc2 );
		  //printf("Tf=%g\n",Tf );
		  //printf(" V0=%g, P0=%g\n",V0,P0 );
		  //printf(" V1=%g, P1=%g, u1=%g\n",V1,P1,u1 );
		  //printf(" Vf=%g, Pf=%g, uf=%g\n",Vf,Pf,uf );
		  
		  //x=[a,-choc1-dx,-choc1,-choc1+dx,-choc2-dx,-choc2+dx,choc2-dx,choc2+dx,choc1-dx,choc1,choc1+dx,b];
		  X[0]=a; X[1]=-choc1-epsilon; X[2]=-choc1;  X[3]=-choc1+epsilon;  X[4]=-choc2-epsilon; X[5]=-choc2;  X[6]=-choc2+epsilon;
		  X[7]=choc2-epsilon; X[8]=choc2;  X[9]=choc2+epsilon; X[10]=choc1-epsilon; X[11]=choc1; X[12]=choc1+epsilon; X[13]=b;
      
		  TAU[0]=V0;   TAU[1]=V0;   TAU[2]=0.5*(V1+V0);     TAU[3]=V1;     TAU[4]=V1;     TAU[5]=0.5*(Vf+V1);     TAU[6]=Vf;   
	    RHO[0]=rho0; RHO[1]=rho0; RHO[2]=0.5*(rho1+rho0); RHO[3]=rho1;   RHO[4]=rho1;   RHO[5]=0.5*(rhof+rho1); RHO[6]=rhof;  
	    U[0]=uf;     U[1]=uf;     U[2]=0.5*(uf+(uf-u1));  U[3]=(uf-u1);  U[4]=(uf-u1);  U[5]=0.5*(uf-u1);       U[6]=0; 
	    EPS[0]=E0;   EPS[1]=E0;   EPS[2]=0.5*(E1+E0);     EPS[3]=E1;     EPS[4]=E1;     EPS[5]=0.5*(Ef+E1);     EPS[6]=Ef;   
	    P[0]=P0;     P[1]=P0;     P[2]=0.5*(P1+P0);       P[3]=P1;       P[4]=P1;       P[5]=0.5*(Pf+P1);       P[6]=Pf;       
	    C[0]=C0; 		 C[1]=C0;     C[2]=0.5*(C1+C0);			  C[3]=C1; 	     C[4]=C1;       C[5]=0.5*(Cf+C1);			  C[6]=Cf;     
	    G[0]=G0; 		 G[1]=G0; 		G[2]=0.5*(G1+G0); 		  G[3]=G1; 	     G[4]=G1;   	  G[5]=0.5*(Gf+G1); 	  	G[6]=Gf;     
	    S[0]=S0; 		 S[1]=S0; 		S[2]=0.5*(S1+S0); 		  S[3]=S1;   	   S[4]=S1;   		S[5]=0.5*(Sf+S1); 	  	S[6]=Sf;     
	    T[0]=T0; 		 T[1]=T0; 		T[2]=0.5*(T1+T0); 		  T[3]=T1;   	   T[4]=T1;   		T[5]=0.5*(Tempf+T1); 	  T[6]=Tempf;
	    lambdaBETA[0]=lambda0[0];  lambdaBETA[1]=lambda0[0];	lambdaBETA[2]=(1-lambda0[0]);  lambdaBETA[3]=lambda1[0];  lambdaBETA[4]=lambda1[0];	lambdaBETA[5]=0.5*(lambdaf[0]+lambda1[0]);  lambdaBETA[6]=lambdaf[0]; 	
	    lambdaGAMMA[0]=lambda0[1]; lambdaGAMMA[1]=lambda0[1];	lambdaGAMMA[2]=(1-lambda0[1]);  lambdaGAMMA[3]=lambda1[1]; lambdaGAMMA[4]=lambda1[1];	lambdaGAMMA[5]=0.5*(lambdaf[1]+lambda1[1]); lambdaGAMMA[6]=lambdaf[1]; 	
	    lambdaLIQ[0]=lambda0[2];   lambdaLIQ[1]=lambda0[2];	  lambdaLIQ[2]=(1-lambda0[2]);    lambdaLIQ[3]=lambda1[2];   lambdaLIQ[4]=lambda1[2];  	lambdaLIQ[5]=0.5*(lambdaf[2]+lambda1[2]);   lambdaLIQ[6]=lambdaf[2]; 	
		}
    
		nhalf=n/2;
    for(int i=0; i<nhalf; i++){
      TAU[n-1-i]=TAU[i];
      RHO[n-1-i]=RHO[i];
      U[n-1-i]=-U[i];
      EPS[n-1-i]=EPS[i];
      P[n-1-i]=P[i];
      C[n-1-i]=C[i];
      G[n-1-i]=G[i];
      S[n-1-i]=S[i];
      T[n-1-i]=T[i];
      lambdaBETA[n-1-i]=lambdaBETA[i];
      lambdaGAMMA[n-1-i]=lambdaGAMMA[i];
      lambdaLIQ[n-1-i]=lambdaLIQ[i];
    }
    for(int i=0; i<n; i++){
      E[i]=EPS[i]+0.5*U[i]*U[i];
    }
  }
  else{
  	//printf("Erreur cas test=%d (sol_reference.c)\n",cas_test );
  	return 1;
  }
  
  // Ecriture de la solution de reference
  /*  On ecrit toujours la solution de reference dans le meme fichier
        et lance le programme a chaque simulation
        (s'adapte mieux au chgmt de solution pour l'étain !)
  */

  // Ecriture dans un fichier des points (x,u(x)) de la solution sur grille centrée
  // si le schema est centrée i.e u est calculé au centres des mailles, on écrit u dans le fichier,
  // sinon on le remplace par des 0 (tout les autres types de schémas)
  /*
     sol_reference.txt                       sol_reference_u.txt                                                sol_reference_lambda.txt
  L1  : X  maillage                      L1  : Xc  maillage centrée                                         L1  : Xc  maillage centrée
  L2  : rho  densité                     L2  : u   vitesse (Godunov)  rien (autres types de schémas)        L2  : fraction massique de la phase BETA
  L3  : tau  volume spécifique                                                                              L3  : fraction massique de la phase GAMMA
  L4  : e   energie totale                                                                                  L4  : fraction massique de la phase LIQUIDE
  L5  : eps   energie interne
  L6  : P   pression
  L6  : T   temperature
  L7  : C   celerite du son
  L8  : G   derivee fondamentale
  L8  : S   entropie
  */
  FILE *Res1, *Res2, *Res3;
  if((Res1 = fopen("fichiers/GNUPLOT/sol_reference.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  if((Res2 = fopen("fichiers/GNUPLOT/sol_reference_u.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference_u.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  if((Res3 = fopen("fichiers/GNUPLOT/sol_reference_lambda.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference_lambda.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  for(int i=0; i<n; i++){ 
    fprintf(Res1,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",X[i],RHO[i],TAU[i],E[i],EPS[i],P[i],T[i],C[i],G[i],S[i]);
    fprintf(Res2,"%.15e %.15e\n",X[i],U[i]);
  }

  if(cas_test>=101 && cas_test<=112){
  	for(int i=0; i<n; i++){ 
	    fprintf(Res3,"%.15e %.15e %.15e %.15e\n",X[i],lambdaBETA[i],lambdaGAMMA[i],lambdaLIQ[i]);
	   }
  }
  
  return 0;
}

