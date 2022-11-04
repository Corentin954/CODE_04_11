#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "head_multi.h"
#include "EOS_etain.h"
#include "sol_head.h"


/*
Fonction qui lie P et U
l'indice k corespond aux valeurs a gauche (k=1) ou à droite (k=2)
*/
double f(double P, double pk, double rhok, double Gamma){
  double res, ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>= pk){
    res = (P-pk)/( rhok*ck*sqrt( (Gamma+1)*pik/(2*Gamma) + (Gamma-1)/(2*Gamma) ) );
  }
  else{
  	res = 2*ck*(pow(pik,(Gamma-1)/(2*Gamma)) -1 )/(Gamma-1);
  }
  return res;
}

/* Dérivee de f */
double fd(double P, double pk, double rhok, double Gamma){
  double num, den;
  double ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>= pk){
    num =  (Gamma+1)*pik + 3*Gamma-1;
    den =  4*Gamma*rhok*ck*sqrt( pow( (Gamma+1)*pik/(2*Gamma) + (Gamma-1)/(2*Gamma), 3) ) ;
  }
  else{
  	num = ck*pow(pik,(Gamma-1)/(2*Gamma));
  	den = Gamma*P;
  }
  return num/den;
}

// Fonction flux massique a
double a(double P, double pk, double rhok, double Gamma){
  double num, den;
  double ck, pik;
  pik=P/pk;
  ck=sqrt(Gamma*pk/rhok);
  if (P>=pk){
  	num=sqrt( rhok*( (Gamma+1)*P/2 + (Gamma-1)*pk/2 ) );  // formule 13.4 / 13.6
  	den=1.0;
  }
  else{
  	num = ((Gamma-1)/(2*Gamma))*rhok*ck*(1 - pik);
  	den = 1 - pow( P/pk,(Gamma-1)/(2*Gamma) ) ;            // formule 13.11 / 13.12
  }
  return num/den;
}

int cal_PU(int cas_test, double* pp, double* pu){
  
  double epsilon =1e-12;
  int n=0, Nmax=1e4;
  
  double Gamma;

  // donnes initiale
  double rhoG, rhoD;
  double uG, uD;
  double pG, pD;
  double cG, cD;

  double P0, delta;
  double num, den;
  double aG, aD; // equivalent du flux massique 

  // valeur à trouver 
  double P, U;
  

  // Sod
  if(cas_test==0){
    Gamma=1.4;
    rhoG=1.0; rhoD=0.125;
    uG=0; uD=0;
    pG=1.0; pD=0.1;
  }  
  // LeBlanc
  else if(cas_test==1){
    Gamma=5./3;
    rhoG=1; rhoD=1e-3;
    uG=0; uD=0;
    pG=(Gamma-1)/10; pD=(Gamma-1)*1e-12;
  }
  else{printf("Erreur cas_test inconu =%d (cal_PU/sol_err_exacte.c)\n",cas_test ); return 1;}
  

  cG=sqrt(Gamma*pG/rhoG); 
  cD=sqrt(Gamma*pD/rhoD);

  // Initialisation P0 
  P0 = (pG*rhoD*cD + pD*rhoG*cG + (uG-uD)*rhoG*cG*rhoD*cD )/( rhoG*cG + rhoD*cD );
  P=P0;
  delta=1;
  // Algo de Newton
  while(( fabs(delta)>=epsilon )&&(n<Nmax)){
    num = f(P,pG,rhoG, Gamma) + f(P,pD,rhoD, Gamma) - (uG - uD);
    den = fd(P,pG,rhoG, Gamma) + fd(P,pD,rhoD, Gamma);
    delta = num/den;
    P-=delta;
    n++;
  }
  /*
  printf("\n Newton :\n");
  printf("delta = %lf et n= %d\n",delta,n );
  */

  
  // Calcul de U
  U = (a(P,pG, rhoG, Gamma)*uG + a(P,pD, rhoD, Gamma)*uD + pG - pD)/(a(P,pG, rhoG, Gamma) + a(P,pD, rhoD, Gamma));
  
  *pp=P; 
  *pu=U;
}


// Valeur des états de choc | Solution exacte HUGONIOT etain
int etat(int cas_test, double* Etat0, double* Etat1, double* EtatFinal){

  //printf(" etat from func_sol_reference.c\n");

	double P,u,rho,V,E,T,S,Gibbs,G,D,C;
	double* lambda=malloc(3*sizeof(double));

	double rho0= 7287,             E0= 3.1279315986251e-10, P0=0,                  T0= 300,             S0=2.6926731550247e-05, Gibbs0=3.1279315986251e-10, G0=2.240097804912,  C0=2740.5538450273;
	double rho1= 8076.00681801076, E1= 51443.002832343,     P1= 7.6739859076668e9, T1= 404.14787213729, S1= 16.007216133979,    Gibbs1=995194.05593801,     G1=2.019367150304,  C1=3442.4648018356, D1=3283.1753189253, u1=320.75848494575;
	double rho2= 8288.3703764792,  E2= 77064.927809945,     P2= 8.4780616126608e9, T2= 354.88066418806, S2= 16.016695936291,    Gibbs2=1094267.3454656,     G2=1.7223706495414, C2=3430.2925202527, D2=3283.1753189253, u2=371.26594786928;
	double rho3= 8518.2842765723,  E3= 112608.3412405,      P3= 11.353849215428e9, T3= 388.4355077502,  S3= 25.068708561843,    Gibbs3=1435750.5849554,     G3=1.6743445292304, C3=3581.513727052,  D3=3283.1753189253, u3=474.56999745137;
	double rho4= 10247.863012022,  E4= 889466.82426922,     P4= 44.866636800633e9, T4= 1965.6747151802, S4= 305.21638123923,    Gibbs4=4667656.3444466,     G4=1.3493141278849, C4=4624.8467142991, D4=4616.3090031665, u4=1333.7667144364;
	double rho5= 10609.008312322,  E5= 1222604.7150191,     P5= 56.903490402702e9, T5= 2231.3107527235, S5= 383.61978969721,    Gibbs5=5730325.2776706,     G5=1.2127819931885, C5=4829.4181283246, D5=4993.8110485992, u5=1563.7165440188;
  
  double V0=1./rho0, V1=1./rho1, V2=1./rho2, V3=1./rho3, V4=1./rho4, V5=1./rho5;
  double u0=0;
  
  // Remplissage des etats 0 et 1
  Etat0[0]=rho0;
  Etat0[1]=1./rho0;
  Etat0[2]=E0 + 0.5*u0*u0;
  Etat0[3]=E0;
  Etat0[4]=P0;
  Etat0[5]=T0;
  Etat0[6]=C0;
  Etat0[7]=G0;
  Etat0[8]=S0;
  Etat0[9]=1;
  Etat0[10]=0;
  Etat0[11]=0;
  Etat0[11]=u0;

  Etat1[0]=rho1;
  Etat1[1]=1./rho1;
  Etat1[2]=E1 + 0.5*u1*u1;
  Etat1[3]=E1;
  Etat1[4]=P1;
  Etat1[5]=T1;
  Etat1[6]=C1;
  Etat1[7]=G1;
  Etat1[8]=S1;
  Etat1[9]=1;
  Etat1[10]=0;
  Etat1[11]=0;
  Etat1[12]=u1;
  Etat1[13]=D1;
  
  /*
	Etat1[0]=rho1;
  Etat1[1]=E1;
  Etat1[2]=P1;
  Etat1[3]=T1;
  Etat1[4]=S1;
  Etat1[5]=Gibbs1;
  Etat1[6]=G1;
  Etat1[7]=C1;
  Etat1[8]=u1;
  Etat1[9]=1;
  Etat1[10]=0;
  Etat1[11]=0;
  Etat1[12]=D1;
  */

	if(cas_test==101){
		lambda[0]=1; lambda[1]=0; lambda[2]=0;
		P=P1; u=u1; rho=rho1; V=V1; E=E1; T=T1; S=S1; Gibbs=Gibbs1; G=G1; D=D1; C=C1;
  }
	// ZONE MIXTE beta/gamma (1/2)
	// lambda=0.2
	if(cas_test==102){
		lambda[0]=0.8; lambda[1]=0.2; lambda[2]=0;
		P= 7837109295.2405; u= 331.046034536; rho= 8118.5453550428; V= 0.00012317477531598; E= 56474.765160832; T= 394.02227652445; S=16.007357604096; Gibbs=1021625.0939154; G=1.9570238747402; D=1963.3946949723; C=3442.6037758954;
  }
	// lambda=0.5 
	else if(cas_test==103){
		lambda[0]=0.5; lambda[1]=0.5; lambda[2]=0;
	  P= 8079944236.04; u= 346.31015149475; rho= 8182.2810751876; V= 0.00012221530778654; E= 64111.232454084; T= 379.07885359786; S=16.009038155771; Gibbs=1051485.1118034; G=1.8680657939052; D=1967.2771487405; C=3440.3507108391;
  }
	//lambda=0.8
  else if(cas_test==104){
		lambda[0]=0.2; lambda[1]=0.8; lambda[2]=0;
	  P= 8319940699.7071; u= 361.35917843614; rho= 8245.9493636798; V= 0.00012127166392806; E= 71850.5064501; T= 364.44976977374; S=16.013105491071; Gibbs=1080763.317928; G=1.7803046192314; D=1970.0261716925; C=3435.2332470752;
	}  

	// P2
	else if(cas_test==105){
		lambda[0]=0; lambda[1]=1; lambda[2]=0;
    rho= 8288.3703764792; E= 77064.927809945; P= 8.4780616126608e9; T= 354.88066418806; D= 3283.1753189253; u= 371.26594786928; S=16.016695936291; Gibbs=1094267.3454656; G=1.7223706495414; C=3430.2925202527;
  }
  // P3
  else if(cas_test==106){
		lambda[0]=0; lambda[1]=1; lambda[2]=0;
  	rho=8518.2842765723; E=112608.3412405; P=11.353849215428e9; T=388.4355077502; D=3283.1753189253; u=474.56999745137; S=25.068708561843; Gibbs=1435750.5849554; G=1.6743445292304; C=3581.513727052;
  }
  // P4
  else if(cas_test==107){
		lambda[0]=0; lambda[1]=1; lambda[2]=0;
	  rho=10247.863012022; E=889466.82426922; P=44.866636800633e9; T=1965.6747151802; D=4616.3090031665; u=1333.7667144364; S=305.21638123923; Gibbs=4667656.3444466; G=1.3493141278849; C=4624.8467142991;
  }

	// ZONE MIXTE gamma/liquide (2/3)
	// lambda=0.2
	else if(cas_test==108){
		lambda[0]=0; lambda[1]=0.8; lambda[2]=0.2;
  	P= 47365759158.054; u= 1383.1267402852; rho= 10326.095622402; V= 9.6842023991191e-05; E= 956519.78984598; T= 2023.7388172682; S=321.74209642994; Gibbs=4965690.1012973; G=1.3211817565948; D=4699.5227372271; C=4669.9996611123;
  }
	// lambda=0.5
	else if(cas_test==109){
		lambda[0]=0; lambda[1]=0.5; lambda[2]=0.5;
  	P= 51018537544.691; u= 1453.6692554144; rho= 10437.177697807; V= 9.5811341816103e-05; E= 1056577.1520685; T= 2105.7330981251; S=345.65235546068; Gibbs=5291859.4850195; G=1.2788079077289; D=4816.3011068108; C=4733.3837556223;
  }
	// lambda=0.8
	else if(cas_test==110){
		lambda[0]=0; lambda[1]=0.2; lambda[2]=0.8;
  	P= 54576171339.318; u= 1520.6957243034; rho= 10542.025726704; V= 9.4858429103141e-05; E= 1156257.7429573; T= 2182.5586666218; S=368.68221712734; Gibbs=5605149.5720574; G=1.2373877250962; D=4925.0650514303; C=4792.2995940202;
	}

	// P5 
	else if(cas_test==111){
		lambda[0]=0; lambda[1]=0; lambda[2]=1;
	  rho=10609.008312322; E=1222604.7150191; P=56.903490402702e9; T=2231.3107527235; D=4993.8110485992; u=1563.7165440188; S=383.61978969721; Gibbs=5730325.2776706; G=1.2127819931885; C=4829.4181283246;
  }
	else if(cas_test==112){
		lambda[0]=0; lambda[1]=1; lambda[2]=0;
  	u=900; rho=1./0.000106217; E=405000; P=2.61173e+10;
	  T=0; D=0; S=0; Gibbs=0; G=0; C=0;
  }

	EtatFinal[0]=rho;
  EtatFinal[1]=1./rho;
  EtatFinal[2]=E + 0.5*u*u;
  EtatFinal[3]=E;
  EtatFinal[4]=P;
  EtatFinal[5]=T;
  EtatFinal[6]=C;
  EtatFinal[7]=S;
  EtatFinal[8]=G;
  EtatFinal[9]=lambda[0];
  EtatFinal[10]=lambda[1];
  EtatFinal[11]=lambda[2];
  EtatFinal[12]=u;
  EtatFinal[13]=D;


  return 0;
}
