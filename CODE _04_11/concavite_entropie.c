#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EOS_etain.h"



double test(int phase, double V){ 
	double K0, N0, gamma0, Cvr, theta0, T0, P0,rho0, v0, E0, Sr;
	double ur=1.2;
	coeff(phase, &K0, &N0, &gamma0, &Cvr, &theta0, &T0, &P0, &rho0, &v0, &E0, &Sr);
  
  //return gamma0*rho0*ur*THETA_prime(V, phase) - (1./Cvr)*Ps_prime(V, phase);
  return  - (1./Cvr)*Es_prime(V, phase) + ur*THETA_prime(V, phase);
}

void test_cst(int phase){ 
	double K0, N0, gamma0, Cvr, theta0, T0, P0,rho0, v0, E0, Sr;
	double ur=1.2;
	coeff(phase, &K0, &N0, &gamma0, &Cvr, &theta0, &T0, &P0, &rho0, &v0, &E0, &Sr);
  
  double tau1=-log(Cvr*ur*theta0*rho0*gamma0*gamma0) + log(K0);
  //printf("	log(Cvr*ur*theta0*rho0*gamma0*gamma0)=%g  log(K0)=%g\n",log(Cvr*ur*theta0*rho0*gamma0*gamma0),log(K0));
  tau1*=1./(N0+1- gamma0);
  //printf("	1./(N0+1- gamma0)=%g\n",1./(N0+1- gamma0));
  tau1+=+1.;
  tau1*=1./rho0;

  printf("	-ur*theta0*(rho0*gamma0)*(rho0*gamma0) + rho0*K0/Cvr : %g\n", -ur*theta0*(rho0*gamma0)*(rho0*gamma0) + rho0*K0/Cvr);
  printf("	1 - gamma0/(N0+1) : %g\n", 1 - gamma0/(N0+1) );
  printf(" tau_1 = (1./rho0)((ln(Cvr.ur.theta0.rho0.Gamma0²)-ln(K0))/(N0+1-gamma0) -1): %g\n", tau1);

  //Majoration de \frac{-1}{C_{Vr}}E_s'(\tau) + u_r\Theta'(\tau)  :
  double TAUstar=log(ur*(rho0*gamma0)*(rho0*gamma0)*theta0*Cvr) - log(rho0*K0);
  TAUstar*=1./(N0+1- gamma0);
  TAUstar = 1 - TAUstar;
  TAUstar*=1./rho0;

  double Fstar=(-1./Cvr)*Es_prime(TAUstar,phase) + ur*THETA_prime(TAUstar,phase);

  printf(" tau*=%g  et f(tau*)=%g\n",TAUstar,Fstar );

  
  // Condition 1   -ur.Theta.rho0.K0 + E's(tau_1)f(tau*) <? 0
  printf("Condition 1   -ur.Theta.rho0.K0 + E's(tau_1)f(tau*) (<? 0) =%g\n",ur*theta0*rho0*K0 + Es_prime(tau1,phase)*Fstar );
}

double cond1(int phase, double V, double E){
  double T=fT(V,E,phase);
  double P=fP(V,E,phase);
  double dvP=dPdV(V,phase);
  double dvT=dTdV(V, phase);

  return T*dvP - P*dvT;
}

double cond2(int phase){
  return dTdE(phase);
}

double cond3(int phase, double V){
  double deT=dTdE(phase);
  double deP=dPdE(phase);
  double dvP=dPdV(V,phase);
  double dvT=dTdV(V,phase);

  return dvT*deP - dvP*deT;
}

double def_entropie(int phase, double V){
	double K0, N0, gamma0, Cvr, theta0, T0, P0,rho0, v0, E0, Sr;
	double ur=1.2;
	coeff(phase, &K0, &N0, &gamma0, &Cvr, &theta0, &T0, &P0, &rho0, &v0, &E0, &Sr);
  
  return fEs(V, phase) - ur*Cvr*THETA(V,phase);
}

double sgn(double x){
	return x/fabs(x);
}


int main(int argc, char const *argv[]){
	double Nv=2e2, Ne=3e2;
	double V,E;

  printf("phase beta :\n");
  test_cst(1);
  printf("phase gamma :\n");
  test_cst(2);
  printf("phase liquide :\n");
  test_cst(3);
  
  /*
  double Vmin=110e-6, Vmax=145e-6;
  double Emin=-5e3,   Emax=3e5;
  */

  //double Vmin=70e-6, Vmax=230e-6;
  double Vmin=70e-6, Vmax=7000e-6;
  double Emin=-10e4, Emax=3e7;

  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);

  FILE* fileCOND_V, *fileCOND_E;
  if((fileCOND_V = fopen("fichiers/cond_V.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  if((fileCOND_E = fopen("fichiers/cond_E.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  
	
  FILE *fileCOND3, *fileCOND3_test;
  if((fileCOND3 = fopen("fichiers/cond3_concavite.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  if((fileCOND3_test = fopen("fichiers/cond3_test_concavite.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  
  FILE *fileDEF_ENTROPIE;
  if((fileDEF_ENTROPIE = fopen("fichiers/def_entropie.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  

  
  FILE *fileCOND1_beta, *fileCOND1_gamma, *fileCOND1_liq;
  if((fileCOND1_beta = fopen("fichiers/cond1_beta_concavite.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  if((fileCOND1_gamma = fopen("fichiers/cond1_gamma_concavite.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  if((fileCOND1_liq = fopen("fichiers/cond1_liq_concavite.txt", "w+")) == NULL){printf("erreur ouverture fichier\n");return 1;}
  
  
  V=Vmin;
  for(int i=0; i<Nv+1; i++){
  	fprintf(fileCOND3, "%.15e %.15e %.15e\n",cond3(1,V),cond3(2,V),cond3(3,V) );
  	fprintf(fileCOND3_test, "%.15e %.15e %.15e\n",test(1,V),test(2,V),test(3,V) );
  	fprintf(fileDEF_ENTROPIE, "%.15e %.15e %.15e\n",def_entropie(1,V),def_entropie(2,V),def_entropie(3,V) );
  	V+=hv;
  }
  
  
  E=Emin;
  V=Vmin;
  for(int i=0; i<Nv+1; i++){
  	for(int j=0; j<Ne+1; j++){
  	  fprintf(fileCOND1_beta, "%.15e ",sgn(cond1(1,V,E)) );
  	  fprintf(fileCOND1_gamma, "%.15e ",sgn(cond1(2,V,E)) );
  	  fprintf(fileCOND1_liq, "%.15e ",sgn(cond1(3,V,E)) );
  	  E+=he;
  	}
  	fprintf(fileCOND1_beta, "\n" );
  	fprintf(fileCOND1_gamma, "\n");
  	fprintf(fileCOND1_liq, "\n");

  	V+=hv;
  }
  
  

	V=Vmin; 
  for(int i=0; i<Nv+1; i++){
  	fprintf(fileCOND_V, "%.15e\n",V );
  	V+=hv;
  }

  E=Emin;
  for(int i=0; i<Ne+1; i++){
  	fprintf(fileCOND_E, "%.15e\n",E );
  	E+=he;
  }


  fclose(fileCOND3); fclose(fileCOND3_test); 
  //fclose(fileCOND1_beta);  fclose(fileCOND1_gamma);  fclose(fileCOND1_liq);   
  fclose(fileCOND_V);  fclose(fileCOND_E);
  fclose(fileDEF_ENTROPIE); 
  

	return 0;
}