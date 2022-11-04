#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_multi.h"

/*  /////////////////////////////////////////////////////////////////////////////////

           Equation d'etat du GAZ PARFAIT

p = (gamma-1)*epsilon/tau
epsilon = tau/((gamma-1)*epsilon)
*/


//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonGP(double Gamma, double tau, double p){
  return tau*p/(Gamma-1) ;
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double egsGP(double Gamma, double tau, double epsilon, double* pP, double* pC){
  *pP=(Gamma-1)*epsilon/tau;
  *pC= sqrt(Gamma*(*pP)*tau);
}

// renvoie ¨Pres et Cres
double G_GP(double Gamma, double tau, double epsilon, double* pG){
  *pG= (Gamma+1)/2;
}

/* /////////////////////////////////////////////////////////////////////////////////
    
    Equation d'etat du BIZARIUM 

p = p_k0 + (sigma0/tau0).(epsilon - epsilon_k0)
epsilon = epsilon_k0 + (tau0/sigma0).(p - p_k0)  
*/

// données 
double rho0 = 1e4; double tau0=1./1e4;
double K0 = 1e11;
double Cv0 = 1e3;
double T0 = 300;
double eps0 = 0;
double sigma0 = 1.5;
double s = 1.5;
double q = -(42080895.0)/(14941154.0);
double r = (727668333.0)/(149411540.0);



//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonBIZ(double tau, double p){
  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);


  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;


  return epsilon_k0 + (tau0/sigma0)*(p - p_k0);
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double egsBIZ(double tau, double epsilon, double* pP, double* pC){

  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);

  double f2 = 2.0*q + 6.0*r*x + 2.0*s*f1;
  f2=f2/(1.0-s*x);

  double f3 = 6.0*r + 3.0*s*f2;
  f3=f3/(1.0-s*x);

  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

  double fac=-K0*rho0*(1.0+x)*(1.0+x)*(1.0+x)/2.0;
  double pd_k0= 2.0*(1.0+3.0*x)*f0;
  pd_k0+= 2.0*x*(2.0+3.0*x)*f1;
  pd_k0+= x*x*(1.0+x)*f2;
  pd_k0=pd_k0*fac;

  double p= p_k0 + (sigma0/tau0)*(epsilon - epsilon_k0);
 
  // resultats
  *pP= p;

  *pC= tau*sqrt( (sigma0/tau0*(p-p_k0)-pd_k0) );
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double egsBIZ_C(double tau, double epsilon, double p, double* pC){

  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);

  double f2 = 2.0*q + 6.0*r*x + 2.0*s*f1;
  f2=f2/(1.0-s*x);

  double f3 = 6.0*r + 3.0*s*f2;
  f3=f3/(1.0-s*x);

  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

  double fac=-K0*rho0*(1.0+x)*(1.0+x)*(1.0+x)/2.0;
  double pd_k0= 2.0*(1.0+3.0*x)*f0;
  pd_k0+= 2.0*x*(2.0+3.0*x)*f1;
  pd_k0+= x*x*(1.0+x)*f2;
  pd_k0=pd_k0*fac;


  *pC= tau*sqrt( (sigma0/tau0*(p-p_k0)-pd_k0) );
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
double G_BIZ(double tau, double epsilon, double* pG){

  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);

  double f2 = 2.0*q + 6.0*r*x + 2.0*s*f1;
  f2=f2/(1.0-s*x);

  double f3 = 6.0*r + 3.0*s*f2;
  f3=f3/(1.0-s*x);

  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

  double fac=-K0*rho0*(1.0+x)*(1.0+x)*(1.0+x)/2.0;
  double pd_k0= 2.0*(1.0+3.0*x)*f0;
  pd_k0+= 2.0*x*(2.0+3.0*x)*f1;
  pd_k0+= x*x*(1.0+x)*f2;
  pd_k0=pd_k0*fac;

  double p= p_k0 + (sigma0/tau0)*(epsilon - epsilon_k0);
  if (p==0){ p=1e-10; }
  // Calcul de la dérivée fondamentale
  double gamma_adia=(tau/p)*((sigma0/tau0)*(p-p_k0) - pd_k0);
  //printf("   tau=%g, p=%g, sigma0=%g, tau0=%g, p_k0=%g, pd_k0=%g, epsilon=%g, epsilon_k0=%g \n",tau,p,sigma0,tau0,p_k0,pd_k0, epsilon, epsilon_k0 );
  
  fac=K0*(1.0+x)*(1.0+x)*(1.0+x)*(1.0+x)*rho0*rho0/2;
  double pdd_k0= 12*(1.0+2*x)*f0;
  pdd_k0+= 6*(1.0 + 6*x + 6*x*x)*f1;
  pdd_k0+= 6*x*(1.0+x)*(1.0+2*x)*f2;
  pdd_k0+= x*x*(1.0+x)*(1.0+x)*f3;
  pdd_k0*=fac;
  
  fac=tau*tau/(2*gamma_adia*p);
  double G=pdd_k0;
  G+=(sigma0/tau0)*(sigma0/tau0)*(p - p_k0);
  G*=fac;

  //printf(" f_3=%g, gamma_adia=%g, pdd_k0=%g, G=%g\n",f3,gamma_adia,pdd_k0,G );
 
  *pG= G;
}


//   fonction PRESSION p + VITESSE DU SON c(tau, p)  AVEC RELAXATION PAR EVOLUTION DE LA PRESSION
// renvoie ¨Pres et Cres
double G_BIZ_cin_phase(double tau, double epsilon, double* pG,
                       int CIN_PHASE, double dx, double dt, double Ttau, 
                       double* Un, double Pn, double Vn, double Cn ){

  double x= tau0/tau - 1.0;

  double gEOS= sigma0*(1.0 - tau/tau0);

  double f0 = 1.0 + (s/3.0 - 2.0)*x + q*x*x + r*x*x*x;
  f0=f0/(1.0-s*x);

  double f1 = s/3.0 - 2.0 + 2.0*q*x + 3.0*r*x*x + s*f0;
  f1=f1/(1.0-s*x);

  double f2 = 2.0*q + 6.0*r*x + 2.0*s*f1;
  f2=f2/(1.0-s*x);

  double f3 = 6.0*r + 3.0*s*f2;
  f3=f3/(1.0-s*x);

  double epsilon_k0= eps0;
  epsilon_k0-=Cv0*T0*(1.0+ gEOS);
  epsilon_k0+=K0*tau0*x*x*f0/2.0;


  double p_k0= -Cv0*T0*sigma0*rho0;
  p_k0+= K0*x*(1.0+x)*(1.0+x)*(2.0*f0 + x*f1)/2.0;

  double fac=-K0*rho0*(1.0+x)*(1.0+x)*(1.0+x)/2.0;
  double pd_k0= 2.0*(1.0+3.0*x)*f0;
  pd_k0+= 2.0*x*(2.0+3.0*x)*f1;
  pd_k0+= x*x*(1.0+x)*f2;
  pd_k0=pd_k0*fac;

  double p= p_k0 + (sigma0/tau0)*(epsilon - epsilon_k0);
  if (p==0){ p=1e-10; }

  // Relaxation de P   optionel
  if(CIN_PHASE==3){ // dtP + rho.c2.div(u) = (-1/tau)(P-Peq)
    double Peq=p;
    double NUM= (-Cn*Cn/Vn)*(Un[1]-Un[0])/dx + Peq/Ttau;
    double DEN= 1+dt/Ttau;
    p= (Pn + dt*NUM)/DEN;
  }

  // Calcul de la dérivée fondamentale
  double gamma_adia=(tau/p)*((sigma0/tau0)*(p-p_k0) - pd_k0);
  //printf("   tau=%g, p=%g, sigma0=%g, tau0=%g, p_k0=%g, pd_k0=%g, epsilon=%g, epsilon_k0=%g \n",tau,p,sigma0,tau0,p_k0,pd_k0, epsilon, epsilon_k0 );
  
  fac=K0*(1.0+x)*(1.0+x)*(1.0+x)*(1.0+x)*rho0*rho0/2;
  double pdd_k0= 12*(1.0+2*x)*f0;
  pdd_k0+= 6*(1.0 + 6*x + 6*x*x)*f1;
  pdd_k0+= 6*x*(1.0+x)*(1.0+2*x)*f2;
  pdd_k0+= x*x*(1.0+x)*(1.0+x)*f3;
  pdd_k0*=fac;
  
  fac=tau*tau/(2*gamma_adia*p);
  double G=pdd_k0;
  G+=(sigma0/tau0)*(sigma0/tau0)*(p - p_k0);
  G*=fac;

  //printf(" f_3=%g, gamma_adia=%g, pdd_k0=%g, G=%g\n",f3,gamma_adia,pdd_k0,G );
 
  *pG= G;
}




//////////////////////////////////////////////////////////////////////
//                     fonction qui "sortent" du fichier
/////////////////////////////////////////////////////////////////////


//   fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
double epsilonEOS(int tst, double tau, double p){
  double Gamma;
  if(tst==0|| tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    return epsilonGP(Gamma,tau,p);
  }
  else if(tst==2){
    return epsilonBIZ(tau,p);
  }
  else if(tst==1){
    Gamma=5./3;
    return epsilonGP(Gamma,tau,p);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (epsilonEOS)\n");
  }
}

//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
void EGS(int tst, double tau, double epsilon, double* pP, double* pC){
  double Gamma; 
  if(tst==0 || tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    egsGP(Gamma, tau,epsilon, pP, pC);
  }
  else if(tst==2){
    egsBIZ(tau,epsilon, pP, pC);
  }
  else if(tst==1){
    Gamma=5./3;
    egsGP(Gamma, tau,epsilon, pP, pC);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (EGS)\n");
  }
}


//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
void EGS_cin_phase(int tst, double tau, double epsilon, double* pP, double* pC, 
                   int CIN_PHASE, double dx, double dt, double Ttau, 
                   double* Un, double Pn, double Vn, double Cn ){
  double Gamma; 
  if(tst==0 || tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    egsGP(Gamma, tau,epsilon, pP, pC);
  }
  else if(tst==2){
    egsBIZ(tau,epsilon, pP, pC); // pression et célérité à l'équilibre
    if(CIN_PHASE==3){ // dtP + rho.c2.div(u) = (-1/tau)(P-Peq)
      double Peq=*pP;
      double NUM= (-Cn*Cn/Vn)*(Un[1]-Un[0])/dx + Peq/Ttau;
      double DEN= 1+dt/Ttau;
      *pP= (Pn + dt*NUM)/DEN;

      // calcul de c
      egsBIZ_C(tau, epsilon, *pP, pC);
    }
  }
  else if(tst==1){
    Gamma=5./3;
    egsGP(Gamma, tau,epsilon, pP, pC);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (EGS)\n");
  }
}


//   fonction PRESSION p + VITESSE DU SON c(tau, p)
// renvoie ¨Pres et Cres
void G_EGS(int tst, double tau, double epsilon, double* pG){
  double Gamma; 
  if(tst==0 || tst==4 || tst==3 || tst==5 || tst==6){
    Gamma=1.4;
    G_GP(Gamma, tau, epsilon, pG);
  }
  else if(tst==2){
    G_BIZ(tau, epsilon, pG);
  }
  else if(tst==1){
    Gamma=5./3;
    G_GP(Gamma, tau, epsilon, pG);
  }
  else{
    printf("Erreur de choix de l'equation d'etat (EGS)\n");
  }
}