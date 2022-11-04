#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrice.h"
#include <gsl/gsl_multimin.h>

// sh  recalage_coeff.sh

//########  Recalage du point triple grace à dP/dT   #########//

// Nouveau coeff
/*
// Phase 1
double K01=53.597e9;
double N01=4.996;
double gamma01=2.27;
double Cv1=210.;
double theta01=250.;
double T01=300.;
double P01=0.;
double rho01=7287.;
double v01=1./7287.;
double E01=0.;
double Sr1=-38.2875;


// Phase 2
double K02=45.505e9;
double N02=5.239;
double gamma02=1.96;
double Cv2=210.;
double theta02=300.;
double T02=300.;
double P02=9.4e9;
double Deltav12=-2.247e-6;
double dPdT12=-2.13e7;
double Dv12=-2.247e-6;
 

// Phase 3
double K03=42.7e9;
double N03=5.;
double gamma03=2.25;
double Cv3=200.;
double theta03=350.;
double T03=505.;
double P03=0.;
double Deltav13=3.87e-6;
double dPdT13=3.125e7;
double Dv13=3.87e-6;
*/

// Ancien coeff
//##################################################################################
// Phase 1

  double K01=5.473e10;
  double N01=5.75;
  double gamma01=2.27;
  double Cv1=210.;
  double theta01=250.;
  double T01=300.;
  double P01=0.;

  double rho01=7287.;
  double v01=1./7287.;
  double E01=0.; 
  double Sr1=-38.2875;
   

  // Phase 2
  double K02=5.9e10;
  double N02=4.8;
  double gamma02=1.96;
  double Cv2=210.;
  double theta02=300.;
  double T02=300.;
  double P02=9.4e9;

  double Dv12=-1.9e-6;
  double dPdT12=-1.7e7;


  // Phase 3
  double K03=4.8e10;
  double N03=5.;
  double gamma03=2.25;
  double Cv3=200.;
  double theta03=350.;
  double T03=505.;
  double P03=0.;
  double Dv13=4e-6;
  double dPdT13=3.125e7;

//##################################################################################



double ur=1.2;


// Point triple de référence 
double PTref=3e9, TTref=580;


// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

// fonction thermo
double fx(double rho0, double V);
double THETA(double rho0, double V, int phase);// Es
double fEs(double rho0,double E0,double V, int phase);
double fPs(double rho0, double V, int phase);
double u(double rho0,double E0, double V, double E, int phase);
double fS(double rho0,double E0, double Sr, double V, double E, int phase);
double fP(double rho0,double E0, double V, double E, int phase);
double fT(double rho0,double E0, double V, double E, int phase);
double Es_prime(double rho0, double V, int phase);
double Ps_prime(double rho0, double V, int phase);
double dPdE(double rho0, int phase);
double dPdV(double rho0, double V, int phase);
double dTdE(int phase);
double dTdV(double rho0, double V, int phase);
double d2PdV2(double rho0, double V, int phase);
double d2TdV2(double rho0, double V, int phase);
void coeff(int phase, double* K0, double* N0, double* gamma0, double* Cvr, double* theta0, double* ur0);
double Ps_prime2(double rho0, double V, int phase);
double Es_prime2(double rho0, double V, int phase);
double THETA_prime(double rho0, double V, int phase);
double THETA_prime2(double rho0, double V, int phase);
double Ps_prime_rho0(double rho0, double V, int phase);
double THETA_prime_rho0(double rho0, double V, int phase);
double h(double rho0, double V, int phase);
double g(double rho0, double E0, double V, double E, int phase);


int newton_coeff(double P, double T, double V, double E, double S, double* pcoeff, int phase);
int init_coeff(double *prho0_gamma, double* pE0_gamma, double *prho0_liq, double* pE0_liq);
int print_double(int phaseB, double P0, double T0);
void Newton_double(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseA, int phaseB, double* pcoeffA, double* pcoeffB, double Va, double Ea, double Vb, double Eb, double* VA, double* EA, double* VB, double* EB);
void algo_double(int Nmax, double epsilon, int phaseA, int phaseB, double P, double* pcoeffA, double* pcoeffB, double Va0, double Ea0, double Vb0, double Eb0, double* pVa, double* pEa, double* pVb, double* pEb);
void point_double(int phaseA, int phaseB, double P, double* pcoeffA, double* pcoeffB, double* pT, double* pVa, double* pEa, double* pVb, double* pEb);
int fVE(int phase, double rho0, double E0, double P, double T, double* pV, double* pE);
void coeff(int phase, double* K0, double* N0, double* gamma0, double* Cvr, double* theta0, double* ur0);

// Optimisation du point triple
int verif_coeff();
int min_simplex_triple_gsl();
void Newton_triple(double* pcoeffgamma, double* pcoeffliq, double* A, double* InvA, double* Delta, double* dX, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC);
void triple(double* pcoeffgamma, double* pcoeffliq, int Nmax, double epsilon, double* Va0, double* Ea0, double* Vb0, double* Eb0, double* Vc0, double* Ec0);
void point_triple(double* pcoeffgamma, double* pcoeffliq, double* Pt, double* Tt, double* Vat, double* Eat, double* Vbt, double* Ebt, double* Vct, double* Ect);
int print_triple();


// hugoniot
int hugoniot(double*, double*, double*,FILE*);
int lancement_hugoniot();

int main(){
  int err;
  
  // calcul d'un point de la ligne de changement de phase
  /*
  int phaseB=2;
  double P0=9.4e9, T0=300;
  err=print_double(phaseB, P0, T0);
  */

  //err=print_triple(); if(err){return err;}
  
  //err= verif_coeff();
  
  //err= min_simplex_triple_gsl();
  
  err= lancement_hugoniot();


  return 0;
}


/*   Fonction pour obtenir les coeff (rho0, E0)
  en resolvant le système P(V,E;rho0,E0)=P et T(V,E;rho0,E0)=T  par un Newton
*/
int newton_coeff(double P, double T, double V, double E, double S, double* pcoeff, int phase){
  int Nmax=1e2;
  double epsilon=1e-15;
  int err; 
  int n=0;
  double f, fprime;
  double rho0, E0, Sr, drho0;
  
  double K0, N0, gamma0, Cvr, theta0, ur0;
  coeff(phase, &K0, &N0, &gamma0, &Cvr, &theta0, &ur0);
  
  if(phase==2){
    rho0=7300; E0=43000;
  }
  else if(phase==3){
    rho0=7036; E0=88480;
  }
  else{ printf("err phase=%d (newton_coeff/recalage_coeff_dPdT.c)\n",phase ); return 1;  }


  double critere=1.;
  
  // Newton sur E
  while(critere>epsilon && n<Nmax){
    //printf("n=%d\n",n );
    f=fPs(rho0, V, phase) - gamma0*rho0*Cvr*(ur*THETA(rho0, V, phase) - T) - P;
    fprime=Ps_prime_rho0(rho0,V,phase) - Cvr*gamma0*(ur*THETA(rho0, V, phase) - T) - Cvr*ur*gamma0*rho0*THETA_prime_rho0(rho0, V, phase);
    
    //printf("  f=%g  fprime=%g\n",f,fprime );

    drho0=f/fprime;
    rho0-=drho0;
    critere=fabs(drho0/rho0);
    //critere=fabs(f);
    n++;
    
    //E0 = Cvr*(ur*THETA(rho0,V,phase) - T) + E - h(rho0,V,phase);

    //printf(" n=%d rho0=%.15g  E0=%.15g  critere=%.3g f=%.3g  fprime=%.3g\n",n,rho0,E0,critere,f,fprime );
  }

  E0 = Cvr*(ur*THETA(rho0,V,phase) - T) + E - h(rho0,V,phase);
  Sr = S - g(rho0,E0,V,E,phase);

  pcoeff[0]=rho0; pcoeff[1]=E0; pcoeff[2]=Sr;
  
  // verif
  //printf("   newton : fP(..) - P0 =%g  fT(..) - T0=%g \n",fP(rho0,E0,V,E,2)-P,fT(rho0,E0,V,E,2)-T );

  return 0;
}

// Calcul des coefficients (rho0,E0,Sr) à partir du jeu de donnée (dPdT12,Dv12,dPdT13,Dv13) fixé
/*
   permet not. de retrouver les valeurs de Heuzé
*/
int verif_coeff(){
  int err; 
  
  // subdivision des coeffs
  double dPdT12=-2.13e7, dPdT13=3.125e7;
  double Dv12=-2.247e-6, Dv13=3.87e-6;
  

  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double P0gamma=9.4e9, T0gamma=300;
  double V_beta_gamma,E_beta_gamma,S_beta_gamma;
  err=fVE(1,rho01,E01,P0gamma,T0gamma,&V_beta_gamma,&E_beta_gamma); if(err){printf("err print_double\n"); return err;}
  S_beta_gamma=fS(rho01, E01, Sr1, V_beta_gamma, E_beta_gamma, 1);

  double P0liq=0, T0liq=505;
  double V_beta_liq,E_beta_liq,S_beta_liq;
  err=fVE(1,rho01,E01,P0liq,T0liq,&V_beta_liq,&E_beta_liq); if(err){printf("err print_double\n"); return err;}
  S_beta_liq=fS(rho01, E01, Sr1, V_beta_liq, E_beta_liq, 1);

  // point de la phase gamma et liquide
  double V_gamma,E_gamma,S_gamma;
  double V_liq,E_liq,S_liq;
 
  // coeff
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffgamma=malloc(3*sizeof(double));
  double* pcoeffliq=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;


  // Calcul du point (P0,T0) de la phase gamma (Vgamma,Egamma) et de la phase liq (Vliq,Eliq)    ! Clapeyron !
  V_gamma=V_beta_gamma + Dv12;
  S_gamma=S_beta_gamma + dPdT12*Dv12;
  E_gamma=E_beta_gamma + T0gamma*dPdT12*Dv12 - P0gamma*Dv12;

  V_liq=V_beta_liq + Dv13;
  S_liq=S_beta_liq + dPdT13*Dv13;
  E_liq=E_beta_liq + T0liq*dPdT13*Dv13 - P0liq*Dv13;
  
  // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
  err= newton_coeff(P0gamma, T0gamma, V_gamma, E_gamma, S_gamma, pcoeffgamma, 2); if(err){ printf("err (print_double)\n"); return err;}
  err= newton_coeff(P0liq, T0liq, V_liq, E_liq, S_liq, pcoeffliq, 3);             if(err){ printf("err (print_double)\n"); return err;}
  
  printf("Gamma  rho0=%.15g E0=%.15g Sr=%.15g \n",pcoeffgamma[0],pcoeffgamma[1],pcoeffgamma[2]);
  printf("Liqui  rho0=%.15g E0=%.15g Sr=%.15g \n",pcoeffliq[0],pcoeffliq[1],pcoeffliq[2]);

  return 0;
}



// Calcul du point de changement de phase entre la phase beta et B à la pression P0
/*  
     on fait varier le coeff dPdT et Dv de la phase B (gamma ou liquide)
     desquels on déduit les coeffs (rho0,E0,Sr) puis par un Newton le pt de chgmt de phase
     on imprime T0, (Va,Ea) et (Vb,Eb)
*/
int print_double(int phaseB, double P0, double T0){
  int err;
  
  // ouverture du fichier
  FILE *fdouble;
  if((fdouble = fopen("fichiers/PARAMS_double.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_double.txt \n"); return 1;}

  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double V_beta,E_beta,S_beta;
  err=fVE(1,rho01,E01,P0,T0,&V_beta,&E_beta); if(err){printf("err print_double\n"); return err;}
  S_beta=fS(rho01, E01, Sr1, V_beta, E_beta, 1);

  // point de la phase B
  double V_b,E_b,S_b;

  // coeff 
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffB=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;

  // Resultat Newton double
  double TN, VbetaN, EbetaN, VbN, EbN;

  //
  double Pbeta,Tbeta, Pb,Tb;
  
  // subdivision des coeffs
  double dPdT, Dv;
  double dPdT_min, dPdT_max;
  double Dv_min, Dv_max;
  int ndPdT=20, nDv=20;

  dPdT_min=-2.10e7; dPdT_max=-2.16e7;
  Dv_min=-2.230e-6;  Dv_max=-2.250e-6;

  double hdPdT=(dPdT_max-dPdT_min)/(ndPdT-1);
  double hDv=(Dv_max-Dv_min)/(nDv-1);
  

  // boucle sur les valeurs de dPdT12 et DeltaV
  // on calcule le point de changement de phase
  for(int i=0; i<ndPdT; i++){
    dPdT=dPdT_min + i*hdPdT;
    //printf("i=%d\n",i);

    for(int j=0; j<nDv; j++){
      printf("(i,j)=(%d,%d) | ",i,j);
      
      Dv=Dv_min + j*hDv;
      
      // Calcul du point (P0,T0) de la phase B (Vb,Eb) et S_b   ! Clapeyron !
      V_b=V_beta + Dv;
      S_b=S_beta + dPdT*Dv;
      E_b=E_beta + T0*dPdT*Dv - P0*Dv;
      

      // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
      err= newton_coeff(P0, T0, V_b, E_b, S_b, pcoeffB, phaseB); if(err){ printf("err (print_double)\n"); return err;}
            
      // Calcul de T0 (Newton)
      point_double(1, phaseB, P0, pcoeffbeta, pcoeffB, &TN, &VbetaN, &EbetaN, &VbN, &EbN);

      Pbeta=fP(rho01,E01,V_beta,E_beta,1);     Tbeta=fT(rho01,E01,V_beta,E_beta,1);
      Pb=fP(pcoeffB[0],pcoeffB[1],V_b,E_b,2);  Tb=fT(pcoeffB[0],pcoeffB[1],V_b,E_b,2);
      
      printf("TN-T0=%.15g \n",TN-T0);

      /*
      printf("dPdT=%.15g Dv=%.15g TN=%.15g \n",dPdT, Dv, TN);
      printf("       Pbeta=%.5g Pb=%.5g | Tbeta=%.5g Tb=%.5g\n",Pbeta*1e-9,Pb*1e-9,Tbeta,Tb);
      
      printf("       rho01=%g E01=%g Sr1=%g\n",rho01,E01,Sr1 );
      printf("       rho0b=%g E0b=%g Srb=%g\n",pcoeffB[0],pcoeffB[1],pcoeffB[2] );

      printf("       (V_b-VbN)=%g  (E_b-EbN)=%g\n",V_b-VbN,E_b-EbN );
      printf("\n"); 
      */

      // ecriture dans le fichier texte
      fprintf(fdouble, "%.15g %.15g %.15g \n",dPdT, Dv, TN);


    }
  }


  fclose(fdouble);
  return 0;
}



// Calcul du point triple en fonction de (dPdT12, dPdT13, Dv12, Dv13)
/*  
    on fait varier les coeffs dPdT et Dv des phases gamma et liquide (sur un quadrillage avec 2 boucles)
    desquels on déduit les coeffs (rho0,E0,Sr) 
    puis par un Newton on calcule le pt triple
    on imprime PT, TT, (VaT,EaT), (VbT,EbT) et (VcT,EcT)
*/
int print_triple(){
  int err;
  double distance;

  // ouverture du fichier
  FILE *ftriple;
  if((ftriple = fopen("fichiers/PARAMS_triple.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_triple.txt \n"); return 1;}
  FILE *fdPdT12;
  if((fdPdT12 = fopen("fichiers/dPdT12.txt", "w+")) == NULL){printf("erreur ouverture fichier fdPdT12.txt \n"); return 1;}
  FILE *fdPdT13;
  if((fdPdT13 = fopen("fichiers/dPdT13.txt", "w+")) == NULL){printf("erreur ouverture fichier fdPdT13.txt \n"); return 1;}
  

  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double P0gamma=P02, T0gamma=T02;
  double V_beta_gamma,E_beta_gamma,S_beta_gamma;
  err=fVE(1,rho01,E01,P0gamma,T0gamma,&V_beta_gamma,&E_beta_gamma); if(err){printf("err print_double\n"); return err;}
  S_beta_gamma=fS(rho01, E01, Sr1, V_beta_gamma, E_beta_gamma, 1);

  double P0liq=P03, T0liq=T03;
  double V_beta_liq,E_beta_liq,S_beta_liq;
  err=fVE(1,rho01,E01,P0liq,T0liq,&V_beta_liq,&E_beta_liq); if(err){printf("err print_double\n"); return err;}
  S_beta_liq=fS(rho01, E01, Sr1, V_beta_liq, E_beta_liq, 1);

  // point de la phase gamma et liquide
  double V_gamma,E_gamma,S_gamma;
  double V_liq,E_liq,S_liq;
 
  // coeff
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffgamma=malloc(3*sizeof(double));
  double* pcoeffliq=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;

  // Resultat Newton double
  double TN, VbetaN, EbetaN, VbN, EbN;

  //
  double Pbeta, Tbeta, Pb, Tb;

  // Point triple
  double VaT, VbT, VcT;
  double EaT, EbT, EcT;
  double PT, TT;
  
  // subdivision des coeffs
  double dPdT12, dPdT13;
  double dPdT12_min, dPdT12_max;
  double dPdT13_min, dPdT13_max;
  
  int ndPdT12=400, ndPdT13=400;

  /*
  dPdT12_min=-2.30e7; dPdT12_max=-2.10e7;
  dPdT13_min=3.00e7; dPdT13_max=3.40e7;
 */

  dPdT12_min=-2.40e7; dPdT12_max=-2.00e7;
  dPdT13_min=2.2e7; dPdT13_max=4.5e7;

  double hdPdT12=(dPdT12_max-dPdT12_min)/(ndPdT12-1);
  double hdPdT13=(dPdT13_max-dPdT13_min)/(ndPdT13-1);

  printf("hdPdT12=%g hdPdT13=%g\n",hdPdT12,hdPdT13 );
  
  //double hDv=(Dv_max-Dv_min)/(nDv-1);


  // impression des coeff variable
  for(int i=0; i<ndPdT12; i++){
    dPdT12=dPdT12_min + i*hdPdT12;
    fprintf(fdPdT12, "%.15g\n",dPdT12 );
  }
  for(int i=0; i<ndPdT13; i++){
    dPdT13=dPdT13_min + i*hdPdT13;
    fprintf(fdPdT13, "%.15g\n",dPdT13 );
  }
  fclose(fdPdT12);
  fclose(fdPdT13);

  
  double min=1e7;
  double mindPdT12=1e9;
  double mindPdT13=1e9;
  double minPT=100e9;
  double minTT=1e5;


  // boucle sur les valeurs de dPdT12 et DeltaV
  // on calcule le point de changement de phase
  for(int i=0; i<ndPdT12; i++){
    dPdT12=dPdT12_min + i*hdPdT12;
    printf("i=%d\n",i);

    for(int j=0; j<ndPdT13; j++){
      //printf("(i,j)=(%d,%d) | ",i,j);
      
      dPdT13=dPdT13_min + j*hdPdT13;
      
      // Calcul du point (P0,T0) de la phase gamma (Vgamma,Egamma) et de la phase liq (Vliq,Eliq)    ! Clapeyron !
      V_gamma=V_beta_gamma + Dv12;
      S_gamma=S_beta_gamma + dPdT12*Dv12;
      E_gamma=E_beta_gamma + T0gamma*dPdT12*Dv12 - P0gamma*Dv12;

      V_liq=V_beta_liq + Dv13;
      S_liq=S_beta_liq + dPdT13*Dv13;
      E_liq=E_beta_liq + T0liq*dPdT13*Dv13 - P0liq*Dv13;
      
      // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
      err= newton_coeff(P0gamma, T0gamma, V_gamma, E_gamma, S_gamma, pcoeffgamma, 2); if(err){ printf("err (print_double)\n"); return err;}
      err= newton_coeff(P0liq, T0liq, V_liq, E_liq, S_liq, pcoeffliq, 3);             if(err){ printf("err (print_double)\n"); return err;}
      
      // Calcul du point triple  (Newton)
      point_triple(pcoeffgamma, pcoeffliq, &PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);

      
      // ecriture dans le fichier texte
      distance=fabs( (PT-PTref)/PTref ) + fabs( (TT-TTref)/TTref );
      fprintf(ftriple, "%.15g %.15g %.15g \n",PT, TT, distance);

      //printf("  PT=%.15g  TT=%.15g distance=%g\n",PT*1e-9, TT, distance);

      if(distance<min){
        min=distance;
        mindPdT12=dPdT12;
        mindPdT13=dPdT13;
        minPT=PT; 
        minTT=TT;
      }


    }
  }
  
  printf(" minimum=%g dPdT12=%g dPdT13=%g PT=%g  TT=%g\n",min,mindPdT12,mindPdT13,minPT*1e-9,minTT );
  
  free(pcoeffbeta);
  free(pcoeffgamma);
  free(pcoeffliq);
  fclose(ftriple);
  return 0;
}


void Newton_triple(double* pcoeffgamma, double* pcoeffliq, double* A, double* InvA, double* Delta, double* dX, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC){
  int dim=6;

  double rho0gamma=pcoeffgamma[0], E0gamma=pcoeffgamma[1], Srgamma=pcoeffgamma[2];
  double rho0liq=pcoeffliq[0],     E0liq=pcoeffliq[1],     Srliq=pcoeffliq[2];

  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  double Pc, Tc, Sc, Gc, dPvc, dPec, dTvc, dTec;

  // phase beta
  Pa=fP(rho01,E01,Va,Ea, 1); Ta=fT(rho01,E01,Va,Ea, 1);
  Sa=fS(rho01,E01,Sr1,Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(rho01,Va,1); dPea=dPdE(rho01,1);
  dTva=dTdV(rho01,Va,1); dTea=dTdE(1);
  //phase gamma
  Pb=fP(rho0gamma,E0gamma,Vb,Eb, 2); Tb=fT(rho0gamma,E0gamma,Vb,Eb, 2);
  Sb=fS(rho0gamma,E0gamma,Srgamma,Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(rho0gamma,Vb,2); dPeb=dPdE(rho0gamma,2);
  dTvb=dTdV(rho0gamma,Vb,2); dTeb=dTdE(2);
  // phase liquide
  Pc=fP(rho0liq,E0liq,Vc,Ec, 3); Tc=fT(rho0liq,E0liq,Vc,Ec, 3);
  Sc=fS(rho0liq,E0liq,Srliq,Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
  dPvc=dPdV(rho0liq,Vc,3); dPec=dPdE(rho0liq,3);
  dTvc=dTdV(rho0liq,Vc,3); dTec=dTdE(3);

  double dGva, dGea, dGvb, dGeb, dGvc ,dGec;
  dGva=Va*dPva - Sa*dTva; 
  dGea=Va*dPea - Sa*dTea;

  dGvb=Vb*dPvb - Sb*dTvb; 
  dGeb=Vb*dPeb - Sb*dTeb;

  dGvc=Vc*dPvc - Sc*dTvc; 
  dGec=Vc*dPec - Sc*dTec;


  Delta[0]=Pb-Pa; Delta[1]=Pc-Pa;
  Delta[2]=Tb-Ta; Delta[3]=Tc-Ta;
  Delta[4]=Gb-Ga; Delta[5]=Gc-Ga; 
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPva;  A[dim*0+1]=dPea; A[dim*0+2]=-dPvb;  A[dim*0+3]=-dPeb;  A[dim*0+4]=0.;     A[dim*0+5]=0.;  

  A[dim*1+0]=dPva;  A[dim*1+1]=dPea; A[dim*1+2]=0.;     A[dim*1+3]=0.;     A[dim*1+4]=-dPvc;  A[dim*1+5]=-dPec; 

  A[dim*2+0]=dTva;  A[dim*2+1]=dTea; A[dim*2+2]=-dTvb;  A[dim*2+3]=-dTeb;  A[dim*2+4]=0.;     A[dim*2+5]=0.; 

  A[dim*3+0]=dTva;  A[dim*3+1]=dTea; A[dim*3+2]=0.;     A[dim*3+3]=0.;     A[dim*3+4]=-dTvc;  A[dim*3+5]=-dTec; 

  A[dim*4+0]=dGva;  A[dim*4+1]=dGea; A[dim*4+2]=-dGvb;  A[dim*4+3]=-dGeb;  A[dim*4+4]=0.;     A[dim*4+5]=0.; 

  A[dim*5+0]=dGva;  A[dim*5+1]=dGea; A[dim*5+2]=0.;     A[dim*5+3]=0.;     A[dim*5+4]=-dGvc;  A[dim*5+5]=-dGec;
  

  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if( fabs(det_a)<1e-10 ){printf("A non inversible : detA= %.8g\n", det_a);}
  
  /*
  printf("   InvA = \n");
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  
  
  // Mutiplcation matricielle de A et InvA
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      res=0;
      for(int k=0; k<6; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */

 
  // Mutiplcation matricielle
  for(int i=0; i<6; i++){
    dX[i]=0;
    for(int j=0; j<6; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
  

  double correc=0.5;

  *VA=Va+correc*dX[0];
  *EA=Ea+correc*dX[1];
  *VB=Vb+correc*dX[2];
  *EB=Eb+correc*dX[3];
  *VC=Vc+correc*dX[4];
  *EC=Ec+correc*dX[5];
}


// *** Calcul le point triple en lancant le newton  ****
/*
   Renvoie les 3 couples (V,E) du triangle
*/
void triple(double* pcoeffgamma, double* pcoeffliq, int Nmax, double epsilon, double* Va0, double* Ea0, double* Vb0, double* Eb0, double* Vc0, double* Ec0){
  int n=0;
  int dim=6;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta = (double *) malloc(sizeof(double) * (dim));
  double* dX = (double *) malloc(sizeof(double) * (dim));

  double VA, EA, VB, EB, VC, EC;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;
  
  
  double rho0gamma=pcoeffgamma[0], E0gamma=pcoeffgamma[1], Srgamma=pcoeffgamma[2];
  double rho0liq=pcoeffliq[0],     E0liq=pcoeffliq[1],     Srliq=pcoeffliq[2];
  
  
  // phase beta
  Ea=*Ea0; Va=*Va0;
  Pa=fP(rho01,E01,Va,Ea, 1); Ta=fT(rho01,E01,Va,Ea, 1);
  Sa=fS(rho01,E01,Sr1,Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase gamma
  Eb=*Eb0; Vb=*Vb0;
  Pb=fP(rho0gamma,E0gamma,Vb,Eb, 2); Tb=fT(rho0gamma,E0gamma,Vb,Eb, 2);
  Sb=fS(rho0gamma,E0gamma,Srgamma,Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase liquide
  Ec=*Ec0; Vc=*Vc0;
  Pc=fP(rho0liq,E0liq,Vc,Ec, 3); Tc=fT(rho0liq,E0liq,Vc,Ec, 3);
  Sc=fS(rho0liq,E0liq,Srliq,Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
  
  /*
  printf("  (triple) init Pa=%g Pb=%g Pc=%g\n",Pa,Pb,Pc);
  printf("                Ta=%g Tb=%g Tc=%g\n",Ta,Tb,Tc);
  */

  double critere=1.;
  
  /*
  printf("init" );
  printf("  P_beta=%.15g P_gamma=%.15g P_liq=%.15g (GPa) ",fP(dPdT12,dPdT13,Va,Ea, 1)*1e-9,fP(dPdT12,dPdT13,Vb,Eb, 2)*1e-9,fP(dPdT12,dPdT13,Vc,Ec, 3)*1e-9 );
  printf("| T_beta=%.15g  T_gamma=%.15g T_liq=%.15g (K)\n",fT(dPdT12,dPdT13,Va,Ea, 1),fT(dPdT12,dPdT13,Vb,Eb, 2),fT(dPdT12,dPdT13,Vc,Ec, 3) );
  */

  while(critere>=epsilon && n<Nmax){
    Newton_triple(pcoeffgamma, pcoeffliq, A, InvA, Delta, dX, Va, Ea, Vb, Eb, Vc, Ec, &VA, &EA, &VB, &EB, &VC, &EC);

    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    Vc=VC; Ec=EC;

    //printf("   Va=%g Vb=%g Vc=%g\n",Va,Vb,Vc);
    //printf("   Ea=%g Eb=%g Ec=%g\n",Ea,Eb,Ec);
    
    // phase beta
    Pa=fP(rho01,E01,Va,Ea, 1); Ta=fT(rho01,E01,Va,Ea, 1);
    Sa=fS(rho01,E01,Sr1,Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
    //phase gamma
    Pb=fP(rho0gamma,E0gamma,Vb,Eb, 2); Tb=fT(rho0gamma,E0gamma,Vb,Eb, 2);
    Sb=fS(rho0gamma,E0gamma,Srgamma,Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
    // phase liquide
    Pc=fP(rho0liq,E0liq,Vc,Ec, 3); Tc=fT(rho0liq,E0liq,Vc,Ec, 3);
    Sc=fS(rho0liq,E0liq,Srliq,Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
    
    

    //critere=fabs(Ga-Gb)+fabs(Ga-Gc);
    critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb) + fabs(dX[4]/Vc) + fabs(dX[5]/Ec);
    /*
    printf("      n=%d",n );
    printf(" critere= %.3g",critere );
    printf("   Pa=%g Pb=%g Pc=%g",Pa*1e-9,Pb*1e-9,Pc*1e-9);
    printf("   Ta=%g Tb=%g Tc=%g \n",Ta,Tb,Tc);
    */
    n++;
  }

  *Va0=Va; *Ea0=Ea; 
  *Vb0=Vb; *Eb0=Eb; 
  *Vc0=Vc; *Ec0=Ec;

  free(A); free(InvA);
  free(Delta);
  free(dX);
}


/*  ****  Calcul du Point TRIPLE à partir des coeff (rho0,E0,Sr) des phases gamma et liquide  ****
   utilisé dans  main_calfront.c
*/
void point_triple(double* pcoeffgamma, double* pcoeffliq, double* Pt, double* Tt, double* Vat, double* Eat, double* Vbt, double* Ebt, double* Vct, double* Ect){
  int Nmax=100;
  double epsilon=1e-10;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;
  double Pc, Tc, Sc, Gc, Ec, Vc;

  //double Va0=129e-6, Ea0=100e3;
  
  double Va0=125e-6, Ea0=100e3;

  double Vb0=0.98*Va0;
  double Eb0=1.5*Ea0;
  
  double Vc0=1.02*Va0;
  double Ec0=2*Ea0;

  /*
  double Va0,Ea0,Vb0,Eb0,Vc0,Ec0;
  Va0=1./7287;
  Vb0=7301.02; Vc0=7036.77;

  Ea0=0;
  Eb0=43080.8; Ec0=88479.8;
  */

  triple(pcoeffgamma, pcoeffliq, Nmax, epsilon, &Va0, &Ea0, &Vb0, &Eb0, &Vc0, &Ec0);
  
  double rho0gamma=pcoeffgamma[0], E0gamma=pcoeffgamma[1], Srgamma=pcoeffgamma[2];
  double rho0liq=pcoeffliq[0],     E0liq=pcoeffliq[1],     Srliq=pcoeffliq[2];

  // phase beta
  Ea=Ea0; Va=Va0;
  Pa=fP(rho01,E01,Va,Ea, 1); Ta=fT(rho01,E01,Va,Ea, 1);
  Sa=fS(rho01,E01,Sr1,Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase gamma
  Eb=Eb0; Vb=Vb0;
  Pb=fP(rho0gamma,E0gamma,Vb,Eb, 2); Tb=fT(rho0gamma,E0gamma,Vb,Eb, 2);
  Sb=fS(rho0gamma,E0gamma,Srgamma,Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase liquide
  Ec=Ec0; Vc=Vc0;
  Pc=fP(rho0liq,E0liq,Vc,Ec, 3); Tc=fT(rho0liq,E0liq,Vc,Ec, 3);
  Sc=fS(rho0liq,E0liq,Srliq,Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;

  //sortie des donnes du point triple
  *Vat=Va; *Eat=Ea; 
  *Vbt=Vb; *Ebt=Eb; 
  *Vct=Vc; *Ect=Ec; 
  *Pt=(Pa+Pb+Pc)/3.0;
  *Tt=(Ta+Tb+Tc)/3.0;

  /*
  printf("Resultats :\n");
  printf(" Ga= %.5lf, Gb= %.5lf,  Gc= %.5lf\n",Ga, Gb, Gc );
  printf(" Pa= %.5lf, Pb= %.5lf,  Pc= %.5lf  (GPa)\n",Pa*1e-9, Pb*1e-9, Pc*1e-9 );
  printf(" Ta= %.5lf, Tb= %.5lf,  Tc= %.5lf\n",Ta, Tb, Tc );

  printf("Va=%.5g Ea=%.5g\n",Va,Ea);
  printf("Vb=%.5g Eb=%.5g\n",Vb,Eb);
  printf("Vc=%.5g Ec=%.5g\n",Vc,Ec);
  */

}


/*  ****  Calcul du Point TRIPLE à partir des coeff (dP/dT) des phases gamma et liquide  ****
   contient le calcul des coeff (rho0,E0,Sr) correspondant (Newton)
   utilisé dans  main_calfront.c
*/
int f_triple_dPdT (double dPdT12, double dPdT13, double* pPT, double* pTT){

  int err;
  double distance;

  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double P0gamma=P02, T0gamma=T02;
  double V_beta_gamma,E_beta_gamma,S_beta_gamma;
  err=fVE(1,rho01,E01,P0gamma,T0gamma,&V_beta_gamma,&E_beta_gamma); if(err){printf("err print_double\n"); return err;}
  S_beta_gamma=fS(rho01, E01, Sr1, V_beta_gamma, E_beta_gamma, 1);

  double P0liq=P03, T0liq=T03;
  double V_beta_liq,E_beta_liq,S_beta_liq;
  err=fVE(1,rho01,E01,P0liq,T0liq,&V_beta_liq,&E_beta_liq); if(err){printf("err print_double\n"); return err;}
  S_beta_liq=fS(rho01, E01, Sr1, V_beta_liq, E_beta_liq, 1);

  // point de la phase gamma et liquide
  double V_gamma,E_gamma,S_gamma;
  double V_liq,E_liq,S_liq;
 
  // coeff
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffgamma=malloc(3*sizeof(double));
  double* pcoeffliq=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;

  // Point triple
  double VaT, VbT, VcT;
  double EaT, EbT, EcT;
  double PT, TT;
  
      
  // Calcul du point (P0,T0) de la phase gamma (Vgamma,Egamma) et de la phase liq (Vliq,Eliq)    ! Clapeyron !
  V_gamma=V_beta_gamma + Dv12;
  S_gamma=S_beta_gamma + dPdT12*Dv12;
  E_gamma=E_beta_gamma + T0gamma*dPdT12*Dv12 - P0gamma*Dv12;

  V_liq=V_beta_liq + Dv13;
  S_liq=S_beta_liq + dPdT13*Dv13;
  E_liq=E_beta_liq + T0liq*dPdT13*Dv13 - P0liq*Dv13;
      
  // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
  err= newton_coeff(P0gamma, T0gamma, V_gamma, E_gamma, S_gamma, pcoeffgamma, 2); if(err){ printf("err (print_double)\n"); return err;}
  err= newton_coeff(P0liq, T0liq, V_liq, E_liq, S_liq, pcoeffliq, 3);             if(err){ printf("err (print_double)\n"); return err;}
      
  // Calcul du point triple  (Newton)
  point_triple(pcoeffgamma, pcoeffliq, &PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);

  *pPT=PT; 
  *pTT=TT; 


  //distance=fabs( (PT-PTref)/PTref ) + fabs( (TT-TTref)/TTref );

  //params[0]=PT;
  //params[1]=TT;
 
  free(pcoeffbeta);
  free(pcoeffgamma);
  free(pcoeffliq);

  return 0;
}




////////   Méthode de minimisation sur la valeur du point triple à partir de dP/dTgamma et dP/dTliquide grâce à la méthode du simplex de Nelder-Mead

/* 
    fonction qui à partir des valeurs de dP/dTgamma et dP/dTliquide
    calcule les coeff (rho0,E0,Sr) des phases gamma et liquide (2 Newtons)
    puis calcule le point triple (Newton) 
    et renvoie l'erreur au point triple de référence
    formatée pour être utilisée avec les fonctions de la librairie GSL
*/
double f_triple_gsl (const gsl_vector *v, void *params){

  //double *p = (double *)params;

  double dPdT12 = gsl_vector_get(v, 0);
  double dPdT13 = gsl_vector_get(v, 1);

  int err;
  double distance;

  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double P0gamma=P02, T0gamma=T02;
  double V_beta_gamma,E_beta_gamma,S_beta_gamma;
  err=fVE(1,rho01,E01,P0gamma,T0gamma,&V_beta_gamma,&E_beta_gamma); if(err){printf("err print_double\n"); return err;}
  S_beta_gamma=fS(rho01, E01, Sr1, V_beta_gamma, E_beta_gamma, 1);

  double P0liq=P03, T0liq=T03;
  double V_beta_liq,E_beta_liq,S_beta_liq;
  err=fVE(1,rho01,E01,P0liq,T0liq,&V_beta_liq,&E_beta_liq); if(err){printf("err print_double\n"); return err;}
  S_beta_liq=fS(rho01, E01, Sr1, V_beta_liq, E_beta_liq, 1);

  // point de la phase gamma et liquide
  double V_gamma,E_gamma,S_gamma;
  double V_liq,E_liq,S_liq;
 
  // coeff
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffgamma=malloc(3*sizeof(double));
  double* pcoeffliq=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;

  // Point triple
  double VaT, VbT, VcT;
  double EaT, EbT, EcT;
  double PT, TT;
  
      
  // Calcul du point (P0,T0) de la phase gamma (Vgamma,Egamma) et de la phase liq (Vliq,Eliq)    ! Clapeyron !
  V_gamma=V_beta_gamma + Dv12;
  S_gamma=S_beta_gamma + dPdT12*Dv12;
  E_gamma=E_beta_gamma + T0gamma*dPdT12*Dv12 - P0gamma*Dv12;

  V_liq=V_beta_liq + Dv13;
  S_liq=S_beta_liq + dPdT13*Dv13;
  E_liq=E_beta_liq + T0liq*dPdT13*Dv13 - P0liq*Dv13;
      
  // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
  err= newton_coeff(P0gamma, T0gamma, V_gamma, E_gamma, S_gamma, pcoeffgamma, 2); if(err){ printf("err (print_double)\n"); return err;}
  err= newton_coeff(P0liq, T0liq, V_liq, E_liq, S_liq, pcoeffliq, 3);             if(err){ printf("err (print_double)\n"); return err;}
  
  printf("     rho0gamma=%.15g  E0gamma=%.15g   Srgamma=%.15g\n",pcoeffgamma[0],pcoeffgamma[1],pcoeffgamma[2] );
  printf("     rho0liq=%.15g  E0liq=%.15g   Srliq=%.15g\n",pcoeffliq[0],pcoeffliq[1],pcoeffliq[2] );

  // Calcul du point triple  (Newton)
  point_triple(pcoeffgamma, pcoeffliq, &PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);

      
  distance=fabs( (PT-PTref)/PTref ) + fabs( (TT-TTref)/TTref );

  //params[0]=PT;
  //params[1]=TT;
 
  free(pcoeffbeta);
  free(pcoeffgamma);
  free(pcoeffliq);

  return distance;
}



/* 

*/
int min_simplex_triple_gsl(){
  
  double par[2] = {0.0 , 0.0};
  double PT,TT;
  int err;

  const gsl_multimin_fminimizer_type *T =
  gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;


  double dPdT12=-2.19e7, dPdT13=3.11e7;
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, dPdT12);
  gsl_vector_set (x, 1, dPdT13);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1e6);

  /* Initialize method and iterate */
  minex_func.n = 2;
  minex_func.f = f_triple_gsl;
  minex_func.params = par;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 4e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
          printf(" Ptriple=%.10g Ttriple=%.10g\n",par[0],par[1] );
        }

      printf ("%5ld %15.15g %15.15g f() = %g size = %.3f\n",
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              s->fval, size);

      if (status == GSL_SUCCESS)
        {
          err=f_triple_dPdT (gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1), &PT, &TT);

          printf(" Ptriple=%.15g (GPa)  Ttriple=%.15g (K) \n",PT*1e-9,TT );
        }
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}


// *******  Newton sur les lignes de chgmt de phase  ********
/* 
  Calcul d'un point de la ligne de changement de phase entre de 2 phase dont on transmet les coeffs
  Valeur d'entree :
   Va,Ea, Vb,Eb, Vc,Ec

   corespodances :  1-a   2-b   3-c

  Valeur sortie
   VA,EA, VB,EB, VC,EC

*/
void Newton_double(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseA, int phaseB, double* pcoeffA, double* pcoeffB, double Va, double Ea, double Vb, double Eb, double* VA, double* EA, double* VB, double* EB){
  int dim=4;
 
  double rho0a=pcoeffA[0], E0a=pcoeffA[1], Sra=pcoeffA[2];
  double rho0b=pcoeffB[0], E0b=pcoeffB[1], Srb=pcoeffB[2];

  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  // phase A
  Pa=fP(rho0a, E0a, Va, Ea, phaseA); Ta=fT(rho0a, E0a, Va, Ea, phaseA);
  Sa=fS(rho0a, E0a, Sra, Va, Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(rho0a,Va,phaseA); dPea=dPdE(rho0a, phaseA);
  dTva=dTdV(rho0a,Va,phaseA); dTea=dTdE(phaseA);
  // phase B
  Pb=fP(rho0b, E0b, Vb,Eb, phaseB); Tb=fT(rho0b, E0b, Vb,Eb, phaseB);
  Sb=fS(rho0b, E0b, Srb, Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(rho0b, Vb,phaseB); dPeb=dPdE(rho0b, phaseB);
  dTvb=dTdV(rho0b, Vb,phaseB); dTeb=dTdE(phaseB);


  double dGva, dGea, dGvb, dGeb;
  dGva=Va*dPva - Sa*dTva; 
  dGea=Va*dPea - Sa*dTea;

  dGvb=Vb*dPvb - Sb*dTvb; 
  dGeb=Vb*dPeb - Sb*dTeb;


  Delta[0]=Pa-Pcst; Delta[1]=Pb-Pcst;
  Delta[2]=Tb-Ta;   Delta[3]=Gb-Ga;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  }
  printf("\n");
  */


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=-dPva; A[dim*0+1]=-dPea; A[dim*0+2]=0.;    A[dim*0+3]=0.;  

  A[dim*1+0]=0.;    A[dim*1+1]=0.;    A[dim*1+2]=-dPvb; A[dim*1+3]=-dPeb;      

  A[dim*2+0]=dTva;  A[dim*2+1]=dTea;  A[dim*2+2]=-dTvb; A[dim*2+3]=-dTeb;  

  A[dim*3+0]=dGva;  A[dim*3+1]=dGea;  A[dim*3+2]=-dGvb; A[dim*3+3]=-dGeb; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.5g\n", det_a);}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  /*
  // Mutiplcation matricielle de A et InvA
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Multiplication matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1./2;

  *VA=Va+correc*dX[0];
  *EA=Ea+correc*dX[1];
  *VB=Vb+correc*dX[2];
  *EB=Eb+correc*dX[3];

}

//  *******  Fonction qui calcule le points (V,E) d'un point d'une ligne de chgmt de phases
/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - [P : pression de la ligne de changement de phase
  - Va0, Ea0, Vb0, Eb0 : Valeur initiale de (Va,Ea) et (Vb,Eb)
   (typiquement les coord du point triple pour chaque phase)

Arguments de sortie :
  - T : temperature de la ligne de changement de phase
  - (Va, Ea) : coordonées pour la phase A du point (P,T)
  - (Vb, Eb) : coordonées pour la phase B du point (P,T)

*/
void algo_double(int Nmax, double epsilon, int phaseA, int phaseB, double P, double* pcoeffA, double* pcoeffB, double Va0, double Ea0, double Vb0, double Eb0, double* pVa, double* pEa, double* pVb, double* pEb){
  int n=1;
  int dim=4;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta = (double *) malloc(sizeof(double) * (dim));
  double* dX = (double *) malloc(sizeof(double) * (dim));

  double rho0a=pcoeffA[0], E0a=pcoeffA[1], Sra=pcoeffA[2];
  double rho0b=pcoeffB[0], E0b=pcoeffB[1], Srb=pcoeffB[2];

  double VA, EA, VB, EB;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;

  Ea=Ea0; Va=Va0;
  Eb=Eb0; Vb=Vb0;


  double critere=1.;
  n=1;
  while(critere>=epsilon && n<Nmax){
    Newton_double(A, InvA, Delta, dX, P, phaseA, phaseB, pcoeffA, pcoeffB, Va, Ea, Vb, Eb, &VA, &EA, &VB, &EB);
    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb);

    Pa=fP(rho0a, E0a, Va, Ea, phaseA); Ta=fT(rho0a, E0a, Va, Ea, phaseA);
    Sa=fS(rho0a, E0a, Sra, Va, Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
    Pb=fP(rho0b, E0b, Vb, Eb, phaseB); Tb=fT(rho0b, E0b, Vb, Eb, phaseB);
    Sb=fS(rho0b, E0b, Srb, Vb, Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;

    //printf(" n=%d T=%.15g  critere= %.3g ",n,fT(Va,Ea,phaseA),critere );
    //printf(" Ga-Gb=%.3g  Pa-Pb=%.3g  Ta-Tb=%.3g\n",Ga-Gb,Pa-Pb,Ta-Tb );
    n++;
  }
  *pVa=Va;
  *pEa=Ea;
  *pVb=Vb;
  *pEb=Eb;

  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
}


//  *******  Fonction qui calcule les points (Va,Ea) et (Vb,Eb) d'un point d'une ligne de chgmt de phases
// On rentre une pression et les deux phases que l'on souhaite
//   utilisé dans  main_calfront.c

void point_double(int phaseA, int phaseB, double P, double* pcoeffA, double* pcoeffB, double* pT, double* pVa, double* pEa, double* pVb, double* pEb){
  int Nmax=100;
  double epsilon=1e-10;

  double rho0a=pcoeffA[0], E0a=pcoeffA[1], Sra=pcoeffA[2];
  double rho0b=pcoeffB[0], E0b=pcoeffB[1], Srb=pcoeffB[2];

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;

  double Va0, Ea0, Vb0, Eb0;

  double Vbeta0=125e-6, Ebeta0=100e3;
  

  if(phaseA==2 || phaseB==2){
    Va0=0.98*Vbeta0;       Ea0=1.5*Ebeta0;
    Vb0=0.98*Vbeta0;  Eb0=1.5*Ebeta0;
  }  
  if(phaseA==1){
    Va0=Vbeta0; Ea0=Ebeta0;
  }
  if(phaseB==3){
    Vb0=1.02*Vbeta0;
    Eb0=2*Ebeta0;
  }


  algo_double(Nmax, epsilon, phaseA, phaseB, P, pcoeffA, pcoeffB, Va0, Ea0, Vb0, Eb0, pVa, pEa, pVb, pEb);

  
  // phase a
  Ea=*pEa; Va=*pVa;
  Pa=fP(rho0a, E0a, Va, Ea, phaseA); Ta=fT(rho0a, E0a, Va, Ea, phaseA);
  Sa=fS(rho0a, E0a, Sra, Va, Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
  // phase b
  Eb=*pEb; Vb=*pVb;
  Pb=fP(rho0b, E0b, Vb, Eb, phaseB); Tb=fT(rho0b, E0b, Vb, Eb, phaseB);
  Sb=fS(rho0b, E0b, Srb, Vb, Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;

  //sortie des donnes du point triple
  *pVa=Va; *pEa=Ea; 
  *pVb=Vb; *pEb=Eb; 
  *pT=(Ta+Tb)/2.0;


  /*
  printf("\nResultats double P=%g:\n",P);
  printf(" Ga= %g, Gb= %g,  Ga-Gb= %g\n",Ga, Gb, Ga-Gb );
  printf(" Pa= %g, Pb= %g,  Pa-Pb= %g  (GPa)\n",Pa*1e-9, Pb*1e-9, (Pa-Pb)*1e-9 );
  printf(" Ta= %g, Tb= %g,  Ta-Tb= %g (K)\n",Ta, Tb, Ta-Tb );

  printf("Va=%.5g Ea=%.5g\n",Va,Ea);
  printf("Vb=%.5g Eb=%.5g\n",Vb,Eb);
  */
}


// ***************************************************************************************************************************************************************************************************************
//                        Courbes d'Hugoniot  en fonction de (dP/dT gamma et liquide)
// ***************************************************************************************************************************************************************************************************************


/* Valeur d'entree :
    _ phaseX, phaseY
    - Vx,Ex, Vy,Ey
    - Vpts, Epts : plot d'Hogoniot
    - x : la fraction massique de phase Y est imposé (ici 0 ou 1)

    - A,InvA,Delta,dX : pour limiter les allocations 

  Valeur sortie
    - Vx,Ex, Vy,Ey
*/
int  Newton_hugniot_pt12(int cas, double* A, double* InvA, double* Delta, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* pVX, double* pEX, double* pVY, double* pEY, double* dX){
  int dim=4;

  double rho0x=pcoeffX[0], E0x=pcoeffX[1], Srx=pcoeffX[2];
  double rho0y=pcoeffY[0], E0y=pcoeffY[1], Sry=pcoeffY[2];


  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(rho0x, E0x, Vx, Ex, phaseX); Tx=fT(rho0x, E0x, Vx, Ex, phaseX);
  Sx=fS(rho0x, E0x, Srx, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(rho0x,Vx,phaseX); dPex=dPdE(rho0x, phaseX);
  dTvx=dTdV(rho0x,Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(rho0y, E0y, Vy, Ey, phaseY); Ty=fT(rho0y, E0y, Vy, Ey, phaseY);
  Sy=fS(rho0y, E0y, Sry, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(rho0y, Vy, phaseY); dPey=dPdE(rho0y, phaseY);
  dTvy=dTdV(rho0y, Vy, phaseY); dTey=dTdE(phaseY);
  
  /*
  printf("Px= %.10lf, Tx= %.10lf, Sx= %.10lf, Gx= %.10lf\n",Px,Tx,Sx,Gx );
  printf("Py= %.10lf, Ty= %.10lf, Sy= %.10lf, Gy= %.10lf\n",Py,Ty,Sy,Gy );
  */
  
  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx; 
  dGex=Vx*dPex - Sx*dTex;
  
  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  if (cas==1){      Delta[3]=2*(Ept-Ex) + (Ppt+Px)*(Vpt-Vx) ; }
  else if(cas==2){  Delta[3]=2*(Ept-Ey) + (Ppt+Py)*(Vpt-Vy) ; }
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  // Dérive des equation Hugoniot
  double dHvx, dHex;
  double dHvy, dHey;

  if(cas==1){
    dHvx= Ppt - Vpt*dPvx + Px + Vx*dPvx;
    dHex= 2.0 - Vpt*dPex + Vx*dPex;
    dHvy= 0;
    dHey= 0;
  }
  else if(cas==2){   
    dHvx= 0;
    dHex= 0;
    dHvy= Ppt - Vpt*dPvy + Py + Vy*dPvy;
    dHey= 2.0 - Vpt*dPey + Vy*dPey;
  }
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey;   
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; 
  
  A[dim*3+0]=dHvx; A[dim*3+1]=dHex; A[dim*3+2]=dHvy;  A[dim*3+3]=dHey;     
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){  printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );  }  printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Print de la mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];

  
  //printf("  Vx= %g, Ex= %g ,  Vy= %g, Ey= %g \n",*pVX,*pEX,*pVY,*pEY );

  return 0;
}


/*  Calcul du point 1 et 2 sur le courbe d'hugoniot
      point 1 : frontière entre la phase beta et la zone de mélange beta/gamma
      point 2 : sortie de la zone de mélange beta/gamma
*/
int hugoniot_pts12(int cas, int Nmax, double epsilon, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double* pVx, double* pEx, double* pVy, double* pEy, double* pP, double* pT, double* pS, double* pG){
  int n=1;
  int dim=4;
  int err;
  
  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double Vx, Ex, Vy, Ey;
  double VX, EX, VY, EY;

  double P,T,G,S;
  
  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  double rho0_x=pcoeffX[0], E0_x=pcoeffX[1], Sr_x=pcoeffX[2];
  double rho0_y=pcoeffY[0], E0_y=pcoeffY[1], Sr_y=pcoeffY[2];

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(rho0_x, E0_x, Vx, Ex, phaseX); Tx=fT(rho0_x, E0_x, Vx, Ex, phaseX);
  Sx=fS(rho0_x, E0_x, Sr_x, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(rho0_y, E0_y, Vy, Ey, phaseY); Ty=fT(rho0_y, E0_y, Vy, Ey, phaseY);
  Sy=fS(rho0_y, E0_y, Sr_y, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  
  
  double critere=1.;
  
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    err= Newton_hugniot_pt12(cas, A, InvA, Delta, phaseX, phaseY, pcoeffX, pcoeffY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, dX); 
    if(err){return 1;}    

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/Vx) +fabs(dX[1]/Ex) +fabs(dX[2]/Vy) +fabs(dX[3]/Ey);
    //critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    //printf(" --n= %d  critere= %g\n",n,critere );

    
    n++;
  }
  printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d , critere= %g\n",n,critere );return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;

  // phase X
  Px=fP(rho0_x, E0_x, Vx, Ex, phaseX); Tx=fT(rho0_x, E0_x, Vx, Ex, phaseX);
  Sx=fS(rho0_x, E0_x, Sr_x, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(rho0_y, E0_y, Vy, Ey, phaseY); Ty=fT(rho0_y, E0_y, Vy, Ey, phaseY);
  Sy=fS(rho0_y, E0_y, Sr_y, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;


  *pP=(Px+Py)/2;
  *pT=(Tx+Ty)/2;
  *pG=(Gx+Gy)/2;
  
  if (cas==1){     *pS=Sx ; }
  else if(cas==2){ *pS=Sy ; }


  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
}


int  Newton_hugniot_pt3(double* A, double* InvA, double* Delta, double* dX, double* pcoeff, double* Vpt, double* Ept, double* Ppt, double V, double E, double* pV, double* pE){
  int dim=2;

  double rho0=pcoeff[0], E0=pcoeff[1], Sr=pcoeff[2];

  // Grandeurs thermodynamiques
  double P, dPv, dPe;
  // phase X
  P=fP(rho0, E0, V, E, 2); 
  dPv=dPdV(rho0, V, 2); dPe=dPdE(rho0, 2);


  double D1=Vpt[0]*sqrt((Ppt[1]-Ppt[0])/(Vpt[0]-Vpt[1]));
  //double D1=Vpt[0]*sqrt((Ppt[1]-Ppt[0])/(Vpt[0]-Vpt[1]));
  
  /*
  printf("Px= %.10lf, Tx= %.10lf, Sx= %.10lf, Gx= %.10lf\n",Px,Tx,Sx,Gx );
  printf("Py= %.10lf, Ty= %.10lf, Sy= %.10lf, Gy= %.10lf\n",Py,Ty,Sy,Gy );
  */
  
  double u1=Vpt[0]*(Ppt[1]-Ppt[0])/D1;
  double DuV2=((D1-u1)/Vpt[1])*((D1-u1)/Vpt[1]);
  double DPDV=((Ppt[1]-Ppt[0])/(Vpt[1]-Vpt[0]));
  //printf("u1= %g\n",u1 );
  
  Delta[0]=2*(Ept[0]-E) + (Ppt[0]+P)*(Vpt[0]-V) ;
  
  //Delta[1]=Ppt[1]-P + DPDV*(V-Vpt[1]) ;
  Delta[1]=Ppt[0]-P + DPDV*(V-Vpt[0]) ;


  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]= (V-Vpt[0])*dPv + Ppt[0] + P;  /*A[dim*0+1]= 2.0 + (V- Vpt[0])*dPe;*/   A[dim*0+1]= 2.0 - (V- Vpt[0])*dPe ; 
  A[dim*1+0]= dPv-DPDV ;                    A[dim*1+1]= dPe ; 

  //A[dim*1+0]= dPv-DuV2 ;  A[dim*1+1]= dPe ;
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){  printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );  }  printf("\n");
  }
  printf("\n");
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Print de la mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *pV=V+correc*dX[0];
  *pE=E+correc*dX[1];
  
  //printf("  V= %g, E= %g \n",*pV,*pE );
  
  return 0;
}


int hugniot_pt3(int Nmax, double epsilon, double* pcoeff, double* Vpt, double* Ept, double* Ppt, double V0, double E0, double* pV, double* pE){
  int n=0;
  int dim=2;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double P;

  double V, E;
  double Vin, Ein;

  Vin=V0; Ein=E0; 

  double rho_0=pcoeff[0], E_0=pcoeff[1], S_r=pcoeff[2];
  
  // phase X
  P=Ppt[0]; 

  double P3,P0,P1,test,a;

  double critere;
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=1;
  while(critere>=epsilon && n<Nmax){
    //printf("  --n= %d\n",n );

    err= Newton_hugniot_pt3(A, InvA, Delta, dX, pcoeff, Vpt, Ept, Ppt, Vin, Ein, &V, &E);
    if(err){return 1;}    

    Vin=V; Ein=E;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/V) +fabs(dX[1]/E);
    //critere=fabs(dX[0]) +fabs(dX[1]) +fabs(dX[2]) +fabs(dX[3])+fabs(dX[4]);
    //printf(" --n= %d  critere= %g, dV/V= %g, dE/E= %g\n",n,critere,dX[0]/V,dX[1]/E );
    n++;
  }
  //printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d , critere= %g\n",n,critere );return 1;}

  *pV=V; *pE=E;
  //printf(" P3= %g, T3= %g\n", fP(V,E,2), fT(V,E,2) );


  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
}


/* Algo Newton-Raphson à fraction massique X fixé pour calculer un point de la courbe d'Hugoniot dans la zone de mélange X/Y
*/
int Newton_hugoniot_mixte_Xfixe(double* A, double* InvA, double* Delta, double* dX, double Xcst, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* VX, double* EX, double* VY, double* EY){
  int dim=4;
  
  double rho0x=pcoeffX[0], E0x=pcoeffX[1], Srx=pcoeffX[2];
  double rho0y=pcoeffY[0], E0y=pcoeffY[1], Sry=pcoeffY[2];

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(rho0x, E0x, Vx, Ex, phaseX); Tx=fT(rho0x, E0x, Vx, Ex, phaseX);
  Sx=fS(rho0x, E0x, Srx, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(rho0x,Vx,phaseX); dPex=dPdE(rho0x, phaseX);
  dTvx=dTdV(rho0x,Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(rho0y, E0y, Vy, Ey, phaseY); Ty=fT(rho0y, E0y, Vy, Ey, phaseY);
  Sy=fS(rho0y, E0y, Sry, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(rho0y, Vy, phaseY); dPey=dPdE(rho0y, phaseY);
  dTvy=dTdV(rho0y, Vy, phaseY); dTey=dTdE(phaseY);


  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx;
  dGex=Vx*dPex - Sx*dTex;

  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  double P=(Px+Py)/2;
  P=Px;
  double P2=(P+Ppt)/2;

  double x=Xcst;

  double Vmel=(1-x)*Vx + x*Vy;
  double Emel=(1-x)*Ex + x*Ey;

  Delta[0]=Py-Px;  Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  Delta[3]=Ept - Emel + P2*(Vpt - Vmel);
   
  
  //for(int j=0; j<dim; j++){
  //  printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  //}
  //printf("\n");
  


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;     A[dim*0+1]=dPex;  A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey;

  A[dim*1+0]=dTvx;     A[dim*1+1]=dTex;  A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey;

  A[dim*2+0]=dGvx;     A[dim*2+1]=dGex;  A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; 

  A[dim*3+0]=P2*(1-x) + ((Vmel-Vpt)/2)*dPvx; A[dim*3+1]=(1-x) + ((Vmel-Vpt)/2)*dPex; A[dim*3+2]=P2*x;  A[dim*3+3]=x;     
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g ",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1.;

  *VX=Vx+correc*dX[0];
  *EX=Ex+correc*dX[1];
  *VY=Vy+correc*dX[2];
  *EY=Ey+correc*dX[3];

  return 0;
}

/*  Calcul des points de l'Hugoniot dans la zone de mélange X/Y 
    à partir d'une discrétisation en fraction massique contenu dans tabX
    à N points 
*/
int Pt_hugoniot_mixte_Xfixe(int Nmax, double epsilon, int N, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double* tabX, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double* tabV, double* tabE, double* tabP, double* tabT, double* tabS, double* tabG){
  // portion de courbe 0-1
  int n=0;
  
  int dim=4;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double rho0_x=pcoeffX[0], E0_x=pcoeffX[1], Sr_x=pcoeffX[2];
  double rho0_y=pcoeffY[0], E0_y=pcoeffY[1], Sr_y=pcoeffY[2];

  double Vx,Ex,Vy,Ey;
  double VX,EX,VY,EY;

  double Gx,Gy, Sx,Sy;
  double Px,Py, Tx,Ty;
  double P,T,G,S;

  double Vmel, Emel, x;
  double critere;
  
  Vx=V0x;   Vy=V0y; 
  Ex=E0x;   Ey=E0y; 
  
  for (int i = 0; i <N; i++){
    x=tabX[i];
    Vx=V0x;   Vy=V0y; 
    Ex=E0x;   Ey=E0y; 
    critere=1.; n=0;
    while(critere>epsilon && n<Nmax){
      // Ici dim=4   // MODIFER Newton_double pour liliter les allocations
      err=Newton_hugoniot_mixte_Xfixe(A, InvA, Delta, dX, x, phaseX, phaseY, pcoeffX, pcoeffY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY);
      if(err){return 1;}
      //printf(" --n= %d , critere= %g\n",n,critere );
      Vx=VX;   Vy=VY; 
      Ex=EX;   Ey=EY; 

      critere=(fabs(dX[0]/Vx) + fabs(dX[1]/Ex) + fabs(dX[2]/Vy) + fabs(dX[3]/Ey) )/dim;
      //printf("    -n= %d, x= %g, critere= %g, Vx= %g, Ex= %g, Vy= %g, Ey= %g\n",n,x,critere, Vx,Ex,Vy,Ey );
      n++;
    }
    if(n==Nmax){printf("Newton non convergent : n=%d, critere= %g (courbe_hugoniot_mixte)\n",n,critere ); return 1;}     
    
    Vmel=(1.-x)*VX + x*VY;
    Emel=(1.-x)*EX + x*EY;
    

    // phase X
    Px=fP(rho0_x, E0_x, Vx, Ex, phaseX); Tx=fT(rho0_x, E0_x, Vx, Ex, phaseX);
    Sx=fS(rho0_x, E0_x, Sr_x, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
    // phase Y
    Py=fP(rho0_y, E0_y, Vy, Ey, phaseY); Ty=fT(rho0_y, E0_y, Vy, Ey, phaseY);
    Sy=fS(rho0_y, E0_y, Sr_y, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;

    P=(Px+Py)/2;
    T=(Tx+Ty)/2;
    G=(Gx+Gy)/2;
    S=(1.-x)*Sx + x*Sy;

    
    tabV[i]=Vmel;
    tabE[i]=Emel;
    tabP[i]=P;
    tabT[i]=T;
    tabG[i]=G;
    tabS[i]=S;


    //printf(" --i= %d, x= %g : n_iter= %d, critere= %g, Vmel= %g, Emel= %g, P= %g , T= %g\n",i,x,n,critere,Vmel,Emel,tabP[i],tabT[i] );
  }

  free(A); free(InvA);
  free(Delta); free(dX);
  
  return 0;
}

// ************************************************************************************************************
//    Calcul de la celerite du son et de la derivee fondamentale 

/* Valeur d'entree :
    - Va,Ea, Vb,Eb
    - Vmel, Emel

    - A,InvA,dX : pour limiter les allocations 

  Valeur sortie
    - Vx,Ex, Vy,Ey, 
    - x : fraction massique de la phase Y

   double* A et InvA permet une seule allocation 
*/
int Newton_mixte(double* A, double* InvA, double* Delta, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Vmel, double Emel, double Vx, double Ex, double Vy, double Ey, double x, double* pVX, double* pEX, double* pVY, double* pEY, double* pX, double* dX){
  int dim=5;

  double rho0x=pcoeffX[0], E0x=pcoeffX[1], Srx=pcoeffX[2];
  double rho0y=pcoeffY[0], E0y=pcoeffY[1], Sry=pcoeffY[2];

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(rho0x, E0x, Vx, Ex, phaseX); Tx=fT(rho0x, E0x, Vx, Ex, phaseX);
  Sx=fS(rho0x, E0x, Srx, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(rho0x,Vx,phaseX); dPex=dPdE(rho0x, phaseX);
  dTvx=dTdV(rho0x,Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(rho0y, E0y, Vy, Ey, phaseY); Ty=fT(rho0y, E0y, Vy, Ey, phaseY);
  Sy=fS(rho0y, E0y, Sry, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(rho0y, Vy, phaseY); dPey=dPdE(rho0y, phaseY);
  dTvy=dTdV(rho0y, Vy, phaseY); dTey=dTdE(phaseY);
  
  /*
  printf("    phaseX=%d  phaseY=%d\n",phaseX,phaseY );
  printf("    rho0x=%g E0x=%g Srx=%g\n", rho0x, E0x, Srx);
  printf("    rho0y=%g E0y=%g Sry=%g\n", rho0y, E0y, Sry);
  printf("    Vx=%g Ex=%g\n",Vx,Ex );
  printf("    Vy=%g Ey=%g\n",Vy,Ey );
  printf("    Px= %.3g, Tx= %.3g, Sx= %.3g, Gx= %.3g\n",Px*1e-9,Tx,Sx,Gx );
  printf("    Py= %.3g, Ty= %.3g, Sy= %.3g, Gy= %.3g\n",Py*1e-9,Ty,Sy,Gy );
  printf("    Px-Py=%.3g, Tx-Ty=%.3g, Sx-Sy=%.3g, Gx-Gy= %.3g\n",(Px-Py)*1e-9,Tx-Ty,Sx-Sy,Gx-Gy );
  */
  
  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx; 
  dGex=Vx*dPex - Sx*dTex;
  
  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;   Delta[2]=Gy-Gx;
  Delta[3]=Vmel - (1.-x)*Vx - x*Vy;
  Delta[4]=Emel - (1.-x)*Ex - x*Ey;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx; A[dim*0+1]=dPex; A[dim*0+2]=-dPvy; A[dim*0+3]=-dPey; A[dim*0+4]=0.;  
  
  A[dim*1+0]=dTvx; A[dim*1+1]=dTex; A[dim*1+2]=-dTvy; A[dim*1+3]=-dTey; A[dim*1+4]=0.; 
  
  A[dim*2+0]=dGvx; A[dim*2+1]=dGex; A[dim*2+2]=-dGvy; A[dim*2+3]=-dGey; A[dim*2+4]=0.; 
  
  A[dim*3+0]=1.-x; A[dim*3+1]=0.;   A[dim*3+2]=x;     A[dim*3+3]=0;     A[dim*3+4]=Vy-Vx; 
  
  A[dim*4+0]=0.;   A[dim*4+1]=1.-x; A[dim*4+2]=0;     A[dim*4+3]=x;     A[dim*4+4]=Ey-Ex; 
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %g  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g (\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %g  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Mutiplcation matricielle de A et InvA
  /*
  double res;
  //printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      //printf("A*InvA[%d][%d]= %g ",i,j,res );
    }
    //printf("\n");
  }
  //printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %.15lf\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1./2;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];
  *pX=x+correc*dX[4];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */
}


/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère
  - Vmel, Emel : Point dont on sait qu'il est dans la zone de mélange (alpha/beta) (X/Y)
  - phaseX, phaseY : les deux phases dans la zone de mélange

Arguments de sortie :
  - Vx, Ex : valeur de la phase X
  - Vy, Ey : valeur de la phase Y
  -   x    : fraction massique de Y
*/
int VE_mixte(int Nmax, double epsilon, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, 
             double* pVx, double* pEx, double* pVy, double* pEy, double* px){
  int err, n=1, dim=5;
  double errV, errE;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double Vx, Ex, Vy, Ey, x;
  double VX, EX, VY, EY, X;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  x=1./2;
  
  /*
  double rho0x=pcoeffX[0], E0x=pcoeffX[1], Srx=pcoeffX[2];
  double rho0y=pcoeffY[0], E0y=pcoeffY[1], Sry=pcoeffY[2];

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx;
  double Py, Ty, Sy, Gy;
  // phase X
  Px=fP(rho0x, E0x, Vx, Ex, phaseX);      Tx=fT(rho0x, E0x, Vx, Ex, phaseX);
  Sx=fS(rho0x, E0x, Srx, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(rho0y, E0y, Vy, Ey, phaseY);      Ty=fT(rho0y, E0y, Vy, Ey, phaseY);
  Sy=fS(rho0y, E0y, Sry, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  */

  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
  
  n=0;
  while(critere>=epsilon && n<Nmax){

    err=Newton_mixte(A, InvA, Delta, phaseX, phaseY, pcoeffX, pcoeffY, Vmel, Emel, Vx, Ex, Vy, Ey, x, &VX, &EX, &VY, &EY, &X, dX);

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    x=X;
    errV=((1-x)*Vx + x*Vy-Vmel)/Vmel;
    errE=((1-x)*Ex + x*Ey-Emel)/Emel;

    //critere=fabs(dX[0]/Vx) +fabs(dX[1]/Ex) +fabs(dX[2]/Vy) +fabs(dX[3]/Ey)+fabs(dX[4]/x);
    critere=fabs(dX[4]);

    //printf(" --n= %d, critere= %g (cG_equilibre)\n",n,critere );


    n++;
  }
  //printf(" --n= %d, critere= %g (cG_equilibre)\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if( isnan(critere) ){printf("Erreur critere=%g , V=%g; E=%g; zone: %d%d (VE_mixte)\n ",critere,Vmel,Emel,phaseX, phaseY ); return 1;}
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d (Nmax=%d) , critere= %g V=%g, E=%g (VE_mixte)\n",n,Nmax,critere, Vmel, Emel );return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *px=x;
   
  /*
  if(x<0 || 1<x){
    printf("Erreur fraction massique= %g (VE_mixte)\n",x );
    printf("Newton_mixte : -n= %d , critere= %lf \n",n,critere );
    return 1;
  }
  */

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(dX);
  free(Delta);
  
  return 0;
}

// Discretisation de courbes sur phase pur en cherchant E
int courbe_hugoniot_pure_E(int Nmax, double epsilon, int N, int phase, double* pcoeff, double V0, double E0, double Vpt, double Ept, double Ppt, double Vfin, double *tabV, double* tabE, double* tabP, double* tabT){
  // portion de courbe 0-1
  int n=0;
  double f,f_prime,dE;
  
  double rho_0=pcoeff[0], E_0=pcoeff[1], S_r=pcoeff[2];

  double V,E;
  V=V0; E=E0;
  double hv=(Vfin-V0)/(N-1);
  double critere=1.;
  tabV[0]=V0;  tabE[0]=E0;
  tabP[0]=fP(rho_0, E_0 ,V0,E0,phase);  tabT[0]=fT(rho_0, E_0, V0,E0,phase);
  

  for(int i=1; i<N; i++){
    n=0;  critere=1.;
    V=tabV[i-1];  E=tabE[i-1];
    V+=hv;
    tabV[i]=V;
    while(critere>epsilon && n<Nmax){
      f=E-Ept - fP(rho_0, E_0, V,E,phase)*(Vpt-V)/2;
      f_prime=1-dPdE(rho_0, phase)*(Vpt-V)/2;

      dE=-f/f_prime;
      critere=fabs(dE/E);
      E+=dE;
      //printf(" --n= %d , dE= %g\n",n,dE );
      n++;
    }
    if(n==Nmax){printf("Newton 1D non convergent : n= %d, critere= %g (courbe_hugoniot_pure)\n",n,critere ); return 1;}
    tabE[i]=E;
    tabP[i]=fP(rho_0, E_0,V,E,phase);
    tabT[i]=fT(rho_0, E_0, V,E,phase);
    //printf(" **i= %d, critere= %g, V= %g, E= %g, P= %g\n",i,critere,V,E,tabP[i] );
  }

  
  return 0;
}


/* Algo Newton-Raphson à P fixé pour calculer un point de la courbe d'Hugoniot dans la zone de mélange X/Y
*/
int Newton_hugoniot_mixte(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Vpt, double Ept, double Ppt, double Vx, double Ex, double Vy, double Ey, double* VX, double* EX, double* VY, double* EY, double* px){
  int dim=5;

  double rho0x=pcoeffX[0], E0x=pcoeffX[1], Srx=pcoeffX[2];
  double rho0y=pcoeffY[0], E0y=pcoeffY[1], Sry=pcoeffY[2];

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(rho0x, E0x, Vx, Ex, phaseX); Tx=fT(rho0x, E0x, Vx, Ex, phaseX);
  Sx=fS(rho0x, E0x, Srx, Vx, Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(rho0x,Vx,phaseX); dPex=dPdE(rho0x, phaseX);
  dTvx=dTdV(rho0x,Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(rho0y, E0y, Vy, Ey, phaseY); Ty=fT(rho0y, E0y, Vy, Ey, phaseY);
  Sy=fS(rho0y, E0y, Sry, Vy, Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(rho0y, Vy, phaseY); dPey=dPdE(rho0y, phaseY);
  dTvy=dTdV(rho0y, Vy, phaseY); dTey=dTdE(phaseY);


  double dGvx, dGex, dGvy, dGey;
  dGvx=Vx*dPvx - Sx*dTvx;
  dGex=Vx*dPex - Sx*dTex;

  dGvy=Vy*dPvy - Sy*dTvy; 
  dGey=Vy*dPey - Sy*dTey;
  
  double P2=(Pcst+Ppt)/2;

  double x=*px;

  double Vmel=(1-x)*Vx + x*Vy;
  double Emel=(1-x)*Ex + x*Ey;

  Delta[0]=Px-Pcst; Delta[1]=Py-Pcst;
  Delta[2]=Ty-Tx;   Delta[3]=Gy-Gx;
  Delta[4]=Ept - Emel + P2*(Vpt - Vmel);
   
  
  //for(int j=0; j<dim; j++){
  //  printf("Delta[%d]= %.15lf  ",j, Delta[j] );
  //}
  //printf("\n");
  


  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=-dPvx;     A[dim*0+1]=-dPex; A[dim*0+2]=0.;    A[dim*0+3]=0.;    A[dim*0+4]=0.;
 
  A[dim*1+0]=0.;        A[dim*1+1]=0.;    A[dim*1+2]=-dPvy; A[dim*1+3]=-dPey; A[dim*1+4]=0.;

  A[dim*2+0]=dTvx;      A[dim*2+1]=dTex;  A[dim*2+2]=-dTvy; A[dim*2+3]=-dTey; A[dim*2+4]=0.;

  A[dim*3+0]=dGvx;      A[dim*3+1]=dGex;  A[dim*3+2]=-dGvy; A[dim*3+3]=-dGey; A[dim*3+4]=0.;

  A[dim*4+0]=P2*(1-x);  A[dim*4+1]=(1-x); A[dim*4+2]=P2*x;  A[dim*4+3]=x;     A[dim*4+4]=Ey - Ex + P2*(Vy-Vx);
  
  /*
  printf("  * A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("A[%d][%d]= %.15lf  ",i,j, A[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");  
  */


  double det_a=inverse_matrice_pivot(A, dim, InvA);
  //printf("detA= %.15lf\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a); return 1;}
  
  /*
  printf("  * InvA = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      printf("InvA[%d][%d]= %.15lf  ",i,j, InvA[dim*i+j] );
    }
    printf("\n");
  }
  printf("\n");
  */
  
  
  // Mutiplcation matricielle de A et InvA
  /*
  double res;
  printf("  * InvA*A = *\n");
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      res=0;
      for(int k=0; k<dim; k++){
        res+=InvA[dim*i+k]*A[dim*k+j];
      }
      printf("A*InvA[%d][%d]= %.15lf ",i,j,res );
    }
    printf("\n");
  }
  printf("\n");
  */


  // Mutiplcation matricielle
  for(int i=0; i<dim; i++){
    dX[i]=0;
    for(int j=0; j<dim; j++){
      dX[i]+=InvA[dim*i+j]*Delta[j];
    }
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
   
  double correc=1./2;

  *VX=Vx+correc*dX[0];
  *EX=Ex+correc*dX[1];
  *VY=Vy+correc*dX[2];
  *EY=Ey+correc*dX[3];
  *px=x +correc*dX[4];

  return 0;
}


/*  Calcul des points de l'Hugoniot dans la zone de mélange X/Y 
    à partir d'une discrétisation en pression [Pdeb, Pfin] à N points 
*/
int courbe_hugoniot_mixte(int Nmax, double epsilon, int N, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double Pdeb, double Pfin, double Vpt, double Ept, double Ppt, double V0x, double E0x, double V0y, double E0y, double *tabV, double* tabE, double* tabP, double* tabT, double* tabX){
  // portion de courbe 0-1
  int n=0;
  
  int dim=5;
  int err;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));

  double Vx,Ex,Vy,Ey;
  double VX,EX,VY,EY;
  double x;

  double rho0_x=pcoeffX[0], E0_x=pcoeffX[1], Sr_x=pcoeffX[2];
  double rho0_y=pcoeffY[0], E0_y=pcoeffY[1], Sr_y=pcoeffY[2];

  double Vmel, Emel;
  double hp=(Pfin-Pdeb)/(N-1);
  double critere=1.;
  tabV[0]=V0x;  tabE[0]=E0x;
  tabP[0]=Pdeb; tabT[0]=fT(rho0_x, E0_x, V0x,E0x,phaseX);
  tabX[0]=0;
  double Pcst=Pdeb;
  
  Vx=V0x;   Vy=V0y;
  Ex=E0x;   Ey=E0y;

  for(int i=1; i<N; i++){
    n=0;  critere=1.;
    Pcst+=hp;
    tabP[i]=Pcst;
    x=1./2;
    while(critere>epsilon && n<Nmax){
      // Ici dim=4   // MODIFER Newton_double pour limiter les allocations
      //Newton_double(A, InvA, Delta, dX, Pcst, phaseX, phaseY, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY);
      err=Newton_hugoniot_mixte(A, InvA, Delta, dX, Pcst, phaseX, phaseY, pcoeffX, pcoeffY, Vpt, Ept, Ppt, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, &x);
      if(err){return 1;}
      //printf(" --n= %d , critere= %g\n",n,critere );
      Vx=VX;   Vy=VY;
      Ex=EX;   Ey=EY;

      critere=(fabs(dX[0]/Vx) + fabs(dX[1]/Ex) + fabs(dX[2]/Vy) + fabs(dX[3]/Ey) + fabs(dX[4]/x))/dim;
      //printf("*i= %d -n= %d, critere= %g, Vx= %g, Ex= %g, Vy= %g, Ey= %g, x= %g\n",i,n,critere, Vx,Ex,Vy,Ey,x );
      n++;
    }
    if(n==Nmax){printf("Newton non convergent : n=%d, critere= %g (courbe_hugoniot_pure)\n",n,critere ); return 1;}     
    
    Vmel=(1.-x)*VX + x*VY;
    Emel=(1.-x)*EX + x*EY;
    
    tabV[i]=Vmel;
    tabE[i]=Emel;
    tabT[i]=fT(rho0_x, E0_x, Vx,Ex,phaseX);
    tabX[i]=x;
    //printf(" **i= %d, critere= %g, Vmel= %g, Emel= %g, x= %g,  P= %g\n",i,critere,Vmel,Emel,x,tabP[i] );
  }
  

  free(A); free(InvA);
  free(Delta); free(dX);
  
  return 0;
}

// Fonction pour tester le calcul de la vitesse du son et de la derivée fondamentale en zone pure
int fct_tst_c_G(double* pcoeffgamma, double* pcoeffliq, double* Vpc, double* Epc){
  
  int err; 
  int Nmax=100;
  double epsilon=1e-11;
  
  int ndV=5;
  double pfacdV[5]={1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

  double dV, dE, x, c, G; 
  double Vx, Ex, Vy, Ey;
  double V, E, facdV, rho_0, E_0, Sr;
  int phase; 

  double prho0[6]={rho01, pcoeffgamma[0], pcoeffliq[0]};
  double pE0[6]=  {E01, pcoeffgamma[1], pcoeffliq[1]};
  double pSr[6]=  {Sr1, pcoeffgamma[2], pcoeffliq[2]};

  double V0=Vpc[0], E0=Epc[0];
  double V1=Vpc[1], E1=Epc[1];
  double V2=Vpc[2], E2=Epc[2];
  double V3=Vpc[3], E3=Epc[3];
  double V4=Vpc[4], E4=Epc[4];
  double V5=Vpc[5], E5=Epc[5];
  
  int PHASEpc[6]={1,1,2,2,2,3};

  double P, dPv, dPe, dPvv, Pp, Pm;
  double dPsdV, d2PsdV2;

  for(int i=0; i<6; i++){
    printf("\n\n############  Point %d  ###############\n",i );
    phase=PHASEpc[i];
    V=Vpc[i]; E=Epc[i];
    dV=V*facdV;
    rho_0=prho0[phase-1]; E_0=pE0[phase-1]; Sr=pSr[phase-1];

    P=fP(rho_0, E_0, V, E, phase);

    // analytique  --------------------------------
    dPv=dPdV(rho_0, V, phase);
    dPe=dPdE(rho_0, phase);
    dPvv=d2PdV2(rho_0, V, phase);
    
    c= -V*V*(dPv-P*dPe);
    
    G= +0.5*V*V*V;
    G*=(dPvv - dPv*dPe + P*dPe*dPe);
    G*=1.0/c;
    
    c=sqrt(c);

    printf("analytique   : \n");
    printf("c=%+.5e   G=%+.5e   dPdV=%+.5e   d2PdV2=%+.5e \n", c, G, dPv-P*dPe, dPvv - dPv*dPe + P*dPe*dPe );
    
    // diff. finies  ------------------------------ 
    printf("diff. finies : \n" );
    for(int j=0; j<ndV; j++){
      facdV=pfacdV[j];

      dV=V*facdV;
      dE=-P*dV;

      Pp=fP(rho_0, E_0, V+dV, E+dE, phase);
      Pm=fP(rho_0, E_0, V-dV, E-dE, phase);
      
      dPsdV=0.5*(Pp-Pm)/dV;
      c=sqrt(-V*V*dPsdV);

      d2PsdV2=(Pm+Pp-2*P)/(dV*dV);
      G=-0.5*V*d2PsdV2/dPsdV;
      
      printf("c=%+.5e   G=%+.5e   dPdV=%+.5e   d2PdV2=%+.5e    facdV=%g\n",c,G,dPsdV,d2PsdV2,facdV );
    }
  }

  return 0;
}

/*  Calcul de la célérite du son C et de la dérivée fondamentale G à l'équilibre 
     dans le cas monophase avec une formule analytique 
       et le cas biphasique phase X et Y
*/
int cG_equilibre(int cas, int phaseX, int phaseY, double* pcoeffX, double* pcoeffY, double** VE_init, double V, double E, double P, double* pc, double* pG){
  int err; 
  int Nmax=100;
  double epsilon=1e-11;

  double dV=V*1e-4;
  double dE=-P*dV;
  double x; 

  double Vx, Ex, Vy, Ey;
  double V0x, E0x, V0y, E0y;
  
  double VbetaT, EbetaT, VgammaT, EgammaT, VliqT, EliqT; 
  V0x=VE_init[0][0];  E0x=VE_init[0][1]; 
  V0y=VE_init[1][0];  E0y=VE_init[1][1]; 

  double rho0_x=pcoeffX[0], E0_x=pcoeffX[1], Sr_x=pcoeffX[2];
  double rho0_y=pcoeffY[0], E0_y=pcoeffY[1], Sr_y=pcoeffY[2];

  // Grandeurs thermodynamiques de la phase X
  double Px, dPvx, dPex, dPvvx;
  Px=fP(rho0_x, E0_x, V, E, phaseX);
  dPvx=dPdV(rho0_x,V,phaseX);
  dPex=dPdE(rho0_x, phaseX);
  dPvvx=d2PdV2(rho0_x, V, phaseX);

  // biphase 
  double Pp, Pm; 
  double dPdV,d2PdV2;

  
  // cas monophase  (phase X)
  if (cas==0){

    *pc= -V*V*(dPvx - P*dPex);

    *pG= +0.5*V*V*V;
    *pG*=(dPvvx - dPvx*dPex + P*dPex*dPex);
    *pG*=1.0/(*pc);

    *pc=sqrt(*pc);

    return 0;
  }
  // phase unique avec taux de variation
  if (cas==2){

    Pp=fP(rho0_x, E0_x, V+dV, E+dE, phaseX);
    Pm=fP(rho0_x, E0_x, V-dV, E-dE, phaseX);
    
    dPdV=0.5*(Pp - Pm)/dV;
    printf("    dPdV=%g  et dPv=%g\n",dPdV,dPvx );
    *pc=sqrt(-V*V*dPdV);  

    d2PdV2=(Pm+Pp-2*P)/(dV*dV);
    printf("    d2PdV2=%g  et dPvv=%g\n",d2PdV2,dPvvx );
    *pG=-0.5*V*d2PdV2/dPdV;

    return 0;
  }
  // biphase X et Y
  else if (cas==1){
    err=VE_mixte(Nmax, epsilon, phaseX, phaseY, pcoeffX, pcoeffY, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x); if (err){ printf("(cG_equilibre)\n"); return err;}
    Pp=0.5*(fP(rho0_x, E0_x, Vx, Ex, phaseX)+fP(rho0_y, E0_y, Vy, Ey, phaseY)); 

    err=VE_mixte(Nmax, epsilon, phaseX, phaseY, pcoeffX, pcoeffY, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x); if (err){ printf("(cG_equilibre)\n"); return err;}
    Pm=0.5*(fP(rho0_x, E0_x, Vx, Ex, phaseX)+fP(rho0_y, E0_y, Vy, Ey, phaseY));
    
    
    dPdV=0.5*(Pp - Pm)/dV;
    *pc=sqrt(-V*V*dPdV);  

    d2PdV2=(Pm+Pp-2*P)/(dV*dV);
    *pG=-0.5*V*d2PdV2/dPdV;
     
    return 0;
  }

  printf("erreur cas=%d (cG_equilibre/recalage_coeff_dPdT.c)\n",cas );
  return 1;
}


int calcul_coeff(double dPdT12, double dPdT13, double* pcoeffgamma, double* pcoeffliq ){
  int err; 

  //      Calcul des coefficients (rho0,E0,Sr) en fonction de dPdT gamma et liquide
  // Calcul du point (P0,T0) de la phase beta (Vbeta,Ebeta)
  double P0gamma=P02, T0gamma=T02;
  double V_beta_gamma,E_beta_gamma,S_beta_gamma;
  err=fVE(1,rho01,E01,P0gamma,T0gamma,&V_beta_gamma,&E_beta_gamma); if(err){printf("err print_double\n"); return err;}
  S_beta_gamma=fS(rho01, E01, Sr1, V_beta_gamma, E_beta_gamma, 1);

  double P0liq=P03, T0liq=T03;
  double V_beta_liq,E_beta_liq,S_beta_liq;
  err=fVE(1,rho01,E01,P0liq,T0liq,&V_beta_liq,&E_beta_liq); if(err){printf("err print_double\n"); return err;}
  S_beta_liq=fS(rho01, E01, Sr1, V_beta_liq, E_beta_liq, 1);

  // point de la phase gamma et liquide
  double V_gamma, E_gamma, S_gamma;
  double V_liq, E_liq, S_liq;
  
      
  // Calcul du point (P0,T0) de la phase gamma (Vgamma,Egamma) et de la phase liq (Vliq,Eliq)    ! Clapeyron !
  V_gamma=V_beta_gamma + Dv12;
  S_gamma=S_beta_gamma + dPdT12*Dv12;
  E_gamma=E_beta_gamma + T0gamma*dPdT12*Dv12 - P0gamma*Dv12;

  V_liq=V_beta_liq + Dv13;
  S_liq=S_beta_liq + dPdT13*Dv13;
  E_liq=E_beta_liq + T0liq*dPdT13*Dv13 - P0liq*Dv13;
      
  // Calcul des coeffs (rho0_b,E0_b,Sr_b)  (Newton)
  err= newton_coeff(P0gamma, T0gamma, V_gamma, E_gamma, S_gamma, pcoeffgamma, 2); if(err){ printf("err (print_double)\n"); return err;}
  err= newton_coeff(P0liq, T0liq, V_liq, E_liq, S_liq, pcoeffliq, 3);             if(err){ printf("err (print_double)\n"); return err;}
  
  return 0;
}


// ***********************************************************************************************************************
//  Calcul d'une courbe d'Hugoniot
/*    à partir des coeffcients (rho0,E0,Sr)
      puis on écrit la discretisation de la courbe dans un ficheir texte file
*/
int hugoniot(double* pcoeffbeta, double* pcoeffgamma, double* pcoeffliq, FILE* file){

  int err;
  int N, Nx, Nmax;
  double epsilon, epsilon_EOS=1e-13;
  int phase, phaseX, phaseY;
  double Vfin;
  double Pdeb, Pfin;
  double Vdeb,Edeb; 
  double Vinit, Einit;
  double D,U,u;
  double dV,x;
  
  double C1,C2,C3,C4,C5;
  double G1,G2,G3,G4,G5;
  double S1,S2,S3,S4,S5;
  double* lambda=malloc(3*sizeof(double));

  double** VE_init=alloctabd(2,2);


  double Stemp, Ctemp;
  double G1_diff,G2_diff,G3_diff,G4_diff,G5_diff;

  
  double rho0beta=pcoeffbeta[0],   E0beta=pcoeffbeta[1],   Srbeta=pcoeffbeta[2];
  double rho0gamma=pcoeffgamma[0], E0gamma=pcoeffgamma[1], Srgamma=pcoeffgamma[2];
  double rho0liq=pcoeffliq[0],     E0liq=pcoeffliq[1],     Srliq=pcoeffliq[2];

  // Calcul du point triple  (Newton)
  double PT,TT;
  double VaT,EaT,VbT,EbT,VcT,EcT;
  point_triple(pcoeffgamma, pcoeffliq, &PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);
  double** VE_triple=alloctabd(3,2);
  VE_triple[0][0]=VaT; VE_triple[0][1]=EaT; 
  VE_triple[1][0]=VbT; VE_triple[1][1]=EbT; 
  VE_triple[2][0]=VcT; VE_triple[2][1]=EcT; 


  // ************************************************************************************************
  //      Calcul des points caracteristiques de l'Hugoniot
  
  Nmax=100;  epsilon=4e-14;

  // Point 0
  printf("**************   Point 0   ****************************************** \n");
  double P0,T0, C0,S0,G0,G0_diff; 
  //P0=1.01325e5; T0=300;
  P0=0.; T0=300;

  double V0,E0;
  err=fVE(1, rho0beta, E0beta, P0, T0, &V0, &E0);
  if(err){return 1;}
  
  VE_init[0][0]=VE_triple[0][0];  VE_init[0][1]=VE_triple[0][1];
  VE_init[1][0]=VE_triple[0][0];  VE_init[1][1]=VE_triple[1][1];

  err=cG_equilibre(0, 1, 1, pcoeffbeta, pcoeffgamma, VE_init, V0, E0, P0, &C0, &G0_diff);
  if (err){printf("Erreur =%d (\n",err); return err;}
  
  S0=fS(rho0beta, E0beta, Srbeta, V0,E0,1);
  G0=E0 + P0*V0 - T0*S0;
  
  // comparaison analytique / taux de variation

  err=cG_equilibre(0, 1, 1, pcoeffbeta, pcoeffgamma, VE_init, V0, E0, P0, &C0, &G0_diff);
  printf("analytique c=%g   der_fond=%g \n", C0, G0_diff);
  err=cG_equilibre(2, 1, 1, pcoeffbeta, pcoeffgamma, VE_init, V0, E0, P0, &C0, &G0_diff);
  printf("taux de variation c=%g   der_fond=%g \n", C0, G0_diff);

  printf("  V0= %g,E0= %g, P0= %g, T0= %g, S0=%g, C0=%g, G0=%g\n",V0,E0,P0,T0,S0,C0,G0 );

  double V1,E1, V2,E2, V3,E3, V4,E4, V5,E5;
  double P1,T1, P2,T2, P3,T3, P4,T4, P5,T5;

  double V0x, E0x, V0y, E0y;
  double Vx, Ex, Vy, Ey;

  double E,V;


  // CALCUL DU POINT 1
  printf("\n********************************************************\nPoint 1\n");
  V0x=123.8e-6;  E0x=51.4e3;
  V0y=121.9e-6;  E0y=79e3;
  
  err= hugoniot_pts12(1, Nmax, epsilon, 1, 2, pcoeffbeta, pcoeffgamma, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey, &P1, &T1, &S1, &G1);   if(err){return 1;}
  V1=Vx; E1=Ex;

  VE_init[0][0]=VE_triple[0][0];  VE_init[0][1]=VE_triple[0][1];
  VE_init[1][0]=VE_triple[0][0];  VE_init[1][1]=VE_triple[1][1];
  // comparaison des 2 méthodes
  err=cG_equilibre(0, 1, 1, pcoeffbeta, pcoeffgamma, VE_init, V1, E1, P1, &C1, &G1_diff);
  printf("Point 1 : 1re methode analytique : c=%g   der_fond=%g\n",C1,G1_diff );
  err=cG_equilibre(2, 1, 1, pcoeffbeta, pcoeffgamma, VE_init, V1, E1, P1, &C1, &G1_diff);
  printf("Point 1 : 2re methode taux de variation (monophase) : c=%g der_fond=%g\n",C1,G1_diff );
  err=cG_equilibre(1, 1, 2, pcoeffbeta, pcoeffgamma, VE_init, V1, E1, P1, &C1, &G1_diff);
  printf("Point 1 : 3nd methode taux de variation (biphase): c=%g  der_fond=%g\n",C1,G1_diff );
 
  printf("  V1= %g,E1= %g, P1= %g, T1= %g, S1=%g, C1=%g, G1=%g\n",V1,E1,P1,T1,S1,C1,G1 );


  // CALCUL DU POINT 2
  printf("\n********************************************************\nPoint 2\n");
  V0x=122.5e-6;  E0x=49.6e3;
  V0y=120.6e-6;  E0y=77e3;

  Nmax=1e4;  epsilon=2e-14;
  err= hugoniot_pts12(2, Nmax, epsilon, 1, 2, pcoeffbeta, pcoeffgamma, V1, E1, P1, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey, &P2, &T2, &S2, &G2);  if(err){return 1;}
  V2=Vy; E2=Ey;
  
  VE_init[0][0]=VE_triple[0][0];  VE_init[0][1]=VE_triple[0][1];
  VE_init[1][0]=VE_triple[0][0];  VE_init[1][1]=VE_triple[1][1];
  // comparaison des 2 méthodes
  err=cG_equilibre(0, 2, 2, pcoeffgamma, pcoeffgamma, VE_init, V2, E2, P2, &C2, &G2_diff);
  printf("Point 2 : 1re methode analytique : c=%g G=%g\n",C2, G2_diff );
  err=cG_equilibre(1, 1, 2, pcoeffbeta, pcoeffgamma, VE_init, V2, E2, P2, &C2, &G2_diff);
  printf("Point 2 : 2nd methode taux de variation : c=%g G=%g\n",C2, G2_diff );
  
  printf("  V2= %g, E2= %g , P2= %g, T2= %g, S2=%g, C2=%g, G2=%g\n",V2,E2,P2,T2,S2,C2,G2 );


  // CALCUL DU POINT 3
  printf("\n********************************************************\nPoint 3\n");
  double* Vpt, *Ept, * Ppt;
  Vpt=malloc(2*sizeof(double));  Ept=malloc(2*sizeof(double));
  Ppt=malloc(2*sizeof(double));  
  
  Vpt[0]=V0; Vpt[1]=V1;  
  Ept[0]=E0; Ept[1]=E1;
  Ppt[0]=P0;  Ppt[1]=P1;
  
  phase=2;
  Vinit=126e-6;  Einit=100e3;
  
  Nmax=1e3; epsilon=5e-15;
  err= hugniot_pt3(Nmax, epsilon, pcoeffgamma, Vpt, Ept, Ppt, Vinit, Einit, &V3, &E3);
  P3=fP(rho0gamma, E0gamma, V3,E3,2);  T3=fT(rho0gamma, E0gamma, V3,E3,2);  S3=fS(rho0gamma, E0gamma, Srgamma, V3,E3,2); G3=E3 + P3*V3 - T3*S3;

  err=cG_equilibre(0, 2, 2, pcoeffgamma, pcoeffgamma, VE_init, V3, E3, P3, &C3, &G3_diff);

  printf("  V3= %g, E3= %g , P3= %g, T3= %g, S3=%g, C3=%g, G3=%g\n",V3,E3,P3,T3,S3,C3,G3 );

  free(Vpt);  free(Ept);  free(Ppt);
  


  // CALCUL DU POINT 4  debut du mélange liquide
  printf("\n********************************************************\nPoint 4\n");
  V0x=96e-6;  E0x=8.85e5;
  V0y=95e-6;  E0y=1e6;
  
  Nmax=1e4;  epsilon=4e-14;
  err= hugoniot_pts12(1, Nmax, epsilon, 2, 3, pcoeffgamma, pcoeffliq, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey, &P4, &T4, &S4, &G4);  if(err){return 1;}  
  V4=Vx; E4=Ex;
  
  // méthode analytique
  err=cG_equilibre(0, 2, 2, pcoeffgamma, pcoeffgamma, VE_init, V4, E4, P4, &C4, &G4_diff);
  
  printf("  V4= %g,E4= %g , P4= %g, T4= %g, S4=%g, C4=%g, G4=%g\n",V4,E4,P4,T4,S4,C4,G4 );

  
  // CALCUL DU POINT 5  debut du mélange liquide
  printf("\n********************************************************\nPoint 5\n");
  V0x=95e-6;  E0x=1.02e6;
  V0y=94.2e-6;  E0y=1.22e6;
  
  Nmax=1e4;  epsilon=4e-14;
  err= hugoniot_pts12(2, Nmax, epsilon, 2, 3, pcoeffgamma, pcoeffliq, V0, E0, P0, V0x, E0x, V0y, E0y, &Vx, &Ex, &Vy, &Ey, &P5, &T5, &S5, &G5); if(err){return 1;}
  V5=Vy; E5=Ey;

  // méthode analytique
  err=cG_equilibre(0, 3, 3, pcoeffliq, pcoeffliq, VE_init, V5, E5, P5, &C5, &G5_diff);

  printf("  V5= %g,E5= %g , P5= %g, T5= %g, S5=%g, C5=%g, G5=%g\n",V5,E5,P5,T5,S5,C5,G5 );
  

  // #########################
  // Test du calcul de la celerité du son c et de la dérivée fondamentale G sur les points caracteristiques
  double Vpc[6]={V0, V1, V2, V3, V4, V5};
  double Epc[6]={E0, E1, E2, E3, E4, E5};
  
  err=fct_tst_c_G(pcoeffgamma, pcoeffliq, Vpc, Epc);

  // #########################

  // Calcul de point dans les zones mixtes à fraction massique fixee
  // Point sur la zone mixte 1/2 avec x fixee
  Nx=11;
  double* tabXm12=malloc(Nx*sizeof(double));  
  double* tabVm12=malloc(Nx*sizeof(double));  double* tabEm12=malloc(Nx*sizeof(double));
  double* tabPm12=malloc(Nx*sizeof(double));  double* tabTm12=malloc(Nx*sizeof(double));
  double* tabSm12=malloc(Nx*sizeof(double));  double* tabGm12=malloc(Nx*sizeof(double));
  double* tabCm12=malloc(Nx*sizeof(double));  double* tabGm12_diff=malloc(Nx*sizeof(double));

  for(int i=0; i<Nx; ++i){
    tabXm12[i]=i*(1./(Nx-1));
    lambda[0]=0;   lambda[1]=1;   lambda[2]=0; 
  }

  V0x=123.8e-6;  E0x=51.4e3;
  V0y=121.9e-6;  E0y=79e3;
  
  Nmax=1e4;  epsilon=4e-14;
  phaseX=1; phaseY=2;
  printf("\n\n mixte 12\n");
  err=Pt_hugoniot_mixte_Xfixe(Nmax, epsilon, Nx, phaseX, phaseY, pcoeffbeta, pcoeffgamma, tabXm12, V1, E1, P1, V0x, E0x, V0y, E0y, tabVm12, tabEm12, tabPm12, tabTm12, tabSm12, tabGm12);
 
  VE_init[0][0]=VE_triple[0][0];  VE_init[0][1]=VE_triple[0][1];
  VE_init[1][0]=VE_triple[0][0];  VE_init[1][1]=VE_triple[1][1];
  for(int i=0; i<Nx; ++i){
    x=i*(1./(Nx-1));
    lambda[0]=1.-x;   lambda[1]=x;   lambda[2]=0; 
    
    err=cG_equilibre(1, 1, 2, pcoeffbeta, pcoeffgamma, VE_triple, tabVm12[i], tabEm12[i], tabPm12[i], &tabCm12[i], &tabGm12_diff[i] );
  }

  // fonction pour C

  // Point sur la zone mixte 2/3 avec x fixee
  Nx=11;
  double* tabXm23=malloc(Nx*sizeof(double));  
  double* tabVm23=malloc(Nx*sizeof(double));  double* tabEm23=malloc(Nx*sizeof(double));
  double* tabPm23=malloc(Nx*sizeof(double));  double* tabTm23=malloc(Nx*sizeof(double));
  double* tabSm23=malloc(Nx*sizeof(double));  double* tabGm23=malloc(Nx*sizeof(double));
  double* tabCm23=malloc(Nx*sizeof(double));  double* tabGm23_diff=malloc(Nx*sizeof(double));

  for(int i=0; i<Nx; ++i){  tabXm23[i]=i*(1./(Nx-1));  }

  V0x=V4;  E0x=E4;
  V0y=94.2e-6;  E0y=1.22e6;
  Nmax=1e2;  epsilon=1e-13; 

  phaseX=2; phaseY=3;
  printf("\n\n mixte 23\n");
  err=Pt_hugoniot_mixte_Xfixe(Nmax, epsilon, Nx, phaseX, phaseY, pcoeffgamma, pcoeffliq, tabXm23, V0, E0, P0, V0x, E0x, V0y, E0y, tabVm23, tabEm23, tabPm23, tabTm23, tabSm23, tabGm23);
  
  VE_init[0][0]=V4;  VE_init[0][1]=E4;
  VE_init[1][0]=V5;  VE_init[1][1]=E5;
  for(int i=0; i<Nx; ++i){
    x=i*(1./(Nx-1));
    lambda[0]=0;   lambda[1]=1.-x;   lambda[2]=x; 
    
    err=cG_equilibre(1, 2, 3, pcoeffgamma, pcoeffliq, VE_triple, tabVm23[i], tabEm23[i], tabPm23[i], &tabCm23[i], &tabGm23_diff[i] );
  }

  double D1=V0*sqrt((P1-P0)/(V0-V1));
  double D2=D1;
  double D3=V0*sqrt((P3-P0)/(V0-V3)); //=D1
  double D4=V0*sqrt((P4-P0)/(V0-V4));
  double D5=V0*sqrt((P5-P0)/(V0-V5));

  double u1=V0*(P1-P0)/D1;
  double u2=u1+sqrt((P2-P1)*(V1-V2));
  double u3=V0*(P3-P0)/D3;
  double u4=V0*(P4-P0)/D4;
  double u5=V0*(P5-P0)/D5;

  printf("\n\n\n");
  /*
  printf("  rho0= %.14g; V0= %.14g; E0= %.14g; P0= %.14g; T0= %.14g; S0=%.14g, G0=%.14g; G0_diff=%.14g; C0=%.14g  \n",1./V0, V0*1e6,E0,P0*1e-9,T0, S0,G0,G0_diff,C0 );
  printf("  rho1= %.14g; V1= %.14g; E1= %.14g; P1= %.14g; T1= %.14g; D1= %.14g; u1= %.14g; S1= %.14g; G1= %.14g; G1_diff= %.14g; C1= %.14g\n",1./V1,V1*1e6,E1,P1*1e-9,T1,D1,u1,S1,G1,G1_diff,C1 );
  printf("  rho2= %.14g; V2= %.14g; E2= %.14g; P2= %.14g; T2= %.14g; D2= %.14g; u2= %.14g; S2= %.14g; G2= %.14g; G2_diff= %.14g; C2= %.14g\n",1./V2,V2*1e6,E2,P2*1e-9,T2,D2,u2,S2,G2,G2_diff,C2 );
  printf("  rho3= %.14g; V3= %.14g; E3= %.14g; P3= %.14g; T3= %.14g; D3= %.14g; u3= %.14g; S3= %.14g; G3= %.14g; G3_diff= %.14g; C3= %.14g\n",1./V3,V3*1e6,E3,P3*1e-9,T3,D3,u3,S3,G3,G3_diff,C3 );
  printf("  rho4= %.14g; V4= %.14g; E4= %.14g; P4= %.14g; T4= %.14g; D4= %.14g; u4= %.14g; S4= %.14g; G4= %.14g; G4_diff= %.14g; C4= %.14g\n",1./V4,V4*1e6,E4,P4*1e-9,T4,D4,u4,S4,G4,G4_diff,C4 );
  printf("  rho5= %.14g; V5= %.14g; E5= %.14g; P5= %.14g; T5= %.14g; D5= %.14g; u5= %.14g; S5= %.14g; G5= %.14g; G5_diff= %.14g; C5= %.14g\n",1./V5,V5*1e6,E5,P5*1e-9,T5,D5,u5,S5,G5,G5_diff,C5 );
  */
  
  // format de sortie pour compte-rendu latex
  printf("                 Point 0          Point 1          Point 2          Point 3          Point 4          Point 5\n");
  printf(" P (GPa)         %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",P0*1e-9,P1*1e-9,P2*1e-9,P3*1e-9,P4*1e-9,P5*1e-9 );
  printf(" T (K)           %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",T0,T1,T2,T3,T4,T5 );
  printf(" rho (kg/m³)     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",1./V0,1./V1,1./V2,1./V3,1./V4,1./V5 );
  printf(" tau (cm³/kg)    %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",V0*1e6,V1*1e6,V2*1e6,V3*1e6,V4*1e6,V5*1e6 );
  printf(" E (J/kg)        %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",E0,E1,E2,E3,E4,E5 );
  printf(" u (m/s)         -                %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",u1,u2,u3,u4,u5 );
  printf(" D (m/s)         -                %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",D1,D2,D3,D4,D5 );
  printf(" c (m/s)         %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",C0,C1,C2,C3,C4,C5 );
  printf(" S (J/kg/K)      %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",S0,S1,S2,S3,S4,S5 );
  printf(" G deri fond     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e     %+.5e  \n",G0_diff,G1_diff,G2_diff,G3_diff,G4_diff,G5_diff );
  

  // format de sortie pour Octave
  printf("\n\n [Point 0 , Point 1 , Point 2 , Point 3 , Point 4 , Point 5]\n");
  printf(" Ppt(i,:)=[%g %g %g %g %g %g];  \n",P0*1e-9,P1*1e-9,P2*1e-9,P3*1e-9,P4*1e-9,P5*1e-9 );
  printf(" Tpt(i,:)=[%g %g %g %g %g %g];  \n",T0,T1,T2,T3,T4,T5 );
  printf(" RHOpt(i,:)=[%g %g %g %g %g %g];  \n",1./V0,1./V1,1./V2,1./V3,1./V4,1./V5 );
  printf(" Vpt(i,:)=[%g %g %g %g %g %g];  \n",V0*1e6,V1*1e6,V2*1e6,V3*1e6,V4*1e6,V5*1e6 );
  printf(" Ept(i,:)=[%g %g %g %g %g %g];  \n",E0,E1,E2,E3,E4,E5 );
  printf(" Upt(i,:)=[-99 %g %g %g %g %g]; \n",u1,u2,u3,u4,u5 );
  printf(" Dpt(i,:)=[-99 %g %g %g %g %g];  \n",D1,D2,D3,D4,D5 );
  printf(" Cpt(i,:)=[%g %g %g %g %g %g];   \n",C0,C1,C2,C3,C4,C5 );
  printf(" Spt(i,:)=[%g %g %g %g %g %g];   \n",S0,S1,S2,S3,S4,S5 );
  printf("G_DIFFpt(i,:)=[%g %g %g %g %g %g];  \n",G0_diff,G1_diff,G2_diff,G3_diff,G4_diff,G5_diff );
  

  /*
  printf("\n\n ZONE MIXTE beta/gamma (1/2)\n");
  for(int i=0; i<Nx; i++){
    D=V1*sqrt((tabPm12[i]-P1)/(V1-tabVm12[i]));
    u=u1 + sqrt((tabPm12[i]-P1)*(V1-tabVm12[i]));
    printf("  x= %.14g; P= %.14g; u= %.14g; rho= %.14g; V= %.14g; E= %.14g; T= %.14g; S=%.14g; G=%.14g; G_diff=%.14g; D=%.14g; C=%.14g\n",tabXm12[i], tabPm12[i], u, 1./tabVm12[i], tabVm12[i], tabEm12[i], tabTm12[i], tabSm12[i], tabGm12[i], tabGm12_diff[i], D, tabCm12[i] );
  }
  */

  // format de sortie pour compte-rendu latex
  printf("\n\n ZONE MIXTE beta/gamma (1/2)\n");
  { 
    {
      printf("frac mass");
      for(int i=0; i<Nx; i++){
        printf("            %.2g",tabXm12[i] );
      }
      printf("\n");
      printf("P (GPa)           ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabPm12[i]*1e-9 );
      }
      printf("\n");
      printf("T (K)             ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabTm12[i] );
      }
      printf("\n");
      printf("rho (kg/m³)       ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",1./tabVm12[i] );
      }
      printf("\n");
      printf("tau (cm³/kg)      ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabVm12[i]*1e6 );
      }
      printf("\n");
      printf("E (J/kg)          ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabEm12[i] );
      }
      printf("\n");
      printf("u (m/s)           ");
      for(int i=0; i<Nx; i++){
        D=V1*sqrt((tabPm12[i]-P1)/(V1-tabVm12[i]));
        u=u1 + sqrt((tabPm12[i]-P1)*(V1-tabVm12[i]));
        printf("  %+.5e",u );
      }
      printf("\n");
      printf("D (m/s)           ");
      for(int i=0; i<Nx; i++){
        D=V1*sqrt((tabPm12[i]-P1)/(V1-tabVm12[i]));
        printf("  %+.5e",D );
        if (i==0){ printf("        "); }
      }
      printf("\n");
      printf("c (m/s)           ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabCm12[i] );
      }
      printf("\n");
      printf("S (m/s)           ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabSm12[i] );
      }
      printf("\n");
      printf("G derivé fond     ");
      for(int i=0; i<Nx; i++){
        printf("  %+.5e",tabGm12_diff[i] );
      }
      printf("\n");
    }
  }

  /*
  printf("\n\n ZONE MIXTE gamma/liquide (2/3)\n");
  for(int i=0; i<Nx; i++){
    printf("  x= %.14g; P= %.14g; u= %.14g; rho= %.14g; V= %.14g; E= %.14g; T= %.14g; S=%.14g; G=%.14g; G_diff=%.14g; D=%.14g; C=%.14g\n",tabXm23[i],tabPm23[i],u, 1./tabVm23[i], tabVm23[i], tabEm23[i], tabTm23[i], tabSm23[i], tabGm23[i], tabGm23_diff[i], D, tabCm23[i] );
  }
  for(int i=0; i<Nx; i++){
    D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
    u=V0*(tabPm23[i]-P0)/D;
    printf("  x= %.14g; P= %.14g; u= %.14g; rho= %.14g; V= %.14g; E= %.14g; T= %.14g; S=%.14g; G=%.14g; G_diff=%.14g; D=%.14g; C=%.14g\n",tabXm23[i],tabPm23[i],u, 1./tabVm23[i], tabVm23[i], tabEm23[i], tabTm23[i], tabSm23[i], tabGm23[i], tabGm23_diff[i], D, tabCm23[i] );
  }
  */

  // format de sortie pour compte-rendu latex
  printf("\n\n ZONE MIXTE gamma/liquide (2/3)\n");
  {
    printf("frac mass");
    for(int i=0; i<Nx; i++){
      printf("            %.2g",tabXm23[i] );
    }
    printf("\n");
    printf("P (GPa)           ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabPm23[i]*1e-9 );
    }
    printf("\n");
    printf("T (K)             ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabTm23[i] );
    }
    printf("\n");
    printf("rho (kg/m³)       ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",1./tabVm23[i] );
    }
    printf("\n");
    printf("tau (cm³/kg)      ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabVm23[i]*1e6 );
    }
    printf("\n");
    printf("E (J/kg)          ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabEm23[i] );
    }
    printf("\n");
    printf("u (m/s)           ");
    for(int i=0; i<Nx; i++){
      D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
      u=V0*(tabPm23[i]-P0)/D;
      printf("  %+.5e",u );
    }
    printf("\n");
    printf("D (m/s)           ");
    for(int i=0; i<Nx; i++){
      D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
      printf("  %+.5e",D );
    }
    printf("\n");
    printf("c (m/s)           ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabCm23[i] );
    }
    printf("\n");
    printf("S (m/s)           ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabSm23[i] );
    }
    printf("\n");
    printf("G derivé fond     ");
    for(int i=0; i<Nx; i++){
      printf("  %+.5e",tabGm23_diff[i] );
    }
    printf("\n");
  }
  

  // format de sortie pour compte-rendu latex
  printf("\n\n 2 ZONES MIXTES beta/gamma et gamma/liquide \n");
  {
    printf("frac mass");
    for(int i=2; i<9; i=i+2){
      printf("            %.2g",tabXm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("            %.2g",tabXm23[i] );
    }
    printf("\n");

    printf("P (GPa)           ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabPm12[i]*1e-9 );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabPm23[i]*1e-9 );
    }
    printf("\n");

    printf("T (K)             ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabTm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabTm23[i] );
    }
    printf("\n");

    printf("rho (kg/m³)       ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",1./tabVm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",1./tabVm23[i] );
    }
    printf("\n");

    printf("tau (cm³/kg)      ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabVm12[i]*1e6 );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabVm23[i]*1e6 );
    }
    printf("\n");
    
    printf("E (J/kg)          ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabEm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabEm23[i] );
    }
    printf("\n");

    printf("u (m/s)           ");
    for(int i=2; i<9; i=i+2){
      D=V0*sqrt((tabPm12[i]-P0)/(V0-tabVm12[i]));
      u=V0*(tabPm12[i]-P0)/D;
      printf("  %+.5e",u );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
      u=V0*(tabPm23[i]-P0)/D;
      printf("  %+.5e",u );
    }
    printf("\n");

    printf("D (m/s)           ");
    for(int i=2; i<9; i=i+2){
      D=V0*sqrt((tabPm12[i]-P0)/(V0-tabVm12[i]));
      printf("  %+.5e",D );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      D=V0*sqrt((tabPm23[i]-P0)/(V0-tabVm23[i]));
      printf("  %+.5e",D );
    }
    printf("\n");

    printf("c (m/s)           ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabCm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabCm23[i] );
    }
    printf("\n");

    printf("S (m/s)           ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabSm12[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabSm23[i] );
    }
    printf("\n");

    printf("G derivé fond     ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabGm12_diff[i] );
    }
    printf(" | ");
    for(int i=2; i<9; i=i+2){
      printf("  %+.5e",tabGm23_diff[i] );
    }
    printf("\n");
  }
  
  // ********************************************************************************************************************************
  //         Ecriture de la courbe d'hugoniot dans un fichier texte
  
  {
    // DISCRETISATION DES PORTIONS DE COURBES ENTRE LES =! POINTS
    N=50;
    double* tabV=malloc(N*sizeof(double));  double* tabE=malloc(N*sizeof(double));
    double* tabP=malloc(N*sizeof(double));  double* tabT=malloc(N*sizeof(double));
    double* tabX=malloc(N*sizeof(double));  
    double* tabU=malloc(N*sizeof(double));  double* tabD=malloc(N*sizeof(double));  
    double* tabC=malloc(N*sizeof(double));  double* tabDer_Fond=malloc(N*sizeof(double));  
    double* tabS=malloc(N*sizeof(double));  
    
    double P,C,Der_Fond,S;


    // ZONE ENTRE 0 ET 1  (phase 1)
    printf("ZONE ENTRE 0 ET 1  (phase beta)\n");
    Nmax=2e4; epsilon=1e-13;
    phase=1;
    Vdeb=V0;  Edeb=E0;  
    Vfin=V1;
    err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, pcoeffbeta, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
    if(err){return 1;}

    for(int i=0; i<N; i++){
      V=tabV[i]; E=tabE[i]; P=tabP[i]; 

      D=V0*sqrt((P-P0)/(V0-V));
      U=V0*(P-P0)/D;
      S=fS(rho0beta,E0beta,Srbeta,V,E,1);
      err=cG_equilibre(0, 1, 1, pcoeffbeta, pcoeffbeta, VE_init, V, E, P, &C, &Der_Fond );

      fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,tabT[i],U,D,S,C,Der_Fond,1.0,0.0,0.0 );
    }
    

    // ZONE BIPHASIQUE ENTRE 1 ET 2
    printf("ZONE BIPHASIQUE ENTRE 1 ET 2\n");
    V0x=V1;  E0x=E1;
    V0y=121.9e-6;  E0y=79e3;
    Nmax=50; epsilon=5e-9; 
    phaseX=1;  phaseY=2;
    Pdeb=P1; Pfin=P2;
    err= courbe_hugoniot_mixte(Nmax, epsilon, N, phaseX, phaseY, pcoeffbeta, pcoeffgamma, Pdeb, Pfin, V1, E1, P1, V0x, E0x, V0y, E0y,  tabV, tabE, tabP, tabT, tabX);   if(err){return 1;}
    
    VE_init[0][0]=VE_triple[0][0];  VE_init[0][1]=VE_triple[0][1];
    VE_init[1][0]=VE_triple[0][0];  VE_init[1][1]=VE_triple[1][1];
    for(int i=0; i<N; i++){
      V=tabV[i]; E=tabE[i]; P=tabP[i]; 
      
      D=D1;
      u=u1 + i*(u2-u2)/(N-1);
     
      S = (1-tabX[i])*fS(rho0beta,E0beta,Srbeta,V,E,1) + tabX[i]*fS(rho0gamma,E0gamma,Srgamma,V,E,2);
      err=cG_equilibre(1, 1, 2, pcoeffbeta, pcoeffgamma, VE_init, V, E, P, &C, &Der_Fond );
       

      fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,tabT[i],U,D,S,C,Der_Fond,1-tabX[i],tabX[i],0.0 );
    }
    

    /*
      // ZONE ENTRE 2 ET 3 (phase 2)
      printf("ZONE ENTRE 2 ET 3 (phase gamma)\n");
      Nmax=2e4; epsilon=1e-13; 
      phase=2;
      Vfin=V3;
      Vdeb=V2; Edeb=E2; Pdeb=P2;
      err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, Vdeb, Edeb, V2, E2, P2, Vfin, tabV, tabE, tabP, tabT);
      //err=courbe_hugoniot_pure_P(Nmax, epsilon, N, phase, Vdeb, Pdeb, V2, E2, P2, Vfin, tabV, tabE, tabP, tabT);
      if(err){return 1;}

      for(int i=0; i<N; i++){
        D=u1+V1*sqrt((tabP[i]-P1)/(V1-tabV[i]));
        u=u1+V1*sqrt((tabP[i]-P1)*(V1-tabV[i]));

        fprintf(fV, "%.15lf \n",tabV[i] );
        fprintf(fE, "%.15lf \n",tabE[i] );
        fprintf(fPh, "%.15lf \n",tabP[i] );
        fprintf(fTh, "%.15lf \n",tabT[i] );

        fprintf(fU, "%.15lf \n",U );
        fprintf(fD, "%.15lf \n",D );

        fprintf(fXA, "%g \n",0. );
        fprintf(fXB, "%g \n",1.0 );
        fprintf(fXC, "%g \n",0.); 
      }
    */

    
    // ZONE ENTRE 3 et 4  (phase 2)
    Nmax=2e3; epsilon=1e-13; 
    phase=2;
    Vfin=V4;
    Vdeb=V3; Edeb=E3;
    err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, pcoeffgamma, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
    if(err){return 1;}

    for(int i=0; i<N; i++){
      V=tabV[i]; E=tabE[i]; P=tabP[i]; 

      D=V0*sqrt((P-P0)/(V0-V));
      U=V0*(P-P0)/D;

      S=fS(rho0gamma,E0gamma,Srgamma,V,E,2);
      err=cG_equilibre(0, 2, 2, pcoeffgamma, pcoeffgamma, VE_init, V, E, P, &C, &Der_Fond );

      fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,tabT[i],U,D,S,C,Der_Fond,0.0,1.0,0.0);
    }


    // ZONE BIPHASIQUE ENTRE 4 ET 5
    V0x=V4;  E0x=E4;
    V0y=94.2e-6;  E0y=1.22e6;
    Nmax=2e3;  epsilon=1e-13;
    phaseX=2;  phaseY=3;
    Pdeb=P4;   Pfin=P5;
    err= courbe_hugoniot_mixte(Nmax, epsilon, N, phaseX, phaseY, pcoeffgamma, pcoeffliq, Pdeb, Pfin, V0, E0, P0, V0x, E0x, V0y, E0y,  tabV, tabE, tabP, tabT, tabX);
    if(err){return 1;}

    VE_init[0][0]=V4;  VE_init[0][1]=E4;
    VE_init[1][0]=V5;  VE_init[1][1]=E5;
    for(int i=0; i<N; i++){
      V=tabV[i]; E=tabE[i]; P=tabP[i]; 

      D=V0*sqrt((P-P0)/(V0-V));
      U=V0*(P-P0)/D;

      S = (1-tabX[i])*fS(rho0gamma,E0gamma,Srgamma,V,E,2) + tabX[i]*fS(rho0liq,E0liq,Srliq,V,E,3);
      err=cG_equilibre(1, 2, 3, pcoeffgamma, pcoeffliq, VE_init, V, E, P, &C, &Der_Fond );

      fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,tabT[i],U,D,S,C,Der_Fond,1-tabX[i],tabX[i],0.0 );
    }


    // ZONE APRES 5  (phase 3)
    Nmax=2e3; epsilon=1e-13; 
    phase=3;
    Vfin=50e-6;
    Vdeb=V5;  Edeb=E5;
    err=courbe_hugoniot_pure_E(Nmax, epsilon, N, phase, pcoeffliq, Vdeb, Edeb, V0, E0, P0, Vfin, tabV, tabE, tabP, tabT);
    if(err){return 1;}

    for(int i=0; i<N; i++){
      V=tabV[i]; E=tabE[i]; P=tabP[i]; 

      D=V0*sqrt((P-P0)/(V0-V));
      U=V0*(P-P0)/D;

      S=fS(rho0liq,E0liq,Srliq,V,E,3);
      err=cG_equilibre(0, 3, 3, pcoeffliq, pcoeffliq, VE_triple, V, E, P, &C, &Der_Fond );

      fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,tabT[i],U,D,S,C,Der_Fond,0.0,1.0,0.0 );
    }
  } 
  
  return 0;
}

// ************************************************************************************************************
/*    Permet de lancer le programme hugoniot à  partir des coeff (rho0,E0,Sr) et le fichier texte 
          qu'on lui fournit 
*/
int lancement_hugoniot(){
  int err; 
  double dPdT12, dPdT13;
  
  FILE *fhug;

  // coeff
  double* pcoeffbeta=malloc(3*sizeof(double));
  double* pcoeffgamma=malloc(3*sizeof(double));
  double* pcoeffliq=malloc(3*sizeof(double));
  pcoeffbeta[0]=rho01;
  pcoeffbeta[1]=E01;
  pcoeffbeta[2]=Sr1;


  // #################################################################
  //    calcul coeff (rho0,E0,Sr) de Heuzé du premier jeu de coeff pour comparer avec le rapport de stage
  
    printf("   calcul coeff (rho0,E0,Sr) de Heuzé du premier jeu de coeff pour comparer avec le rapport de stage \n"); 
    pcoeffgamma[0]=7401.45;
    pcoeffgamma[1]=39033.0;
    pcoeffgamma[2]=24.7813;

    pcoeffliq[0]=7026.17;
    pcoeffliq[1]=90660.5;
    pcoeffliq[2]=165.11;
    
    if((fhug = fopen("fichiers/hug0.txt", "w+")) == NULL){printf("erreur ouverture fichier hugoniot (recalage_coeff_dPdT.c)\n"); return 1;}


    err=hugoniot(pcoeffbeta, pcoeffgamma, pcoeffliq, fhug);
    fclose(fhug);
  
  
  // #################################################################
  //    calcul valeur convergés des dPdT de Heuzé
  /*
    printf("\n\n####################################################################################\n");
    printf("            Calcul valeur convergés des dP/dT de Heuzé \n");
    printf("####################################################################################\n");

    dPdT12=-2.13e7; dPdT13=3.125e7;

    if((fhug = fopen("fichiers/hug1.txt", "w+")) == NULL){printf("erreur ouverture fichier hugoniot (recalage_coeff_dPdT.c)\n"); return 1;}
    
    // Calcul des coefficients (rho0,E0,Sr) en fonction de dPdT gamma et liquide
    err=calcul_coeff(dPdT12, dPdT13, pcoeffgamma, pcoeffliq);

    err=hugoniot(pcoeffbeta, pcoeffgamma, pcoeffliq, fhug);
    fclose(fhug);


    // #################################################################
    //    calcul valeur convergés des dPdT pour un point triple optimisé
    printf("\n\n####################################################################################\n");
    printf("            Calcul valeur des dP/dT optimisé sur le point triple \n");
    printf("####################################################################################\n");

    dPdT12= -21906759.5996373; dPdT13= 31164733.1363037;

    if((fhug = fopen("fichiers/hug2.txt", "w+")) == NULL){printf("erreur ouverture fichier hugoniot (recalage_coeff_dPdT.c)\n"); return 1;}
    
    // Calcul des coefficients (rho0,E0,Sr) en fonction de dPdT gamma et liquide
    err=calcul_coeff(dPdT12, dPdT13, pcoeffgamma, pcoeffliq);

    err= hugoniot(pcoeffbeta, pcoeffgamma, pcoeffliq, fhug);
    fclose(fhug);
  */


  return 0;
}


// ***************************************************************************************************************************************************************************************************************
//                      Fonctions thermodynamiques dependantes des coeff (rho0,E0,Sr) et de (V,E)
// ***************************************************************************************************************************************************************************************************************


//  *******  Fonction qui calcule les points (V,E) des lignes de chgmt de phases
/*
Arguments d'entree :
  - Nmax : nombre max d'iterations
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - [Pdeb, Pfin] : segment de presion sur lequel calcule la zone mixte
  - Np : nomnbre de points de discretisatoin du segment
// *********    fonction d'inversion    ************
/*
    Obtention de (V,E) à partir de (P,T)
    NEwton classique 
*/
int fVE(int phase, double rho0, double E0, double P, double T, double* pV, double* pE){
  int Nmax=1e2;
  double epsilon=1e-15;
  int n=0;
  double f, fprime;
  double V,E, dV;
  
  double K0, N0, gamma0, Crv, theta0, ur;
  coeff(phase, &K0, &N0, &gamma0, &Crv, &theta0, &ur);

  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  VaT=0.000132108; EaT=59228.9;
  VbT=0.0001303;   EbT=90321;
  VcT=0.00013427;  EcT=120171;
  PT=3.13049e+09;  TT=582.243;

  // valeur initiale  (point dans la phase en l'occurence un sommet du triangle triple)
  if(phase==1){
    V=VaT;
    E=EaT;
  }
  else if(phase==2){
    V=VbT;
    E=EbT;
  }
  else if(phase==3){
    V=VcT;
    E=EcT;
  }
  else{printf("erreur phase= %d (fVE/recalage_coeff_dPdT.c)\n",phase );}

  
  double critere=1.;
  
  // Newton sur E
  while(critere>epsilon && n<Nmax){
    f=fPs(rho0,V, phase) + gamma0*rho0*Crv*(T-ur*THETA(rho0,V, phase))-P;
    fprime=Ps_prime(rho0,V,phase) - Crv*ur*gamma0*gamma0*rho0*rho0*exp(gamma0*fx(rho0,V));
    
    dV=f/fprime;
    V-=dV;
    critere=fabs(dV/V);
    //printf("  n= %d, V= %.15g, E= %.15g  critere=%.3g\n",n,V,fEs(rho0,E0,V,phase) + (P-fPs(rho0,V,phase))/(rho0*gamma0),critere );
    n++;
  }
  E = fEs(rho0,E0,V,phase) + (P-fPs(rho0,V,phase))/(rho0*gamma0);

  *pV=V;
  *pE=E;

  return 0;
}



void coeff(int phase, double* K0, double* N0, double* gamma0, double* Cvr, double* theta0, double* ur0){
  if(phase==1){
    *gamma0=gamma01; *Cvr=Cv1; *theta0=theta01; *K0=K01; *N0=N01; *ur0=ur;
  }
  else if(phase==2){
    *gamma0=gamma02; *Cvr=Cv2; *theta0=theta02; *K0=K02; *N0=N02; *ur0=ur;
  }
  else if(phase==3){
    *gamma0=gamma03; *Cvr=Cv3; *theta0=theta03; *K0=K03; *N0=N03; *ur0=ur;
  }
  else{printf("Erreur choix phase (coeff)\n");}
}


// x
double fx(double rho0, double V){
  return 1-rho0*V;
}

// THETA
double THETA(double rho0,double V, int phase){
  double gamma0, theta0;
  if(phase==1){
    gamma0=gamma01; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; theta0=theta03;
  }
  else{printf("Erreur choix phase (THETA)\n");}

  
  return theta0*exp( gamma0*fx(rho0,V) );
}

// Es
double fEs(double rho0,double E0,double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (fEs)\n");}

  double x=fx(rho0,V);

  double res=E0;
  double fac=K0/(rho0*(N0+1.0));
  res+= fac*( exp( (N0+1.0)*x ) -1.0 )/(N0+1.0);
  res-= fac*x;

  return res;
}


// Ps
double fPs(double rho0,double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (fPs)\n");}


  double fac=K0/(N0+1);
  return fac*( exp( (N0+1)*fx(rho0,V) ) -1.0 );
}

// u
double u(double rho0, double E0, double V, double E, int phase){
  double Cvr;
  if(phase==1){
  	Cvr=Cv1;
  }
  else if(phase==2){
  	Cvr=Cv2;
  }
  else if(phase==3){
  	Cvr=Cv3;
  }
  else{printf("Erreur choix phase (u)\n");}

  double res=ur;
  res+=( E-fEs(rho0,E0,V,phase) )/( Cvr*THETA(rho0,V,phase) ); 

  return res;
}

// S
double fS(double rho0, double E0, double Sr, double V, double E, int phase){
  double Cvr;
  if(phase==1){
    Cvr=Cv1;
  }
  else if(phase==2){
    Cvr=Cv2;
  }
  else if(phase==3){
    Cvr=Cv3;
  }
  else{printf("Erreur choix phase (fS)\n");}

  double res=Sr;
  res+= Cvr*log(u(rho0,E0,V,E,phase));

  return res;
}


// P
double fP(double rho0,double E0, double V, double E, int phase){
  double gamma0;
  if(phase==1){
    gamma0=gamma01; 
  }
  else if(phase==2){
    gamma0=gamma02;
  }
  else if(phase==3){
    gamma0=gamma03;
  }
  else{printf("Erreur choix phase (fP)\n");}

  double res=fPs(rho0,V,phase);
  res+= gamma0*rho0*( E - fEs(rho0,E0,V,phase) );

  return res;
}


// T
double fT(double rho0,double E0, double V, double E, int phase){
  double Cvr;
  if(phase==1){
    Cvr=Cv1;
  }
  else if(phase==2){
    Cvr=Cv2;
  }
  else if(phase==3){
    Cvr=Cv3;
  }
  else{printf("Erreur choix phase (fT)\n");}

  double res=(E-fEs(rho0, E0, V,phase))/Cvr;
  res+= ur*THETA(rho0,V,phase);

  return res;
}


// Es'
double Es_prime(double rho0, double V, int phase){
  /*
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (Es_prime)\n");}

  double x=fx(rho0,V);

  return (K0/(N0+1))*( 1 - exp( (N0+1)*x ) );
  */
  return -fPs( rho0, V,  phase);

}

// Es''
double Es_prime2(double rho0, double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(rho0,V);

  return rho0*K0*exp( (N0+1)*x );
}

// Es'''
double Es_prime3(double rho0, double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(rho0,V);

  return -rho0*rho0*K0*(N0+1)*exp( (N0+1)*x );
}

// Ps'
double Ps_prime(double rho0, double V, int phase){
  return - Es_prime2(rho0,V, phase);
}

// Ps''
double Ps_prime2(double rho0, double V, int phase){
  return - Es_prime3(rho0, V, phase);
}

double THETA_prime(double rho0, double V, int phase){
  double theta0, gamma0;
  if(phase==1){
    theta0=theta01; gamma0=gamma01;
  }
  else if(phase==2){
    theta0=theta02; gamma0=gamma03;
  }
  else if(phase==3){
    theta0=theta03; gamma0=gamma03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(rho0,V);

  return -theta0*gamma0*rho0*exp(gamma0*x);
}

double THETA_prime2(double rho0, double V, int phase){
  double theta0, gamma0;
  if(phase==1){
    theta0=theta01; gamma0=gamma01;
  }
  else if(phase==2){
    theta0=theta02; gamma0=gamma03;
  }
  else if(phase==3){
    theta0=theta03; gamma0=gamma03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(rho0,V);

  return gamma0*gamma0*rho0*rho0*theta0*exp(gamma0*x);
}

// dP/dE
double dPdE(double rho0, int phase){
  double gamma0;
  if(phase==1){
    gamma0=gamma01;
  }
  else if(phase==2){
    gamma0=gamma02;
  }
  else if(phase==3){
    gamma0=gamma03;
  }
  else{printf("Erreur choix phase (dPdE)\n");}

  return gamma0*rho0;
}


// dP/dV
double dPdV(double rho0, double V, int phase){
  double gamma0;
  if(phase==1){
    gamma0=gamma01;
  }
  else if(phase==2){
    gamma0=gamma02;
  }
  else if(phase==3){
    gamma0=gamma03;
  }
  else{printf("Erreur choix phase (dPdV)\n");}

  return -Es_prime2(rho0,V,phase) - gamma0*rho0*Es_prime(rho0,V,phase);
}

// dT/dE
double dTdE(int phase){
  double Cvr;
  if(phase==1){
    Cvr=Cv1; 
  }
  else if(phase==2){
    Cvr=Cv2; 
  }
  else if(phase==3){
    Cvr=Cv3; 
  }
  else{printf("Erreur choix phase (dTdE)\n");}

  return 1./Cvr;
}


// dTdV
double dTdV(double rho0, double V, int phase){
  double gamma0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (dTdV)\n");}
  
  return (-1./Cvr)*Es_prime(rho0,V,phase) + ur*THETA_prime(rho0,V, phase);

}


// d2P/dE2
double d2PdE2(){
  return 0;
}


// d2P/dV2
double d2PdV2(double rho0, double V, int phase){
  double gamma0;
  if(phase==1){
    gamma0=gamma01;
  }
  else if(phase==2){
    gamma0=gamma02;
  }
  else if(phase==3){
    gamma0=gamma03;
  }
  else{printf("Erreur choix phase (dPdV)\n");}

  return -Es_prime3(rho0,V,phase) - gamma0*rho0*Es_prime2(rho0,V,phase);
}

// d2T/dE2
double d2TdE2(){
  return 0;
}


//    d2T/dV2
double d2TdV2(double rho0, double V, int phase){
  double gamma0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (dTdV)\n");}
  
  return (-1./Cvr)*Es_prime2(rho0,V,phase) + ur*THETA_prime2(rho0,V,phase);
}


double Ps_prime_rho0(double rho0, double V, int phase){
  double gamma0, N0, K0;
  if(phase==1){
    gamma0=gamma01; N0=N01;  K0=K01;
  }
  else if(phase==2){
    gamma0=gamma02; N0=N02;  K0=K02;
  }
  else if(phase==3){
    gamma0=gamma03; N0=N03;  K0=K03;
  }
  else{printf("Erreur choix phase (Ps_prime_rho0)\n");}
  
  double x=fx(rho0,V);

  return -V*K0*exp((N0+1)*x);
}

double THETA_prime_rho0(double rho0, double V, int phase){
  double theta0, gamma0;
  if(phase==1){
    theta0=theta01; gamma0=gamma01;
  }
  else if(phase==2){
    theta0=theta02; gamma0=gamma03;
  }
  else if(phase==3){
    theta0=theta03; gamma0=gamma03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(rho0,V);

  return -V*gamma0*theta0*exp(gamma0*x);
}
    
double h(double rho0, double V, int phase){
  double K0, N0;
  if(phase==1){
    K0=K01; N0=N01;
  }
  else if(phase==2){
    K0=K02; N0=N02;
  }
  else if(phase==3){
    K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (h)\n");}

  double x=fx(rho0,V);

  double res=0;
  double fac=K0/(rho0*(N0+1.0));
  res+= fac*( exp( (N0+1.0)*x ) -1.0 )/(N0+1.0);
  res-= fac*x;

  return res;
}


double g(double rho0, double E0, double V, double E, int phase){
  double Cv;
  if(phase==1){
    Cv=Cv1;
  }
  else if(phase==2){
    Cv=Cv2;
  }
  else if(phase==3){
    Cv=Cv3;
  }
  else{printf("Erreur choix phase (g)\n");}

  return Cv*log( u(rho0, E0, V, E, phase) );

}