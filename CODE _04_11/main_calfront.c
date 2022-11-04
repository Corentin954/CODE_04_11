#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "head_multi.h"
#include "EOS_etain.h"


time_t time (time_t *temps);
double difftime (time_t temps1, time_t temps2);

// ZONE 
// repere les zones (triangle pour le moment)
//err=diag_zone();
//if(err){return err;}

// Frontière VE
//err=routine_VE();
//if(err){return err;}


//**  MAILLAGE  **
//err=maillage();


// Zone mixte
//mapsPT();

// Ecriture de la tabulation
  
int main(){
  int err;

  // Calcul du point triple 
  /*
  double Va,Vb,Vc;
  double Ea,Eb,Ec;
  double Pt,Tt;
  void point_triple(double* Pt, double* Tt, double* Vat, double* Eat, double* Vbt, double* Ebt, double* Vct, double* Ect);
  point_triple(&Pt, &Tt, &Va, &Ea, &Vb, &Eb, &Vc, &Ec);

  printf("Resultats :\n");
  printf(" Pt= %.10lf  (GPa)\n",Pt*1e-9 );
  printf(" Tt= %.10lf\n",Tt );

  printf("Va=%.5g Ea=%.5g\n",Va,Ea);
  printf("Vb=%.5g Eb=%.5g\n",Vb,Eb);
  printf("Vc=%.5g Ec=%.5g\n",Vc,Ec);
  
  int fVE(int phase, double P, double T, double* pV, double* pE);
  */

  /*
    double Vbeta=0.000120369032083173, Ebeta=51821.1842294349, Sbeta=-58.5725135654394;

    double P02=9.4e9, T02=300;
    double Deltav12=-2.247e-6;
    //double Deltav12=-2.24707791549829e-06;
    double dPdT12=-2.13e7;
    //double dPdT12=-21299251.7349429;

    
    // Point de reference BETA/GAMMA
    double Vgamma = Vbeta + Deltav12;
    double Sgamma = Sbeta + dPdT12*Deltav12;
    //double Egamma = Ebeta + P02*(Vbeta - Vgamma) - T02*(Sbeta - Sgamma);
    double Egamma = Ebeta - P02*Deltav12 + T02*dPdT12*Deltav12;

    //double Vgamma=0.000118121954167674;
    //double Egamma=87301.2315927558;
   
    double Pbeta=fP(Vbeta,Ebeta,1);
    double Tbeta=fT(Vbeta,Ebeta,1);
    double Pgamma=fP(Vgamma,Egamma,2);
    double Tgamma=fT(Vgamma,Egamma,2);

    printf("V_beta=%g (rho_beta=%g), S_beta=%g, E_beta=%g\n", Vbeta,1./Vbeta, Sbeta, Ebeta );
    printf("V_gamma=%.10g (rho_gamma=%g), S_gamma=%g, E_gamma=%.10g\n", Vgamma, 1./Vgamma, Sgamma, Egamma );
    printf("Pbeta=%g Tbeta=%g\n", Pbeta, Tbeta );
    printf("Pgamma=%g Tgamma=%g\n", Pgamma, Tgamma );
    printf("(Pgamma-P02)/P02=%g  (Tgamma-T02)/T02=%g\n", (Pgamma-P02)/P02, (Tgamma-T02)/T02 );
     
     
    double V0=1./7287, E0=0;
    double V1=1./7301.02, E1=43080.8;
    double V_beta, E_beta;


    printf("\n");
    printf("Point de reference beta  : V0=%g E0=%g - P=%g  T=%g\n",V0,E0,fP(V0,E0,1),fT(V0,E0,1));
    printf("Point de reference gamma : V1=%g E1=%g - P=%g  T=%g\n",V1,E1,fP(V1,E1,2),fT(V1,E1,2));
    
    printf("\n");
    double P0=0, T0=300;
    printf("Newton phase beta\n");
    err=fVE(1, P0, T0, &V_beta, &E_beta);    if(err){printf("erreur\n"); return 1;}
    printf("P0=%g T0=%g - V=%g (rho=%g) E=%g\n",P0,T0,V_beta,1./V_beta,E_beta );
   
    double T0=300, P0=9.4e9;
    double V_beta, E_beta;
    double V_gamma, E_gamma;
    double V_liq, E_liq;  
    double S_beta,S_gamma,S_liq;

    printf("\n--------------------------------------------------\n");
    printf("Point de ref (P=%g T=%g)    BETA/GAMMA\n",P0,T0);

    printf("Newton phase beta\n");
    err=fVE(1, P0, T0, &V_beta, &E_beta);    if(err){printf("erreur\n"); return 1;}
    printf("V_beta=%.15g  E_beta=%.15g rho_beta=%.15g\n",V_beta,E_beta,1./V_beta );

    printf("\nNewton phase gamma\n");
    err=fVE(2, P0, T0, &V_gamma, &E_gamma);  if(err){printf("erreur\n"); return 1;}
    printf("V_gamma=%.15g  E_gamma=%.15g rho_gamma=%.15g\n",V_gamma,E_gamma,1./V_gamma );
    
    printf("Delta_V=V_gamma - V_beta=%.15g  \n",V_gamma - V_beta );


    printf("\n");
    S_beta=fS(V_beta,E_beta,1);  S_gamma=fS(V_gamma,E_gamma,2);
    printf("S_beta=%.15g   S_gamma=%.15g\n",S_beta,S_gamma );
    printf("Clapeyron : \n(S_gamma - S_beta)/(V_gamma - V_beta) = %.15g\n",(S_gamma - S_beta)/(V_gamma - V_beta) );
    
    // Verif G_beta = G_gamma
    printf("\n -----------------------------------------------------------------------\n");
    double G_beta,G_gamma;
   
    double Prel0, coeffP=1e9;
    double Vrel_gamma, Vrel_beta, coeffV=1e-6;
    Prel0=P0/coeffP;
    Vrel_beta=V_beta/coeffV;
    Vrel_gamma=V_gamma/coeffV;
    G_beta= E_beta + (coeffV*coeffP)*Prel0*Vrel_beta - T0*S_beta;
    G_gamma=E_gamma + (coeffV*coeffP)*Prel0*Vrel_gamma - T0*S_gamma;
   
    G_beta= E_beta  + P0*V_beta  - T0*S_beta;
    G_gamma=E_gamma + P0*V_gamma - T0*S_gamma;

    printf("(V,E) Newton : Pbeta=%g  Tbeta=%g    Pgamma=%g Tgamma=%g\n",fP(V_beta,E_beta,1),fT(V_beta,E_beta,1),fP(V_gamma,E_gamma,2),fT(V_gamma,E_gamma,2) );
    printf("               Pbeta-Pgamma=%g    Tbeta-Tgamma=%g\n",fP(V_beta,E_beta,1)-fP(V_gamma,E_gamma,2),fT(V_beta,E_beta,1)-fT(V_gamma,E_gamma,2) );

    printf("G_beta=%.15g G_gamma=%.15g  G_beta-G_gamma=%g\n",G_beta,G_gamma,G_beta-G_gamma );
    
    double Etest=E_beta + P0*(V_beta-V_gamma) - T0*(S_beta-S_gamma);
    printf("E (G=G) =%.15g  E_gamma=%.15g  DeltaE=%g\n",Etest,E_gamma,Etest-E_gamma);
    
    printf("\nfonction fG\n");
    printf("fG_beta=%.15g  fG_gamma=%.15g  DeltafG=%g\n",fG(V_beta,E_beta,1),fG(V_gamma,E_gamma,2),fG(V_beta,E_beta,1)-fG(V_gamma,E_gamma,2) );
    

    printf("\n--------------------------------------------------\n");
    printf("Newton ligne double \n");
    
    double Tn0, VbetaN, EbetaN, VgammaN, EgammaN;
    
    void point_double(int phaseA, int phaseB, double P, double* pT, double* pVa, double* pEa, double* pVb, double* pEb);
    point_double(1, 2, P0, &Tn0, &VbetaN, &EbetaN, &VgammaN, &EgammaN);

    printf("\n--------------------------------------------------\n");
    T0=505; P0=0;
    printf("Point de ref (P=%g T=%g)   LIQUIDE/BETA\n",P0,T0);

    printf("\nNewton phase beta\n");
    err=fVE(1, P0, T0, &V_beta, &E_beta);    if(err){printf("erreur\n"); return 1;}
    printf("V_beta=%.15g  E_beta=%.15g rho_beta=%.15g \n",V_beta,E_beta,1./V_beta);

    printf("\nNewton phase liquide\n");
    err=fVE(3, P0, T0, &V_liq, &E_liq);      if(err){printf("erreur\n"); return 1;}
    printf("V_liq=%.15g  E_liq=%.15g rho_liq=%.15g\n",V_liq,E_liq,1./V_liq );

    printf("Delta_V=V_liq - V_beta=%.15g  \n",V_liq - V_beta );
    
    printf("\n");
    S_beta=fS(V_beta,E_beta,1);  S_liq=fS(V_liq,E_liq,3);
    printf("S_beta=%g   S_liq=%g\n",S_beta,S_liq );
    printf("Clapeyron : \n(S_liq - S_beta)/(V_liq - V_beta) = %g\n",(S_liq - S_beta)/(V_liq - V_beta) );

    double rho0_beta=7287.0, rho0_gamma=7401.45, rho0_liq=7026.17;
    double delta12=-1.9e-6, delta13=4e-6;
    
    printf("rho0_gamma-rho0_beta=%g  1/Delta_V=%g\n",rho0_gamma-rho0_beta,1./delta12 );
    printf("rho0_liq-rho0_beta=%g  1/Delta_V=%g\n",rho0_liq-rho0_beta,1./delta13 );
  */


  // ZONE 
  // repere les zones (triangle pour le moment)
  
  //int diag_zone();
  //err=diag_zone();

  //int plot_PT_zone3_domaine_validite();
  //err=plot_PT_zone3_domaine_validite();

  // Ecriture de la tabulation
  
  int nb_points=10;
  
  int Nv_dis=2e2;
  int Ne_dis=2e2;

  /*
  double Vmin_dis=81e-6, Vmax_dis=172e-6;
  double Emin_dis=-7e4, Emax_dis=1.4e6;
  */

  double Vmin_dis=100.e-6, Vmax_dis=120e-6;
  double Emin_dis=0,  Emax_dis=1e6;
  
  int ecr_benchmark_PT(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax);
  err=ecr_benchmark_PT(nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis);
  if(err){return err;}
  

  /*
  double fP(double V, double E, int phase);
  double fT(double V, double E, int phase);

  double P=3e9, T=400;
  double V,E;
  err=fVE(1, P, T, &V, &E);
  printf("P=%g T=%g\n",P,T );
  printf("V=%g E=%g\n",V,E );
  double Pnew=fP(V,E,1); 
  double Tnew=fT(V,E,1); 
  printf("Pnew=%g Tnew=%g\n",Pnew,Tnew );
  */


  /*
  int j, temp;
  int n=10;
  int * L=malloc(n*sizeof(int));
  L[0]=2;   L[1]=3;    L[2]=6;    L[3]=1;  L[4]=10;  
  L[5]=-1;  L[6]=-81;  L[7]=100;  L[8]=4;  L[9]=5;

  for(int i=0; i<n; ++i){ printf(" L[%d]=%d",i,L[i] ); } printf("\n");

  for(int i=1; i<n; ++i){
    printf(" -i=%d\n",i );
    while( L[i-1]<L[i] ){
      j=i;
      while( L[j-1]<L[j] && j>0 ){
        printf("    j=%d\n",j );
        temp=L[j];
        L[j]=L[j-1];
        L[j-1]=temp;
        j--;
      }
    }
  }

  printf(" L sorted \n");
  for(int i=0; i<n; ++i){ printf(" L[%d]=%d",i,L[i] ); } printf("\n");
  */

  

  // USE BENCHMARK PT
  // Discrétisation du test
  /*
  int EOS=2;

  int Nv_test=2e2;
  int Ne_test=2e2;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;
  
  
  Vmin_test=110e-6, Vmax_test=140e-6;
  Emin_test=-1e3,   Emax_test=3e5;

  int use_benchmark_PT(int EOS, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  err=use_benchmark_PT(EOS, Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);
  if(err){return err;}
  */
  
  // TEST POUR LA CIN2TIQUE DES PHASES
  // Discrétisation du test
  /*
  int EOS=0; //hybride + sans détection de zone
  double epsilon=1e-13;
  
  int Nv_test=500;//=400;
  int Ne_test=500;//=400;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;

  // ZOOM
  
  Vmin_test=110e-6, Vmax_test=145e-6;
  Emin_test=-5e3,   Emax_test=3e5;
 
  //Vmin_test=82e-6, Vmax_test=180e-6;
  //Emin_test=-5e3,   Emax_test=9e5;


  // MAX DOMAIN
  //Vmin_test=85e-6, Vmax_test=170e-6;
  //Emin_test=-6.5e4,   Emax_test=8e5;

  time_t t1;
  time_t t2;
  int use_benchmark_PT_cin_phase(int EOS, double epsilon, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  t1 = time(NULL);
  err=use_benchmark_PT_cin_phase(EOS, epsilon, Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);  if(err){return err;}
  t2 = time(NULL);
  
  printf(" Pour EOS=%d, Nv=%d et Ne=%d | Temps de calcul= %ld secondes \n",EOS,Nv_test,Ne_test,t2-t1 );
  */


  // TEST FRACTION MASSIQUE DANS LE TRIANGLE TRIPLE
  /*
  int coorbary(double *xy, double **sommets, double **pcoorbar);
  // Valeur d'initialisation
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  double* Xfm=malloc(3*sizeof(double));

  double* VE=malloc(2*sizeof(double));
  double V= 0.000129443; double E= 88522.7;
  VE[0]=V;   VE[1]=E;

  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;  
  err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
  
  printf(" V= %g, E= %g, xA= %f, xB= %g, xC= %g\n",VE[0],VE[1],Xfm[0],Xfm[1],Xfm[2] );
  */
  
  // Calcul d'une isentrope
  //int isentrope();
  //err=isentrope();


  //********   Calcul d'une courbe d'Hugoniot  *********
  //int hugoniot();
  //err=hugoniot();

  //********   onde_compression_isentropique   ********* 
  //int profil_onde_compression_isentropique();
  //err=profil_onde_compression_isentropique();


  // **********  TEST DES DIFFERENTES METHODES  ******************** 
  // Discrétisation du test
  /*
  int Nv_test=3e2;
  int Ne_test=3e2;

  double Vmin_test, Vmax_test;
  double Emin_test, Emax_test;

  Vmin_test=98e-6, Vmax_test=150e-6;
  Emin_test=-5e3,   Emax_test=4e5;

  int methode=1;
  
  int use_fPTC_METH_cin_phase(int methode, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test);

  err=use_fPTC_METH_cin_phase(methode, Nv_test, Ne_test, Vmin_test, Emin_test, Vmax_test, Emax_test);
  if(err){return err;}
  */

  
  


}