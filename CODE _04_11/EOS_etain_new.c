#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//### Nouveau jeu de coefficients  ######//

// Isotherme Inversible

// Phase 1
double K01=53.597e9;
double N01=4.996;
double gamma01=2.27;
double Cv1=210.;
double theta01=250.;
double T01=300.;
double P01=0.;

double rho01=7287.;   // 1./rho01;
double v01=1./7287.;   // 1./rho01;
double E01=0.; 
double Sr1=-38.2875;
//double F03= 300*38.2875;          // E0 - T0*S0;
 

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

double rho02=7301.02;   // 1./rho01;
double v02=1./7301.02;   // 1./rho01;
double E02=43080.8; 
double Sr2=45.9203;
//double F03=39033. - 300*24.7813;          // E0 - T0*S0;


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

double rho03=7036.77;   // 1./rho01;
double v03=1./7036.77;   // 1./rho01;
double E03=88479.8; 
double Sr3=160.786;
//double F03=90660.5 - 505*165.11;          // E0 - T0*S0;

double ur=1.2;

// Coeff MABIRE
/*
K02=94e9;
N02=4.88;

K03=42e9;
*/




void coeff(int phase, double* K0, double* N0, double* gamma0, double* Cvr, double* theta0, double* T0, double* P0, double* rho0, double* v0, double* E0, double* Sr, double* ur0){
  if(phase==1){
    *gamma0=gamma01; *rho0=rho01; *Cvr=Cv1; *theta0=theta01; *E0=E01; *K0=K01; *N0=N01; *ur0=ur;
  }
  else if(phase==2){
    *gamma0=gamma02; *rho0=rho02; *Cvr=Cv2; *theta0=theta02; *E0=E02; *K0=K02; *N0=N02; *ur0=ur;
  }
  else if(phase==3){
    *gamma0=gamma03; *rho0=rho03; *Cvr=Cv3; *theta0=theta03; *E0=E03; *K0=K03; *N0=N03; *ur0=ur;
  }
  else{printf("Erreur choix phase (coeff)\n");}
}


// x
double fx(double V, int phase){
  double rho0;
  if(phase==1){
    rho0=rho01; 
  }
  else if(phase==2){
    rho0=rho02;
  }
  else if(phase==3){
    rho0=rho03;
  }
  else{printf("Erreur choix phase (fx)\n");}

  return 1-rho0*V;
}

// THETA
double THETA(double V, int phase){
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

  
  return theta0*exp( gamma0*fx(V,phase) );
}

// Es
double fEs(double V, int phase){
  double E0, K0, rho0, N0;
  if(phase==1){
    E0=E01; K0=K01; rho0=rho01; N0=N01;
  }
  else if(phase==2){
    E0=E02; K0=K02; rho0=rho02; N0=N02;
  }
  else if(phase==3){
    E0=E03; K0=K03; rho0=rho03; N0=N03;
  }
  else{printf("Erreur choix phase (fEs)\n");}

  double x=fx(V,phase);

  double res=E0;
  double fac=K0/(rho0*(N0+1.0));
  res+= fac*( exp( (N0+1.0)*x ) -1.0 )/(N0+1.0);
  res-= fac*x;

  return res;
}


// Ps
double fPs(double V, int phase){
  double E0, K0, N0;
  if(phase==1){
    E0=E01; K0=K01; N0=N01;
  }
  else if(phase==2){
    E0=E02; K0=K02; N0=N02;
  }
  else if(phase==3){
    E0=E03; K0=K03; N0=N03;
  }
  else{printf("Erreur choix phase (fPs)\n");}


  double fac=K0/(N0+1);
  return fac*( exp( (N0+1)*fx(V,phase) ) -1.0 );
}

// u
double u(double V, double E, int phase){
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
  res+=( E-fEs(V,phase) )/( Cvr*THETA(V,phase) ); 

  return res;
}

// S
double fS(double V, double E, int phase){
  double Sr, Cvr;
  if(phase==1){
    Sr=Sr1; Cvr=Cv1;
  }
  else if(phase==2){
    Sr=Sr2; Cvr=Cv2;
  }
  else if(phase==3){
    Sr=Sr3; Cvr=Cv3;
  }
  else{printf("Erreur choix phase (fS)\n");}

  double res=Sr;
  res+= Cvr*log(u(V,E,phase));

  return res;
}


// P
double fP(double V, double E, int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (fP)\n");}

  double res=fPs(V,phase);
  res+= gamma0*rho0*( E - fEs(V,phase) );

  return res;
}


// T
double fT(double V, double E, int phase){
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

  double res=(E-fEs(V,phase))/Cvr;
  res+= ur*THETA(V,phase);

  return res;
}

double fG(double V, double E, int phase){
  return E + V*fP(V,E,phase) - fT(V,E,phase)*fS(V,E,phase);
}


// Es'
double Es_prime(double V, int phase){
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

  double x=fx(V,phase);

  return (K0/(N0+1))*( 1 - exp( (N0+1)*x ) );
}

// Es''
double Es_prime2(double V, int phase){
  double K0, N0, rho0;
  if(phase==1){
    K0=K01; N0=N01; rho0=rho01;
  }
  else if(phase==2){
    K0=K02; N0=N02; rho0=rho02;
  }
  else if(phase==3){
    K0=K03; N0=N03; rho0=rho03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(V,phase);

  return rho0*K0*exp( (N0+1)*x );
}

// Es'''
double Es_prime3(double V, int phase){
  double K0, N0, rho0;
  if(phase==1){
    K0=K01; N0=N01; rho0=rho01;
  }
  else if(phase==2){
    K0=K02; N0=N02; rho0=rho02;
  }
  else if(phase==3){
    K0=K03; N0=N03; rho0=rho03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(V,phase);

  return -rho0*rho0*K0*(N0+1)*exp( (N0+1)*x );
}

// Ps'
double Ps_prime(double V, int phase){
  return - Es_prime2(V, phase);
}

// Ps''
double Ps_prime2(double V, int phase){
  return - Es_prime3(V, phase);
}

double THETA_prime(double V, int phase){
  double theta0, gamma0, rho0;
  if(phase==1){
    theta0=theta01; gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    theta0=theta02; gamma0=gamma03; rho0=rho02;
  }
  else if(phase==3){
    theta0=theta03; gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(V,phase);

  return -theta0*gamma0*rho0*exp(gamma0*x);
}

double THETA_prime2(double V, int phase){
  double theta0, gamma0, rho0;
  if(phase==1){
    theta0=theta01; gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    theta0=theta02; gamma0=gamma03; rho0=rho02;
  }
  else if(phase==3){
    theta0=theta03; gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (Ps_prime)\n");}

  double x=fx(V,phase);

  return gamma0*gamma0*rho0*rho0*theta0*exp(gamma0*x);
}

// dP/dE
double dPdE(int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (dPdE)\n");}

  return gamma0*rho0;
}


// dP/dV
double dPdV(double V, int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (dPdV)\n");}

  return -Es_prime2(V,phase) - gamma0*rho0*Es_prime(V,phase);
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
double dTdV(double V, int phase){
  double gamma0, rho0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (dTdV)\n");}
  
  return (-1./Cvr)*Es_prime(V,phase) + ur*THETA_prime(V, phase);

}


// d2P/dE2
double d2PdE2(){
  return 0;
}


// d2P/dV2
double d2PdV2(double V, int phase){
  double gamma0, rho0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03;
  }
  else{printf("Erreur choix phase (dPdV)\n");}

  return -Es_prime3(V,phase) - gamma0*rho0*Es_prime2(V,phase);
}

// d2T/dE2
double d2TdE2(){
  return 0;
}


//    d2T/dV2
double d2TdV2(double V, int phase){
  double gamma0, rho0, Cvr, theta0;
  if(phase==1){
    gamma0=gamma01; rho0=rho01; Cvr=Cv1; theta0=theta01;
  }
  else if(phase==2){
    gamma0=gamma02; rho0=rho02; Cvr=Cv2; theta0=theta02;
  }
  else if(phase==3){
    gamma0=gamma03; rho0=rho03; Cvr=Cv3; theta0=theta03;
  }
  else{printf("Erreur choix phase (dTdV)\n");}
  
  return (-1./Cvr)*Es_prime2(V,phase) + ur*THETA_prime2(V,phase);
}


// EPSILON fonction de V et P  (pour l'initailisation dans les schemas)
double fE_VP(double V, double P, int phase){
  double gamma0, rho0;
  if(phase==1){  gamma0=gamma01; rho0=rho01;  }
  else if(phase==2){  gamma0=gamma02; rho0=rho02;  }
  else if(phase==3){  gamma0=gamma03; rho0=rho03;  }
  else{printf("Erreur choix phase (fE_VP)\n");}

  return fEs(V, phase) - (P-fPs(V, phase))/(gamma0*rho0);
}

