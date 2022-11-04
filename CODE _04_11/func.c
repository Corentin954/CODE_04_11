#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "head_multi.h"
#include "EOS_etain.h"
#include "sol_head.h"


/* Retourne la valeur du pas de temps

*/
double pastemps(double T, double t, double dt, double dt_old){
  double res;
  if(dt>1.1*dt_old){
    dt=1.1*dt_old;
  }
  if ((T-t)<dt){ dt=T-t; }
  return dt;
}

// ind_u=0 : vitesse U défini sur maille primal
// ind_u=1 : vitesse U défini sur maille duale
int print_sol(int ind_u, int ideb, int ifin, double* X, double* Xc, double* RHO, double* TAU, double* U, double* E,  double* EPS, double* P, double* T, double* C, double* G, double* S, double** matFM){
  
  FILE* Res;

  // Pour Octave
  // ecriture dans un fichier des points (x,u(x)) de la solution 
  /*
  L1 : Xd  maillage décalée
  L2 : Xc  maillage centrée
  L3 : rho  densité
  L4 : tau  volume spécifique
  L5 : u   vitesse
  L6 : e   energie totale
  L7 : eps   energie interne
  L8 : P   pression
  L9 : T   temperature
  L10 : C   celerite du son
  L11 : G   derivee fondamentale
  L12 : S   entropie
  */
  if((Res = fopen("sol.txt", "w+")) != NULL){
    // x
    for(int i=ideb; i<=ifin+1; i++){ 
      fprintf(Res,"%.15e ",X[i]);
    }
    fprintf(Res, "\n\n");
    // X duale
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ", Xc[i]);
    }
    fprintf(Res, "1e10\n\n");
    // rho
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",RHO[i]);
    }
    fprintf(Res, "1e10\n\n");
    // tau
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",TAU[i]);
    }
    fprintf(Res, "1e10\n\n");
    // u
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",U[i]);
    }
    if(ind_u==0){
      fprintf(Res,"%.15e \n\n",U[ifin+1]);
    }
    else if(ind_u==1){
      fprintf(Res, "1e10\n\n");
    }
    // e
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",E[i]);
    }
    fprintf(Res, "1e20\n\n");
    // epsilon
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15e ",EPS[i]);
    }
    fprintf(Res, "1e20\n\n");
    // P
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",P[i]);
    }
    fprintf(Res, "1e20\n\n");
    // T
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",T[i]);
    }
    fprintf(Res, "1e20\n\n");
    // C
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",C[i]);
    }
    fprintf(Res, "1e20\n\n");
    // G
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15e ",G[i]);
    }
    fprintf(Res, "1e20\n\n");
    // S
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15e ",S[i]);
    }
    fprintf(Res, "1e20\n\n");
    fclose(Res);
  }
  else{
    printf("Erreur d'ouverture de sol.txt (affiche_sol)\n");
    return 1;
  }

  //FILE* Res;
  /*
  L1 : xA
  L2 : xB
  L3 : xC
  */
  if((Res = fopen("frac_mass.txt", "w+")) == NULL){ printf("erreur impression\n"); return 1;}
  // xA
  for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%g ",matFM[i][0]);  }
  fprintf(Res, "1e20\n \n");
  // xB
  for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%g ",matFM[i][1]);  }
  fprintf(Res, "1e20\n \n");
  // xC
  for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%g ",matFM[i][2]);  }
  fprintf(Res, "1e20\n \n");
  fclose(Res);
  
  
  // Pour gnuplot
  // Ecriture dans un fichier des points (x,u(x)) de la solution sur grille centrée
  // si le schema est centrée i.e u est calculé au centres des mailles, on écrit u dans le fichier,
  // sinon on le remplace par des 0 (tout les autres types de schémas)
  /*    
       output_c.txt                        output_d.txt                    output_lambda.txt
  L1  : Xc  maillage centrée          L1  : Xc ou Xd  maillage          L1  : Xc  maillage centrée
  L2  : rho  volume spécifique        L2  : u vitesse                   L2  : fraction massique de la phase beta
  L3  : tau  volume spécifique                                          L3  : fraction massique de la phase gamma
  L4  : e   energie totale                                              L4  : fraction massique de la phase liquide
  L5  : eps   energie interne
  L6  : P   pression
  L7  : T   temperature
  L8  : C   celerite du son
  L9  : G   derivee fondamentale
  L10 : S   entropie
  */
  FILE *Res1,*Res2,*Res3;
  if((Res1 = fopen("output_c.txt", "w+")) == NULL){printf("Erreur ouverture output_c (func.c/print_sol)\n");  return 1;}
  if((Res2 = fopen("output_d.txt", "w+")) == NULL){printf("Erreur ouverture output_d (func.c/print_sol)\n");  return 1;}
  if((Res3 = fopen("output_lambda.txt", "w+")) == NULL){printf("Erreur ouverture output_lambda (func.c/print_sol)\n");  return 1;}
  
  FILE *ResXd;
  if((ResXd = fopen("output_Xd.txt", "w+")) == NULL){printf("Erreur ouverture output_c (func.c/print_sol)\n");  return 1;}
  
  // variables centrées aux mailles
  for(int i=ideb; i<=ifin; i++){
    fprintf(Res1,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",Xc[i],RHO[i],TAU[i],E[i],EPS[i],P[i],T[i],C[i],G[i],S[i]);
    fprintf(Res3,"%.15e %.15e %.15e %.15e\n",Xc[i],matFM[i][0],matFM[i][1],matFM[i][2]);
  }
  // vitesse
  if(ind_u==1){
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res2,"%.15e %.15e\n",Xc[i],U[i]);
    }
  }
  else if(ind_u==0){
    for(int i=ideb; i<=ifin+1; i++){ 
      fprintf(Res2,"%.15e %.15e\n",X[i],U[i]);
    }
  }

  for(int i=ideb; i<=ifin+1; i++){
    fprintf(ResXd,"%.15e\n",X[i]);
  }

  fclose(Res1);  fclose(Res2);  fclose(Res3);  fclose(ResXd);

  return 0;
}

int print_sol_RKint(int ind_u, int ideb, int ifin, double* X, double* Xc, double* RHO, double* TAU, double* U, double* E,  double* EPS, double* P, double* T, double* C, double* G, double* S, double* RTAU, double* RU, double* REPS){
  // ecriture dans un fichier des points (x,u(x)) de la solution
  FILE* Res;
  /*
  L1 : X  maillage
  L2 : tau  volume spécifique
  L3 : u  vitesse
  L4 : e   energie interne
  L5 : p   pression
  L6 : epsilon   energie interne spécifique
  */

  if((Res = fopen("sol.txt", "w+")) != NULL){
    // x
    for(int i=ideb; i<=ifin+1; i++){ 
      fprintf(Res,"%.15e ",X[i]);
    }
    fprintf(Res, "\n\n");
    // X duale
    for (int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ", Xc[i]);
    }
    fprintf(Res, "1e10\n\n");
    // rho
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",RHO[i]);
    }
    fprintf(Res, "1e10\n\n");
    // tau
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",TAU[i]);
    }
    fprintf(Res, "1e10\n\n");
    // u
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",U[i]);
    }
    if(ind_u==0){
      fprintf(Res,"%.15e \n\n",U[ifin+1]);
    }
    else if(ind_u==1){
      fprintf(Res, "1e10\n\n");
    }
    // e
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",E[i]);
    }
    fprintf(Res, "1e20\n\n");
    // epsilon
    for(int i=ideb; i<=ifin; i++){ 
      fprintf(Res,"%.15e ",EPS[i]);
    }    
    fprintf(Res, "1e20\n\n");
    // P
    for(int i=ideb; i<=ifin; i++){
      fprintf(Res,"%.15e ",P[i]);
    }
    fprintf(Res, "1e20\n\n");
    // T
    for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15e ",T[i]);  } 
    fprintf(Res, "1e20\n\n");
    // C
    for(int i=ideb; i<=ifin; i++){   fprintf(Res,"%.15e ",C[i]);  } 
    fprintf(Res, "1e20\n\n");
    // G
    for(int i=ideb; i<=ifin; i++){   fprintf(Res,"%.15e ",G[i]);  }
    fprintf(Res, "1e20\n\n");
    // S
    for(int i=ideb; i<=ifin; i++){   fprintf(Res,"%.15e ",S[i]);  }
    fprintf(Res, "1e20\n\n");

    // RTAU
    for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15e ",RTAU[i]);  }
    fprintf(Res, "1e10\n\n");
    // RU
    for(int i=ideb; i<=ifin+1; i++){  fprintf(Res,"%.15e ",RU[i]);  }
    fprintf(Res, "\n\n");
    // REPS
    for(int i=ideb; i<=ifin; i++){  fprintf(Res,"%.15e ",REPS[i]);  }
    fprintf(Res, "1e20\n\n");

    fclose(Res);
    return 0;
  } 
  else{
    printf("Erreur d'ouverture de sol.txt (affiche_sol)\n");
    return 1;
  }
}


int print_nrj(int n, double* Etot, double* IMPUL, double* Mtot, double* VOL){
  FILE* Res;
  
  // Conservation au cours du temps
  if((Res = fopen("NRJ.txt", "w+")) != NULL){
    // energie totale
    for(int i=0; i<=n-1; i++){
      fprintf(Res,"%.15e ",Etot[i]);
    }
    fprintf(Res, "\n");
    // impulsion
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15e ",IMPUL[i]);
    }
    fprintf(Res, "\n");
    // masse totale
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15e ",Mtot[i]);
    }
    fprintf(Res, "\n");
    // volume total
    for(int i=0; i<=n-1; i++){ 
      fprintf(Res,"%.15e ",VOL[i]);
    }
    fprintf(Res, "\n");
    fclose(Res);
    return 0;
  } 
  else{
    printf("Erreur d'ouverture de NRJ.txt (affiche_nrj)\n");
    return 2;
  }
}


  

// Met a jour les condition aux limites pour la variables Z
/* Si Z est centré aux maille  jfin= ifin
      Z est frontières aux mailles jfin=ifin+1

  int iu = 0 si cond sym
      iu = 1 si cond anti-sym
*/
int condlim(int ind_cond, int jdeb,int jfin, int nbghosts, double* Z, int ind){
  if(ind_cond==0){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){  // U est frontière aux mailles, la valeur à la frontière n'intervient pas
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+1+j];
        // à droite
        Z[jfin+1+j]=-Z[jfin-1-j];
      }
    }
    return 0;
  }
  // Condition périodique 
  else if(ind_cond==1){
    if(ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin-j];
        // à droite
        Z[jfin+1+j]=Z[jdeb+j];
      }
    }
    else if(ind==1){ // pour l'energie (variable au carré)
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin+1+j];
        // à droite
        Z[jfin+1+j]=Z[jdeb-1-j];
      }
    }
    return 0;
  }
  // Condition aux limites pour l'etain : Mur à gauche et libre à droite
  else if(ind_cond==2){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){  // U est frontière aux mailles, la valeur à la frontière n'intervient pas
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+1+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-1-j];
      }
    }
    return 0;
  }
  else{
    printf("Erreur choix de ind_cond (condlim)\n");
    return 5;
  }
}

// Met a jour les condition aux limites pour la variables Z
// Différent de condlim car U est centré aux mailles
/* Si Z est centré aux maille  jfin= ifin
      Z est frontières aux mailles jfin=ifin+1

  int iu = 0 si cond sym
      iu = 1 si cond anti-sym
*/
int condlim_Gtot(int ind_cond, int jdeb,int jfin, int nbghosts, double* Z, int ind){
  //nbghosts=+3;
  if(ind_cond==0){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=-Z[jfin-j];
      }
    }
    return 0;
  }
  // Onde Acoustique 
  else if(ind_cond==1){
    if(ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin-j];
        // à droite
        Z[jfin+1+j]=Z[jdeb+j];
      }
    }
    else if(ind==1){ // pour l'energie (variable au carré)
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jfin+j];
        // à droite
        Z[jfin+1+j]=Z[jdeb-j];
      }
    }
    return 0;
  }
  if(ind_cond==2){
    if (ind==0){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    else if (ind==1){
      for(int j=0; j<nbghosts; j++){
        // à gauche
        Z[jdeb-1-j]=-Z[jdeb+j];
        // à droite
        Z[jfin+1+j]=Z[jfin-j];
      }
    }
    return 0;
  }
  else{
    printf("Erreur choix ind_cond (condlim)\n");
    return 5;
  }
}


// The Minmod limiter
double MinMod(double r){
  return fmax(0.,fmin(1.,r));
}

// The Superbee limiter
double Superbee(double r){
  return fmax(0., fmax(fmin(1.,2.*r),fmin(2.,r)) );
}

// The Christensen limiter
double Christensen(double r, double arg){
  return fmax(0., fmin(2.,fmin(2*r,arg)) );
}

// The van Lerr limiter
double vanLeer(double r_p, double r_m){
  return fmax(0., fmin(1., fmin(2*r_m, fmin(2*r_p, (r_m+r_p)/2))) );
} 


double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM){
  int third=q%10;               // pseudo-viscosity
  int second=(q%100-third)/10;  // limiteur pour l'ordre 2
  int first=q/100;              // ordre


  double du=U[1]-U[0];
  // ordre 1
  //   first==0
  // Ordre 2
  if (first==1){
    // delta vitesse
    double dUm=U[0]-U[-1];
    double dUp=U[2]-U[1];
    // delta masse
    double dM=DM[0];
    double dMm=DM[-1];
    double dMp=DM[1];

    double rplus=(dUp/dMp)*(dM/du);
    double rmoins=(dUm/dMm)*(dM/du);
    
    double phiplus;
    double phimoins;

    if(second==0){
      phiplus=MinMod(rplus);
      phimoins=MinMod(rmoins);
    }
    else if(second==1){
     phiplus=Superbee(rplus);
     phimoins=Superbee(rmoins);
    }
    else if(second==2){
      double argplus=(dM*rplus+dMp)/(dM+dMp);
      double argmoins=(dM*rmoins+dMm)/(dM+dMm);
      phiplus=Christensen(rplus,argplus);
      phimoins=Christensen(rmoins, argmoins);
    }
    else{ printf("Erreur q= %d (qvis)\n",q );  }

    du*=(1-0.5*(phiplus + phimoins));
  }

  
  // von Neumann-Ritchmyer
  if (third == 0){
    return  -Cq*(du)*fabs(du)/tau;
  }
  // Rosenbluth
  else if(third == 1){
  	if (du<0){  return -Cq*(du)*fabs(du)/tau;	}
  	else{ return 0.; }
  }
  // Landschoff
  else if(third == 2){
  	if (du<0){  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;	}
  	else{ return 0; }
  }
  // "Magical" q combination
  else if(third == 3){
  	if (du<0){  return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;	}
  	else{ return -Cl*c*du/tau; }
  }
  // Landschoff
  else if(third == 4){  
    return -Cq*(du)*fabs(du)/tau -Cl*c*du/tau;  
  }
  else if(third == 5){  
    return 0.;  
  }
  else{  printf("Erreur dans le choix de pseudo q= %d (qvis)\n",q);  }
}

// Evaluation de la moyenne (Average) ou ponctuelle (valable jusqu'a l'ordre 4) de la Z
// Le type d'evaluation depend des arguments de la focntion
// k : spatial order
// si . Coef=Ck retourne Point-Wise
//    . Coef=Ckbar retourne Average
// Z = &Z[i]
double phi(int so, double* Z, double** Coef){
  double res=0.0;
  int k=so-2;
  for(int j=-k; j<=k; j++){
    res+=Coef[k][abs(j)]*Z[j];
  }
  return res;
}


// Même principe que la fonction au dessus renvoie la moyenne 
// sauf qu'elle prend en arguments les tableaux : rho.Z et rho 
// il faut faire leur division pour obtenir Z
double phidiv(int so, double* RZ, double* RHO, double** Ckbar){
  double res=0.0;
  int k=so-2;
  for(int j=-k; j<=k; j++){
    res+=Ckbar[k][abs(j)]*(RZ[j]/RHO[j]);
  }
  return res;
}

// Calcul la maille dual (centrée) à partir de la maille primal (frontière) 
double fxc(int so,double* X, double** rk){
  // voir CRAS 2016 p214
  double res=0.0;
  int k=so-2;
  for(int l=0; l<=k; l++){
    res+=rk[k][l]*(X[l+1] + X[-l]);
  }
  return res;
}

// Pour les variables centrées aux mailles :   ind=0
//   delta i+1/2 = Z[i]-Z[i-1]
// Pour les variables frontières aux mailles :   ind=1
//   delta i = Z[i+1]-Z[i]
double delta(int so,int ind, double* Z, double** dk){
  double dz=0.0;
  int k=so-2;
  if (ind==0){
    for(int l=0; l<=k; l++){
      //dz+=dk[k][l]*(Z[l+1] - Z[l]);
      dz+=dk[k][l]*(Z[l] - Z[-l-1]);
    }
  }
  else if (ind==1){
    for(int l=0; l<=k; l++){
      //dz+=dk[k][l]*(Z[l+1] - Z[l]);
      dz+=dk[k][l]*(Z[l+1] - Z[-l]);
    }
  }
  return dz;
}

// Pour rho.eps  , retourne AV(PI.\delta u)
// Retourne la moyenne de YdZ  
// soit : . Y=PI et Z=u  pour EPS
//        . Y=u  et Z=PI pour EKIN
// int    cas==0 : EPS    cas==1 : EKIN
double YdZ(int so, int cas, double* P, double* Q, double* u, double** dk, double** Ckbar){
  double du, res=0.0;
  int k=so-2;
  if(cas==0){ // EPS+=(dt/dx).PIdu
    for(int j=-k; j<=k; j++){
      // calcul de du
      du=0.0;
      for(int l=0; l<=k; l++){
        du+=dk[k][l]*(u[j+(l+1)]-u[j-l]);
      }
      // PI.du
      res+=Ckbar[k][abs(j)]*(P[j]+Q[j])*du;
    }
  }
  else if(cas==1){ // EKIN+=(dt/dx).udPI
    for(int j=-k; j<=k; j++){
      // calcul de dPI
      du=0.0;
      for(int l=0; l<=k; l++){
        du+=dk[k][l]*( (P[j+l]+Q[j+l])-(P[j-l-1]+Q[j-l-1]) );
      }
      // u.dPI
      res+=Ckbar[k][abs(j)]*u[j]*du;
    }
  }
  return res;
}

// Fonction uniquement pour le kinetic energy fix  
// sum Qbar.Delta K
// k : spatial order
// Z = &Z[i]
double phiQK(int so, double* DK, double** Q){
  double res=0.0;
  int k=so-2;
  for(int j=0; j<=k; j++){
    res+=Q[k][j]*(DK[j+1]+DK[-j]);
  }
  return res;
}

void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk){
  // Ck
  Ck[0][0]=1.0;               Ck[0][1]=0.0;               Ck[0][2]=0.0;             Ck[0][3]=0.0;           Ck[0][4]=0.0; 
  Ck[1][0]=13./12;            Ck[1][1]=-1./24;            Ck[1][2]=0.0;             Ck[1][3]=0.0;           Ck[1][4]=0.0; 
  Ck[2][0]=1067./960;         Ck[2][1]=-29./480;          Ck[2][2]=3./640.;         Ck[2][3]=0.0;           Ck[2][4]=0.0; 
  Ck[3][0]=30251./26880;      Ck[3][1]=-7621./107520;     Ck[3][2]=159./17920;      Ck[3][3]=-5./7168;      Ck[3][4]=0.0; 
  Ck[4][0]=5851067./5160960;  Ck[4][1]=-100027./1290240;  Ck[4][2]=31471./2580480;  Ck[4][3]=-425./258048;  Ck[4][4]=35./294912; 
  // Ckbar
  Ckbar[0][0]=1.0;                 Ckbar[0][1]=0.0;                Ckbar[0][2]=0.0;                 Ckbar[0][3]=0.0;              Ckbar[0][4]=0.0; 
  Ckbar[1][0]=11./12;              Ckbar[1][1]=1./24;              Ckbar[1][2]=0.0;                 Ckbar[1][3]=0.0;              Ckbar[1][4]=0.0; 
  Ckbar[2][0]=863./960;            Ckbar[2][1]=77./1440;           Ckbar[2][2]=-17./5760;           Ckbar[2][3]=0.0;              Ckbar[2][4]=0.0;
  Ckbar[3][0]=215641./241920;      Ckbar[3][1]=6361./107520;       Ckbar[3][2]=-281./53760;         Ckbar[3][3]=367./967680;      Ckbar[3][4]=0.0; 
  Ckbar[4][0]=41208059./46448640;  Ckbar[4][1]=3629953./58060800;  Ckbar[4][2]=-801973./116121600;  Ckbar[4][3]=49879./58060800;  Ckbar[4][4]=-27859./464486400;  
  // dk
  dk[0][0]=1.0;           dk[0][1]=0.0;         dk[0][2]=0.0;         dk[0][3]=0.0;           dk[0][4]=0.0; 
  dk[1][0]=9./8;          dk[1][1]=-1./24;      dk[1][2]=0.0;         dk[1][3]=0.0;           dk[1][4]=0.0; 
  dk[2][0]=75./64;        dk[2][1]=-25./384;    dk[2][2]=3./640;      dk[2][3]=0.0;           dk[2][4]=0.0; 
  dk[3][0]=1225./1024;    dk[3][1]=-245./3072;  dk[3][2]=49./5120;    dk[3][3]=-5./7168;      dk[3][4]=0.0; 
  dk[4][0]=19845./16384;  dk[4][1]=-735./8192;  dk[4][2]=567./40960;  dk[4][3]=-405./229376;  dk[4][4]=35./294912; 
  // Qbar
  Qbar[0][0]=1./2;              Qbar[0][1]=0.0;             Qbar[0][2]=0.0;           Qbar[0][3]=0.0;              Qbar[0][4]=0.0; 
  Qbar[1][0]=13./24;            Qbar[1][1]=-1./24;          Qbar[1][2]=0.0;           Qbar[1][3]=0.0;              Qbar[1][4]=0.0; 
  Qbar[2][0]=401./720;          Qbar[2][1]=-31./480;        Qbar[2][2]=11./1440;      Qbar[2][3]=0.0;              Qbar[2][4]=0.0; 
  Qbar[3][0]=68323./120960;     Qbar[3][1]=-353./4480;      Qbar[3][2]=1879./120960;  Qbar[3][3]=-191./120960;     Qbar[3][4]=0.0; 
  Qbar[4][0]=2067169./3628800;  Qbar[4][1]=-40111./453600;  Qbar[4][2]=581./25920;    Qbar[4][3]=-28939./7257600;  Qbar[4][4]=2497./7257600; 
  // rk
  rk[0][0]=1./2;          rk[0][1]=0.0;           rk[0][2]=0.0;         rk[0][3]=0.0;          rk[0][4]=0.0; 
  rk[1][0]=9./16;         rk[1][1]=-1./16;        rk[1][2]=0.0;         rk[1][3]=0.0;          rk[1][4]=0.0; 
  rk[2][0]=75./128;       rk[2][1]=-25./256;      rk[2][2]=3./256.;     rk[2][3]=0.0;          rk[2][4]=0.0; 
  rk[3][0]=1225./2048;    rk[3][1]=-245./2048;    rk[3][2]=49./2048;    rk[3][3]=-5./2048;     rk[3][4]=0.0; 
  rk[4][0]=19845./32768;  rk[4][1]=-2205./16384;  rk[4][2]=567./16384;  rk[4][3]=-405./65536;  rk[4][4]=35./65536; 
}

int u0_multi(int tst, int ind, int ind_air, int ind_etain, double a, double b,double x, double* W){ 
  double eps=1e-8; 
  double k=2*M_PI; //////////////////////////
  double rho0=1.0; double p0=5./7;
  double Gamma;
  if (tst==0){ // Sod
    if(x<=a+(b-a)/2.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    else{
      W[0]=0.125; // rho
      W[1]=0.0; // u
      W[2]=0.1; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==2){ // Bizarrium
    if(x<=a+(b-a)/2.0){
      W[0]=1./0.7e-4; // rho
      W[1]=0.0;       // u
      W[2]=1e11;      // p
    }
    else{
      W[0]=1e4; // rho
      W[1]=250.0; // u
      W[2]=0.0; // p
    }
    return 0;
  }
  else if (tst==3){ // Onde acoustique
    W[0]=rho0 + eps*sin(k*x);  // rho
    W[1]=eps*sin(k*x);  // u
    W[2]=p0 + eps*sin(k*x);  // p

    return 0;
  }
  else if (tst==4){ // Sod symétrisé
    if(x<=a+(b-a)/4.0){
      W[0]=0.125; // rho
      W[1]=0.0;       // u
      W[2]=0.1;      // p
    }
    else if(a+3.0*(b-a)/4.0<=x){
      W[0]=0.125; // rho
      W[1]=0.0;  // u
      W[2]=0.1; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  else if (tst==1){ // LeBlanc
    Gamma=5./3;
    if(x<=a+(b-a)/3.0){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)/10; // p
    }
    else{
      W[0]=1e-3; // rho
      W[1]=0.0; // u
      W[2]=(Gamma-1)*1e-3*1e-9; // p
      //W[2]=1.0;
    }
    return 0;
  }
  else if (tst==5){ // Woodward 3 etats 
    if(x<=a+(b-a)/10){
      W[0]=1.0; // rho
      W[1]=0.0;       // u
      W[2]=1e3;      // p
    }
    else if(a+9.0*(b-a)/10<=x){
      W[0]=1.0; // rho
      W[1]=0.0;  // u
      W[2]=1e2; // p
    }
    else{  // if(2.0*(a+b)/5<x || x<3.0*(a+b)/5){
      W[0]=1.0; // rho
      W[1]=0.0; // u
      W[2]=1e-2; // p
    }
    return 0;
  }
  else if (tst==6){ // Shu_Oscher 
    if(x<=a+(b-a)/10){
      W[0]=27./7;        // rho
      W[1]=4*sqrt(35)/9;        // u
      W[2]=31./3; // p
    }
    else{  
      W[0]=1.0+sin(5.*x)/5; // rho
      W[1]=0; // u
      W[2]=1.0; // p
    }
    return 0;
  }
  double u;
  if(tst==100 ){ u=320.756; }
  else if(tst==101){ u=320.756; } // x=0    beta/gamma      POLE 1
  else if(tst==102){ u=331.046; } // x=0.2  beta/gamma
  else if(tst==103){ u=346.306; } // x=0.5  beta/gamma
  else if(tst==104){ u=361.359; } // x=0.8  beta/gamma
  else if(tst==105){ u=371.266; } // x=1    beta/gamma      POLE 2
  else if(tst==106){ u=474.570; } // phase gamma            POLE 3
  else if(tst==107){ u=1333.77; } // x=0    gamma/liquide   POLE 4
  else if(tst==108){ u=1383.13; } // x=0.2  gamma/liquide
  else if(tst==109){ u=1453.67; } // x=0.5  gamma/liquide
  else if(tst==110){ u=1520.7;  } // x=0.8  gamma/liquide
  else if(tst==111){ u=1563.72; } // x=1    gamma/liquide   POLE 5
  else if(tst==112){ u=400.266; }
  else if(tst==113){ u=1700; }
  else if(tst==114){ u=200; }

  if (tst==100){ // CHOC SIMPLE
    if(ind==ind_air){  W[0]=1.27; W[2]=1e5;} // rho et p
    else if(ind==ind_etain){  W[0]=7287.; W[2]=0;} // rho et p
    W[1]=-u;  //-1100; //-400.761; // u
    W[3]=300.0; // T
    W[4]=1.;  // xA
    W[5]=0.0; // xB
    W[6]=0.0; // xC
    return 0;
  } 
  else if(200>tst && tst>100) { // CHOC SYMETRISE
    //printf(" choc sym\n");
    if(x<=a+(b-a)/2.0){
      if(ind==ind_air){  W[0]=1.27; W[2]=1e5; } // rho et p
      else if(ind==ind_etain){ W[0]=7287.; W[2]=0; } // rho et p
      W[1]=+u;  //-1100; //-400.761; // u
      W[3]=300.0; // T
      W[4]=1.;  // xA
      W[5]=0.0; // xB
      W[6]=0.0; // xC
    }
    else{
      if(ind==ind_air){  W[0]=1.27; W[2]=1e5; } // rho et p
      else if(ind==ind_etain){  W[0]=7287.; W[2]=0; } // rho et p
      W[1]=-u;  //-1100; //-400.761; // u
      W[3]=300.0; // T
      W[4]=1.;  // xA
      W[5]=0.0; // xB
      W[6]=0.0; // xC
    }
    return 0;
  }    
  else if(tst==200) { // CHOC SYMETRISE    ETAIN DANS L'ETAT P4
    //printf(" choc sym\n");
    if(ind==ind_air){  
      W[0]=1.27;
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=1e5; // P
      W[3]=0; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    } // rho
    else if(ind==ind_etain){  
      W[0]=10247.863117; // rho
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=44.866639e9; // P
      W[3]=1965.674770; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    }
    return 0;
  }
  else if(tst==201) { // CHOC SYMETRISE    ETAIN DANS L'ETAT P3
  //printf(" choc sym\n");
    if(ind==ind_air){  
      W[0]=1.27;
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=1e5; // P
      W[3]=0; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    } // rho
    else if(ind==ind_etain){  
      W[0]=8518.284342; // rho
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=11.353850e9; // P
      W[3]=388.435466; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    }
    return 0;
  }
  else if(tst==202) { //Point hugoniot entre P3 et P4
  //printf(" choc sym\n");
    if(ind==ind_air){  
      W[0]=1.27;
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=1e5; // P
      W[3]=0; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    } // rho
    else if(ind==ind_etain){  
      W[0]=1./104.09e-6; // rho
      W[1]=0;  //-1100; //-400.761; // u
      W[2]=30e9; // P
      W[3]=1133.1; // T
      W[4]=0.;  // xA
      W[5]=1.0; // xB
      W[6]=0.0; // xC
    }
    return 0;
  }

  printf("Erreur dans le choix de tst (u0)\n");
  return 3;
}

void onde_acous(double x, double* W){
  double eps=1e-8; 
  double k=2*M_PI; //////////////////////////
  double rho0=1.0, p0=5./7;
  W[0]=rho0 + eps*sin(k*x);  // rho
  W[1]=eps*sin(k*x);         // u
  W[2]=p0 + eps*sin(k*x);    // p
}

int Shu_Oscher_droite(int n, int tst_air, double xg, double xd, double* RES){
  int res; 
  gsl_integration_glfixed_table * TTT = gsl_integration_glfixed_table_alloc(n);
  double xi,wi;
  double rho, eps, u=0;  
  RES[0]=0; RES[1]=0; RES[2]=0; RES[3]=0; RES[4]=0; RES[5]=0; RES[6]=0; RES[7]=0; RES[8]=0;
  for(int i=0; i<n; i++){
    res= gsl_integration_glfixed_point(xg, xd, i, &xi, &wi, TTT); if(res){return 1;}
    rho=1.0 + sin(5.0*xi)/5.0;  eps=epsilonEOS(tst_air, rho, 1.0 );
    RES[0]+= wi*rho;   // rho
    RES[1]+= wi/rho;   // rho
    RES[2]+= wi*u;     // u
    RES[3]+= wi*(eps + u*u/2);  // energie interne  eps
    RES[4]+= wi*eps;     // energie totale  e
    RES[5]+= wi*rho*u;   //rho.e
    RES[6]+= wi*rho*(eps+ u*u/2);   //rho.e
    RES[7]+= wi*rho*eps;   //rho.epsilon
    RES[8]+= wi*rho*u*u/2;   //rho.epsilon
  }
  for(int i=0; i<=8; i++){  RES[i]*=1./(xd-xg); }

  return 0;
}

int u0_bar_multi(int tst, int ind, int ind_air, int ind_etain, int tst_air, double a, double b, double xg, double xd, double* W){ 
  double Gamma;
  double xdis, xdis1, xdis2;
  double rhog, rhod, rhom, rho;
  double ug, ud, um, u;
  double pg, pd, pm, p;
  double eg, ed, em, e;
  double epsg, epsd, epsm, eps;
  double Vg, Vd;

  // integration par la bibliothèque  GNU Scientific Library
  int res;
  double* RES=malloc(9*sizeof(double));
  size_t n=10; 
  gsl_integration_glfixed_table * TTT = gsl_integration_glfixed_table_alloc(n);
  double xi,wi;

  if(tst==0 || tst==1 || tst==2){ // tube à choc à deux états
    if(tst==0){ // Sod
      xdis=a+(b-a)/2.0;
      rhog=1.0; rhod=0.125;
      ug=0.0; ud=0.0;
      pg=1.0; pd=0.1;
    }
    else if(tst==1){ // LeBlanc
      xdis=a+(b-a)/3.0;
      Gamma=5./3;
      rhog=1.0; rhod=1e-3; // rho
      ug=0.0; ud=0.0; // u
      pg=(Gamma-1)/10; pd=(Gamma-1)*1e-12; // p
    }
    else if(tst==2){  // Bizarrium
      xdis=a+(b-a)/2.0;
      rhog=1./7e-5; rhod=1e4; // rho
      ug=0.0;  ud=250.0; // u
      pg=1e11;  pd=0.0;   // p
    }
    else{printf("Erreur u0_bar_multi tst=%d\n",tst ); return 1;}

    epsg=epsilonEOS(tst_air, 1./rhog, pg); epsd=epsilonEOS(tst_air, 1./rhod, pd); // energie interne
    eg=epsg + ug*ug/2; ed=epsd + ud*ud/2;  // energie totale
      
    if(xd<=xdis){  // [xg < xd] < xdis
      W[0]=rhog; // rho
      W[1]=1./rhog; // tau
      W[2]=ug;   // u
      W[3]=eg;   // e
      W[4]=epsg;   // eps
      W[5]=rhog*ug;   // rho.u
      W[6]=rhog*eg;   // rho.e
      W[7]=rhog*epsg;   // rho.eps
      W[8]=rhog*ug*ug/2;   // rho.ekin
    }
    else if(xdis<=xg){  //  xdis < [xg < xd]
      W[0]=rhod; // rho
      W[1]=1./rhod; // tau
      W[2]=ud;   // u
      W[3]=ed;   // e
      W[4]=epsd;   // eps
      W[5]=rhod*ud;   // rho.u
      W[6]=rhod*ed;   // rho.e
      W[7]=rhod*epsd;   // rho.eps
      W[8]=rhod*ud*ud/2;   // rho.eps
    }
    else{  //  xg < xdis < xd
      Vg=rhog; Vd=rhod;
      W[0]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho
      Vg=1./rhog; Vd=1./rhod;
      W[1]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho
      Vg=ug; Vd=ud;
      W[2]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // u
      Vg=eg; Vd=ed;
      W[3]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // e
      Vg=epsg; Vd=epsd;
      W[4]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // eps
      Vg=rhog*ug; Vd=rhod*ud;
      W[5]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.u
      Vg=rhog*eg; Vd=rhod*ed;
      W[6]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.e
      Vg=rhog*epsg; Vd=rhod*epsd;
      W[7]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.eps
      Vg=rhog*ug*ug/2; Vd=rhod*ud*ud/2;
      W[8]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.ekin
    }
    return 0;
  }
  else if (tst==3){ // Onde acoustique
    W[0]=0; W[1]=0; W[2]=0; W[3]=0; W[4]=0; W[5]=0; W[6]=0; W[7]=0; W[8]=0; 
    for(int i=0; i<n; i++){
      res= gsl_integration_glfixed_point(xg, xd, i, &xi, &wi, TTT); if(res){return 1;}
      onde_acous(xi, RES);
      rho=RES[0]; u=RES[1]; p=RES[2]; 

      //printf("i=%d | x=%.15e, wi=%.15e, rho=%.15e, u=%.15e, p=%.15e\n",i,xi,wi,rho,u,p);  

      eps=epsilonEOS(tst_air, 1./rho, p);  // energie interne
      e=eps + u*u/2;   // energie totale
      
      W[0]+=wi*rho; // rho
      W[1]+=wi/rho; // tau
      W[2]+=wi*u;   // u
      W[3]+=wi*e;   // e
      W[4]+=wi*eps; // eps
      W[5]+=wi*rho*u; // rho.u
      W[6]+=wi*rho*e; // rho.e
      W[7]+=wi*rho*eps;  //rho.eps
      W[8]+=wi*rho*u*u/2;  //rho.ekin
    }
    //printf("   rho_bar=%.15e, tau_bar=%.15e\n",W[0],W[1]);  
    free(RES);

    for(int i=0; i<=8; i++){  W[i]*=1./(xd-xg); }

    return 0;
  }
  else if(tst==4 || tst==5){ // tube à choc à 3 états constants
    if (tst==4){ // Sod symétrisé
      xdis1=a+(b-a)/4.0;
      xdis2=a+3.0*(b-a)/4.0;
      rhog=0.125; rhod=0.125;   rhom=1.0; // rho
      ug=0.0;   ud=0.0;     um=0.0;    // u
      pg=0.1;   pd=0.1;     pm=1.0;  // p
    }
    else if (tst==5){ // Woodward 3 etats 
      xdis1=a+(b-a)/10;
      xdis2=a+9.0*(b-a)/10;
      rhog=1.0; rhod=1.0;  rhom=1.0;  // rho
      ug=0.0;   ud=0.0;    um=0.0;    // u
      pg=1e3;   pd=1e2;    pm=1e-2;    // p
    }

    epsg=epsilonEOS(tst_air, 1./rhog, pg); epsd=epsilonEOS(tst_air, 1./rhod, pd); epsm=epsilonEOS(tst_air, 1./rhom, pm); // energie interne
    eg=epsg + ug*ug/2; ed=epsd + ud*ud/2;  em=epsm + um*um/2; // energie totale

    if(xd<=xdis1){  // [xg < xd] < xdis1 < xdis2
      W[0]=rhog; // rho
      W[1]=1./rhog; // tau
      W[2]=ug;   // u
      W[3]=eg;   // e
      W[4]=epsg;   // eps
      W[5]=rhog*ug;   // rho.u
      W[6]=rhog*eg;   // rho.e
      W[7]=rhog*epsg;   // rho.eps
      W[8]=rhog*ug*ug/2;   // rho.ekin
    }
    else if(xd>=xdis2){  //  xdis1 < xdis2 < [xg < xd] 
      W[0]=rhod; // rho
      W[1]=1./rhod; // rho
      W[2]=ud;   // u
      W[3]=ed;   // e
      W[4]=epsd;   // eps
      W[5]=rhod*ud;   // rho.u
      W[6]=rhod*ed;   // rho.e
      W[7]=rhod*epsd;   // rho.eps
      W[8]=rhod*ud*ud/2;   // rho.ekin
    }
    else if (xg>=xdis1 && xd<=xdis2){  // xdis1 < [xg < xd] < xdis2 
      W[0]=rhom; // rho
      W[1]=1./rhom; // rho
      W[2]=um;   // u
      W[3]=em;   // e
      W[4]=epsm;   // eps
      W[5]=rhom*um;   // rho.u
      W[6]=rhom*em;   // rho.e
      W[7]=rhom*epsm;   // rho.eps
      W[8]=rhom*um*um/2;   // rho.eps
    }
    else if(xg<=xdis1 && xdis1<=xd){  // [xg < xdis1 < xd] < xdis2 
      Vg=rhog; Vd=rhom; xdis=xdis1;
      W[0]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho
      Vg=1./rhog; Vd=1./rhom;
      W[1]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // tau
      Vg=ug; Vd=um;
      W[2]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // u
      Vg=eg; Vd=em;
      W[3]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // e
      Vg=epsg; Vd=epsm;
      W[4]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // eps
      Vg=rhog*ug; Vd=rhom*um;
      W[5]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.u
      Vg=rhog*eg; Vd=rhom*em;
      W[6]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.e
      Vg=rhog*epsg; Vd=rhom*epsm;
      W[7]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.eps
      Vg=rhog*ug*ug/2; Vd=rhom*um*um/2;
      W[8]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.ekin
    }
    else if(xg<=xdis2 && xdis2<=xd){  //  xdis1 < [xg < xdis2 < xd]
      Vg=rhom; Vd=rhod; xdis=xdis2;
      W[0]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho
      Vg=1./rhom; Vd=1./rhod;
      W[1]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // tau
      Vg=um; Vd=ud;
      W[2]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // u
      Vg=em; Vd=ed;
      W[3]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // e
      Vg=epsm; Vd=epsd;
      W[4]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // eps
      Vg=rhom*um; Vd=rhod*ud;
      W[5]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.u
      Vg=rhom*em; Vd=rhod*ed;
      W[6]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.e
      Vg=rhom*epsm; Vd=rhod*epsd;
      W[7]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.eps
      Vg=rhom*um*um/2; Vd=rhod*ud*ud/2;
      W[8]=((xdis-xg)*Vg + (xd-xdis)*Vd)/(xd-xg); // rho.eps
    }
    return 0;
  }
  else if (tst==6){ // Shu_Oscher 
    xdis=a+(b-a)/10;
    rhog=27./7;        // rho
    ug=4*sqrt(35)/9;   // u
    pg=31./3;          // p    
    epsg=epsilonEOS(tst_air, 1./rhog, pg); // energie interne
    eg=epsg + ug*ug/2; // energie totale
    
      
    if(xd<=xdis){ //  [xg < xd] < xdis
      W[0]=rhog; // rho
      W[1]=1./rhog; // tau
      W[2]=ug;   // u
      W[3]=eg;   // e
      W[4]=epsg;   // eps
      W[5]=rhog*ug;   // rho.u
      W[6]=rhog*eg;   // rho.e
      W[7]=rhog*epsg;   // rho.eps
      W[8]=rhog*ug*ug/2;   // rho.ekin
    }
    else if(xdis<=xg){ // xdis < [xg < xd] 
      res=Shu_Oscher_droite(n, tst_air, xg, xd, W);
    }
    else{  //  [xg < xdis < xd] 
      res=Shu_Oscher_droite(n, tst_air, xdis, xd, W);

      Vg=rhog;
      W[0]+=(xdis-xg)*Vg/(xd-xg); // rho
      Vg=1./rhog;
      W[1]+=(xdis-xg)*Vg/(xd-xg); // tau
      Vg=ug;
      W[2]+=(xdis-xg)*Vg/(xd-xg); // u
      Vg=eg;
      W[3]+=(xdis-xg)*Vg/(xd-xg); // e
      Vg=epsg; 
      W[4]+=(xdis-xg)*Vg/(xd-xg); // eps
      Vg=rhog*ug; 
      W[5]+=(xdis-xg)*Vg/(xd-xg); // rho.u
      Vg=rhog*eg; 
      W[6]+=(xdis-xg)*Vg/(xd-xg); // rho.e
      Vg=rhog*epsg;
      W[7]+=(xdis-xg)*Vg/(xd-xg); // rho.eps
      Vg=rhog*ug*ug/2;
      W[8]+=(xdis-xg)*Vg/(xd-xg); // rho.ekin
    }
    return 0;
  }

  double uabs;
  if(tst==100 ){ uabs=320.756; }
  else if(tst==101){ uabs=320.756; } // x=0    beta/gamma      POLE 1
  else if(tst==102){ uabs=331.046; } // x=0.2  beta/gamma
  else if(tst==103){ uabs=346.306; } // x=0.5  beta/gamma
  else if(tst==104){ uabs=361.359; } // x=0.8  beta/gamma
  else if(tst==105){ uabs=371.266; } // x=1    beta/gamma      POLE 2
  else if(tst==106){ uabs=474.570; } // phase gamma            POLE 3
  else if(tst==107){ uabs=1333.77; } // x=0    gamma/liquide   POLE 4
  else if(tst==108){ uabs=1383.13; } // x=0.2  gamma/liquide
  else if(tst==109){ uabs=1453.67; } // x=0.5  gamma/liquide
  else if(tst==110){ uabs=1520.7;  } // x=0.8  gamma/liquide
  else if(tst==111){ uabs=1563.72; } // x=1    gamma/liquide   POLE 5
  else if(tst==112){ uabs=400.266; }
  else if(tst==113){ uabs=1700; }
  else if(tst==114){ uabs=200; }


  if (tst==100){ // CHOC SIMPLE
    if(ind==ind_air){  
      rho=1.27; p=1e5;
      eps=epsilonEOS(tst_air, 1./rho, p); // energie interne
      u=-uabs;
    }
    else if(ind==ind_etain){  
      rho=7287.; p=0;
      //eps=-9.371e-3; 
      eps= fE_VP(1./rho, p, 1);
      u=-uabs;
    }
  } 
  else if(100<tst && tst<200) { // CHOC SYMETRISE
    if(xd<=a+(b-a)/2.0){
      if(ind==ind_air){
        rho=1.27; p=1e5;
        eps=epsilonEOS(tst_air, 1./rho, p); // energie interne
        u=uabs;
      }
      else if(ind==ind_etain){
        rho=7287.; p=0;
        //eps=-9.371e-3;
        eps= fE_VP(1./rho, p, 1);
        u=uabs;
      }
    }
    else{
      if(ind==ind_air){
        rho=1.27; p=1e5;
        eps=epsilonEOS(tst_air, 1./rho, p); // energie interne
        u=-uabs;
      }
      else if(ind==ind_etain){
        rho=7287.; p=0;
        //eps=-9.371e-3;
        eps= fE_VP(1./rho, p, 1);        
        u=-uabs;
      }
    }
  }
  else if(tst==200) { // CHOC SYMETRISE    ETAIN DANS L'ETAT P4
    if(ind==ind_air){
      rho=1.27; p=1e5;
      eps=epsilonEOS(tst_air, 1./rho, p); // energie interne
      u=0;
    }
    else if(ind==ind_etain){
      rho=10247.863117; // rho
      u=0;   // u
      p=44.866639e9; // P
      eps=889466.87; // energie interne
    }
  }
  else if(tst==201) { // CHOC SYMETRISE    ETAIN DANS L'ETAT P3
    if(ind==ind_air){
      rho=1.27; p=1e5;
      eps=epsilonEOS(tst_air, 1./rho, p); // energie interne
      u=0;
    }
    else if(ind==ind_etain){
      rho=8518.284342; // rho
      u=0;     // u
      p=11.353850e9; // P
      eps=112608.341; // energie interne
    }
  }
  else{ printf("Erreur dans le choix de tst (u0_bar_multi)\n"); return 3;}

  W[0]=rho;    // rho
  W[1]=1./rho; // tau
  W[2]=u;     // u 
  W[3]=eps+u*u/2; // e
  W[4]=eps;       // eps
  W[5]=rho*u;     // rho.u
  W[6]=rho*(eps+u*u/2); // rho.e
  W[7]=rho*eps;   // rho.eps
  W[8]=rho*u*u/2; // rho.ekin
  return 0;
}

int initButcher(int sch, double** A, double* THETA, double* ALPHA){
  // RK1  (Matsuno forward-backward)
  if (sch==0){
    A[0][0]=1.0;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1.;
    return 0;
  }
  // RK2 (SSP) HEUN
  else if (sch==1){
    A[0][0]=1.0;

    THETA[0]=1.0/2.0; THETA[1]=1.0/2.0;

    ALPHA[0]=1.;
    return 0;
  }
  // RK3 (SSP)
  else if (sch==2){
    A[0][0]=1.0; 
    A[1][0]=1.0/4.0; A[1][1]=1.0/4.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/6.0; THETA[2]=2.0/3.0;

    ALPHA[0]=1.; ALPHA[1]=1./2;
    return 0;
  }
  // KUTTA ORDRE 4 (Classic RK4)
  else if (sch==3){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=1.0/2.0;
    A[2][0]=0.0;     A[2][1]=0.0;      A[2][2]=1.0;

    THETA[0]=1.0/6.0; THETA[1]=1.0/3.0; THETA[2]=1.0/3.0; THETA[3]=1.0/6.0;

    ALPHA[0]=1./2; ALPHA[1]=1./2; ALPHA[2]=1.;
    return 0;
  }
  // Cash-Karp (ordre 5)
  else if (sch==4){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;        A[1][1]=9./40;
    A[2][0]=3./10;        A[2][1]=-9./10;    A[2][2]=6./5;
    A[3][0]=-11./54;      A[3][1]=5./2;      A[3][2]=-70./27;     A[3][3]=-35./27;
    A[4][0]=1631./55296;  A[4][1]=175./512;  A[4][2]=575./13824;  A[4][3]=44275./110592;  A[4][4]=253./4096;

    THETA[0]=37./378; THETA[1]=0.; THETA[2]=250./621; THETA[3]=125./594; THETA[4]=0.; THETA[5]=512./1771;  

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=3./5; ALPHA[3]=1.; ALPHA[4]=7./8;
    return 0;
  }
  // Dormand-Prince
  else if (sch==5){
    A[0][0]=1.0/5.0; 
    A[1][0]=3./40;         A[1][1]=9./40;
    A[2][0]=44./45;        A[2][1]=-56./15;        A[2][2]=32./9;
    A[3][0]=19372./6561;   A[3][1]=-25360./2187;   A[3][2]=64448./6561;   A[3][3]=-212./729;
    A[4][0]=9017./3168;    A[4][1]=-355./33;       A[4][2]=46732./5247;  A[4][3]=49./176;    A[4][4]=-5103./18656;
    A[5][0]=35./384;       A[5][1]=0.;             A[5][2]=500./1113;    A[5][3]=125./192;   A[5][4]=-2187./6784;  A[5][5]=11./84;

    THETA[0]=35./384; THETA[1]=0.; THETA[2]=500./1113; THETA[3]=125./192; THETA[4]=-2187./6784; THETA[5]=11./84; THETA[6]=0.;

    ALPHA[0]=1./5; ALPHA[1]=3./10; ALPHA[2]=4./5; ALPHA[3]=8./9; ALPHA[4]=1.; ALPHA[5]=1.;
    return 0;
  }
  // HEUN  ordre 3
  else if (sch==6){
    A[0][0]=1.0/3; 
    A[1][0]=0; A[1][1]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=0; THETA[2]=3./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
    return 0;
  }
  // RALSTON  ordre 3
  else if (sch==7){
    A[0][0]=1.0/2; 
    A[1][0]=0; A[1][1]=3.0/4;

    THETA[0]=2.0/9; THETA[1]=1./3; THETA[2]=4./9;

    ALPHA[0]=1./2; ALPHA[1]=3./4;
    return 0;
  }
  // Bogacki-Shampine ordre 3 
  else if (sch==8){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=2./9; THETA[1]=1.0/3; THETA[2]=4./9; THETA[3]=0;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
    return 0;
  }
  // Bogacki-Shampine ordre 2 
  else if (sch==9){
    A[0][0]=1.0/2.0; 
    A[1][0]=0.0;     A[1][1]=3.0/4;
    A[2][0]=2./9;    A[2][1]=1./3;      A[2][2]=4./9;

    THETA[0]=7./24; THETA[1]=1.0/4; THETA[2]=1./3; THETA[3]=1./8;

    ALPHA[0]=1./2; ALPHA[1]=3./4; ALPHA[2]=1.; 
    return 0;
  }
  // KUTTA odre 3
  else if (sch==10){
    A[0][0]=1./2; 
    A[1][0]=-1.; A[1][1]=2.0;

    THETA[0]=1./6; THETA[1]=2./3; THETA[2]=1./6;

    ALPHA[0]=1./2; ALPHA[1]=1.;
    return 0;
  }
  // KUTTA ordre 2 (Explicit Midpoint Method) 
  else if (sch==11){
    A[0][0]=1.0/2;

    THETA[0]=0.0; THETA[1]=1.0;

    ALPHA[0]=1./2;
    return 0;
  }
  // KUTTA ordre 2 Ralston
  else if (sch==12){
    A[0][0]=2.0/3;

    THETA[0]=1.0/4; THETA[1]=3./4;

    ALPHA[0]=2./3;
    return 0;
  }
  // EULER forward
  else if (sch==13){
    THETA[0]=0; 
    return 0;
  }
  // SSPRK4(5,4) Spiteri-Ruuth
  else if (sch==14){
    A[0][0]=0.39175222700392; 
    A[1][0]=0.21766909633821;  A[1][1]=0.36841059262959;
    A[2][0]=0.08269208670950;  A[2][1]=0.13995850206999;  A[2][2]=0.25189177424738;
    A[3][0]=0.06796628370320;  A[3][1]=0.11503469844438;  A[3][2]=0.20703489864929;  A[3][3]=0.54497475021237;

    THETA[0]=0.14681187618661; THETA[1]=0.24848290924556; THETA[2]=0.10425883036650; THETA[3]=0.27443890091960; THETA[4]=0.22600748319395; 

    ALPHA[0]=0.39175222700392; ALPHA[1]=0.58607968896779; ALPHA[2]=0.47454236302687; ALPHA[3]=0.93501063100924; 
    return 0;
  }
  // SSPRK3(4,3) Spiteri-Ruuth
  else if (sch==15){
    A[0][0]=0.5; 
    A[1][0]=0.5;   A[1][1]=0.5;
    A[2][0]=1./6;  A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=0.5; 

    ALPHA[0]=0.5; ALPHA[1]=1.0; ALPHA[2]=0.5; 
    return 0;
  }
  // SSPRK3(5,3) Spiteri-Ruuth
  else if (sch==16){
    A[0][0]=0.37726891511710; 
    A[1][0]=0.75453783023419;  A[1][1]=0.37726891511710;
    A[2][0]=0.49056882269314;  A[2][1]=0.16352294089771;  A[2][2]=0.16352294089771;
    A[3][0]=0.78784303014311;  A[3][1]=0.14831273384724;  A[3][2]=0.14831273384724;  A[3][3]=0.34217696850008;

    THETA[0]=0.19707596384481; THETA[1]=0.11780316509765; THETA[2]=0.11709725193772; THETA[3]=0.27015874934251; THETA[4]=0.29786487010104; 

    ALPHA[0]=0.37726891511710; ALPHA[1]=0.75453783023419; ALPHA[2]=0.49056882269314; ALPHA[3]=0.78784303014311; 
    return 0;
  }
  // SPPRK1(3,1) Spiteri-Ruuth
  else if (sch==17){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./3; ALPHA[1]=2./3;
    return 0;
  }
  // SPPRK2(3,2) Spiteri-Ruuth
  else if (sch==18){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;

    THETA[0]=1./3; THETA[1]=1./3; THETA[2]=1./3;

    ALPHA[0]=1./2; ALPHA[1]=1.;
    return 0;
  }
  // SPPRK2(4,2) Spiteri-Ruuth
  else if (sch==19){
    A[0][0]=1./3; 
    A[1][0]=1./3; A[1][1]=1./3;
    A[2][0]=1./3; A[2][1]=1./3;  A[2][2]=1./3;

    THETA[0]=1./4; THETA[1]=1./4; THETA[2]=1./4; THETA[3]=1./4;

    ALPHA[0]=1./3; ALPHA[1]=2./3; ALPHA[2]=1.;
    return 0;
  }
  // SPPRK1(2,1) Spiteri-Ruuth
  else if (sch==20){
    A[0][0]=1./2; 

    THETA[0]=1./2; THETA[1]=1./2; 

    ALPHA[0]=1./2;
    return 0; 
  }
  // SPPRK3(4,3) Spiteri-Ruuth
  else if (sch==21){
    A[0][0]=1./2; 
    A[1][0]=1./2; A[1][1]=1./2;
    A[2][0]=1./6; A[2][1]=1./6;  A[2][2]=1./6;

    THETA[0]=1./6; THETA[1]=1./6; THETA[2]=1./6; THETA[3]=1./2;

    ALPHA[0]=1./2; ALPHA[1]=1.; ALPHA[2]=1./2;
    return 0;
  }
  else{ 
    printf("Erreur choix sch (initButcher)\n");
    return 4;
  }
}

int init_fluides(int tst_multi, int* indFLUIDE, int N, int ind_air, int ind_etain){
  int Ntemp=N/6;

  if(tst_multi==0){ // Gaz parfait air
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_air;
    }
  }
  else if(tst_multi==1){  // Bizaruim
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_air;
    }
  }
  else if(tst_multi==2){
    for (int i = 0; i<N; i++){
      indFLUIDE[i]=ind_etain;
    }
  }
  else if(tst_multi==3){
    for(int i=0; i<Ntemp; i++){
      indFLUIDE[i]=ind_air;
      indFLUIDE[N-1-i]=ind_air;
    }
    for(int i=Ntemp; i<N-Ntemp; i++){
      indFLUIDE[i]=ind_etain ;
    }
  }
  else if(tst_multi==4){
    for(int i=0; i<Ntemp; i++){
      indFLUIDE[N-1-i]=ind_air;
    }
    for(int i=0; i<N-Ntemp; i++){
      indFLUIDE[i]=ind_etain ;
    }
  }
  else{printf("Erreur choix tst_multi= %d\n",tst_multi ); return 1;}
  return 0;
}

// Schémas GoHy ..................................................................................................................................
// fonction pour passer d'une valeur point-wise à une valeur average sur la maille à l'odre 3
double fpw_to_av(double* Ypw){
  return (11./12)*Ypw[0] + (1./24)*(Ypw[1] + Ypw[-1]);
}

// fonction pour passer d'une valeur average à une valeur point-wise sur la maille à l'odre 3
double fav_to_pw(double* Yav){
  return (13./12)*Yav[0] - (1./24)*(Yav[1] + Yav[-1]);
}

// Reconstruction des positions des centres des mailles à l'odre 3
double f_grid(double* Xd){
  return (9./16)*(Xd[1]+Xd[0]) + (-1./16)*(Xd[2]+Xd[-1]);
}

double fdtu(double* RHO0, double* P, double dx){
  return (2./(RHO0[-1]+ RHO0[0]))*(P[-1]-P[0])/dx;
}

double fdtp(double* RHO0, double* RHOC2, double* U, double dx){
  return ((RHOC2[-1] + RHOC2[0])/(RHO0[-1] + RHO0[0]))*(U[-1]-U[0])/dx;
}

double fdttu(double* RHO, double* RHOC2, double* U, double dx){
  double dttu_1= (RHOC2[-1] + RHOC2[0])/(RHO[-1]*RHO[-1] + RHO[0]*RHO[0]);
  dttu_1*= ((U[1]+U[-2])-(U[-1]+U[0]))/(2*dx*dx);

  double dttu_2= 2./dx;
  dttu_2*= 1./(RHO[-1]+RHO[0]); 
  dttu_2*= RHOC2[0]/RHO[0] - RHOC2[-1]/RHO[-1];
  dttu_2*= (U[0] - U[-1])/dx;
 
  //printf(" dRC2=%g\n", (RHOC2[-1] + RHOC2[0]));
  //printf(" dttu_1=%g, dttu_2=%g ",dttu_1,dttu_2 );
  return dttu_1 + dttu_2;
}

double fdttp(double* RHO, double* RHOC2, double* G, double* P, double* U, double dx){
  double dttp_1= (RHOC2[-1]+RHOC2[0])/(RHO[-1]*RHO[-1]*RHO[-1] + RHO[0]*RHO[0]*RHO[0]);
  dttp_1*= (RHO[0]-RHO[-1])/dx;
  dttp_1*= (P[-1]-P[0])/dx;

  double dttp_2= 2*(RHOC2[-1]*G[-1] + RHOC2[0]*G[0])/(RHO[-1] + RHO[0]);
  dttp_2*= ((U[0]-U[-1])/dx)*((U[0]-U[-1])/dx);

  double dttp_3= (RHOC2[-1] + RHOC2[0])/(RHO[-1]*RHO[-1] + RHO[0]*RHO[0]);
  dttp_3*= ((P[1]+P[-2])-(P[-1]+P[0]))/(2*dx*dx);
  
  //printf(" dttp_1=%g, dttp_2=%g, dttp_3=%g ",dttp_1, dttp_2, dttp_3 );
  return dttp_1 + dttp_2 + dttp_3;
}

// Calcul la valeur à l'interface à partir des valeurs aux centres des mailles
double f_inter(int so, double* Yc, double** rk){
  // voir CRAS 2016 p214
  double res=0.0;
  int k=so-2;
  for(int l=0; l<=k; l++){
    res+=rk[k][l]*(Yc[l] + Yc[-l-1]);
  }
  return res;
}



// #################################################################################################################################
//             CALCUL ERREUR 
int err_palier_u(int ind_u, double xa, double xb, int nx, double* Xc, double* Xd, double* U, double U_ref, double* err){
  int conda=0, condb=0, ia, ib;
  double S;
  
  // 
  if(ind_u==0){ // U décentré
    for(int i=0; i<nx; i++){
      if(Xc[i]<xa){ia=i;} else{conda=1;}
      if(Xc[i]<xb){ib=i;} else{condb=1;}
      if(conda*condb){break;}
    }
    ia+=1;
    
    S=(Xc[ia]-xa)*U[ia];
    S+=(xb-Xc[ib])*U[ib+1];
    for(int i=ia; i<=ib-1; i++){
      S+=U[i+1]*(Xc[i+1]-Xc[i]);
    }
    S=S/(xb-xa);
  }
  else if(ind_u==1){ // U centré

    for(int i=0; i<nx+1; i++){
      if(Xd[i]<xa){ia=i;} else{conda=1;}
      if(Xd[i]<xb){ib=i;} else{condb=1;}
      if(conda*condb){break;}
    }
    ia+=1;
    
    S=(Xd[ia]-xa)*U[ia-1];
    S+=(xb-Xd[ib])*U[ib];
    for(int i=ia; i<=ib-1; i++){
      //printf(" i=%d\n",i );
      S+=U[i]*(Xd[i+1]-Xd[i]);
    }
    S=S/(xb-xa);
  }
  
  err[0]=fabs(U_ref-S);
  err[1]=fabs(U_ref-S)/fabs(U_ref);

  return 0;
}


int err_palier_c(int n, double xa, double xb, int nx, double* Xd, double** Sol, double* Etat_ref, double** Err){
  int conda=0, condb=0, ia=-1, ib=-1;
  double S;

  // 
  for(int i=0; i<nx+1; i++){
    if(Xd[i]<xa){ia=i;} else{conda=1;}
    if(Xd[i]<xb){ib=i;} else{condb=1;}
    if(conda*condb){break;}
  }
  ia+=1;
  
  for(int j=0; j<n; ++j){
    S=(Xd[ia]-xa)*Sol[j][ia-1];
    S+=(xb-Xd[ib])*Sol[j][ib];
    for(int i=ia; i<=ib-1; i++){
      S+=Sol[j][i]*(Xd[i+1]-Xd[i]);
    }
    S=S/(xb-xa);
    
    Err[j][0]=fabs(Etat_ref[j]-S);
    Err[j][1]=fabs(Etat_ref[j]-S)/fabs(Etat_ref[j]);
  }

  return 0;
}


int sol_err_exacte(int cas_test, int nx, int Rsch){
  int ind_u;

  if(Rsch==0||Rsch==4){ ind_u=1; }
  else if(Rsch==1||Rsch==2||Rsch==3){ ind_u=0; }
  else{printf(" erreur sch/100=%d\n",Rsch); return 1;}
 
  int n,res; double err;
  double temp, xaf,xbf,xa1,xb1;
  double Gamma,x,a,b,Tf,dx,eps;
  double epsilon=1e-12;

  // #variables - #u + #lambda = 13 - 1 = 12

  double** Sol=alloctabd(12, nx);  // palier
  double *Xc=malloc(nx*sizeof(double));
  double *Xd=malloc((nx+1)*sizeof(double));
  double *U;
  if(ind_u==1){
    U=malloc((nx)*sizeof(double));
  }
  else if(ind_u==0){
    U=malloc((nx+1)*sizeof(double));
  }
  
  double** ErrEP=alloctabd(13, 2); // palier
  double** ErrEF=alloctabd(13, 2); // etat final
  

  //#######################################################################################
  // lecture des données de la simulation
  
  // Ecriture dans un fichier des points (x,u(x)) de la solution sur grille centrée
  // si le schema est centrée i.e u est calculé au centres des mailles, on écrit u dans le fichier,
  // sinon on le remplace par des 0 (tout les autres types de schémas)  
  /*
       output_c.txt                        output_d.txt                    output_lambda.txt
  L1  : Xc  maillage centrée          L1  : Xc ou Xd  maillage          L1  : Xc  maillage centrée
  L2  : rho  volume spécifique        L2  : u vitesse                   L2  : fraction massique de la phase beta
  L3  : tau  volume spécifique                                          L3  : fraction massique de la phase gamma
  L4  : e   energie totale                                              L4  : fraction massique de la phase liquide
  L5  : eps   energie interne
  L6  : P   pression
  L7  : T   temperature
  L8  : C   celerite du son
  L9  : G   derivee fondamentale
  L10 : S   entropie
  */
  FILE *Res1,*Res2,*Res3,*Res4;
  if((Res1 = fopen("output_c.txt", "r")) == NULL){printf("Erreur ouverture output_c (sol_err_exacte)\n");  return 1;}
  if((Res2 = fopen("output_d.txt", "r")) == NULL){printf("Erreur ouverture output_d (sol_err_exacte)\n");  return 1;}
  if((Res3 = fopen("output_lambda.txt", "r")) == NULL){printf("Erreur ouverture output_lambda (sol_err_exacte)\n");  return 1;}
  if((Res4 = fopen("output_Xd.txt", "r")) == NULL){printf("Erreur ouverture output_lambda (sol_err_exacte)\n");  return 1;}
  
  // variables centrées aux mailles
  for(int i=0; i<nx; i++){
    fscanf(Res1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&Xc[i],&Sol[0][i],&Sol[1][i],&Sol[2][i],&Sol[3][i],&Sol[4][i],&Sol[5][i],&Sol[6][i],&Sol[7][i],&Sol[8][i]);
    fscanf(Res3,"%lf %lf %lf %lf\n",&temp,&Sol[9][i],&Sol[10][i],&Sol[11][i]);
  }
  // vitesse
  if(ind_u==1){
    for(int i=0; i<nx; i++){ 
      fscanf(Res2,"%lf %lf\n",&temp,&U[i]);
    }
  }
  else if(ind_u==0){
    for(int i=0; i<nx+1; i++){ 
      fscanf(Res2,"%lf %lf\n",&Xd[i],&U[i]);
    }
  }
  for(int i=0; i<nx+1; i++){ 
    fscanf(Res4,"%lf\n",&Xd[i]);
  }
  fclose(Res1);  fclose(Res2);  fclose(Res3);  fclose(Res4);
  

  //#######################################################################################
  // Calcul de l'erreur sur les solutions exactes
  //#######################################################################################
  if(cas_test==0||cas_test==1){ // cas test de Sod ou LeBlanc
    printf(" prise d'erreur pas implémenté\n"); return 1;
  }
  else if(cas_test==3){ //Onde acoustique
    printf(" prise d'erreur pas implémenté\n"); return 1;
  }
  else if(cas_test>=101 && cas_test<=112){ // Choc étain
    //printf(" Etain cas_test=%d (sol_err_exacte/func.c)\n",cas_test );
    a=-1e-3;   b=1e-3; 
    Tf=0.2e-6;
    // IDEE : lancer le calcul de la vrai solution grace au Newton de L'hugoniot
    int nb_choc, nhalf;
    double delta_choc;
    if(cas_test==101 || cas_test==106 || cas_test==107 || cas_test==108 || cas_test==109 || cas_test==110 || cas_test==111 || cas_test==112){
      nb_choc=1;
    }
    else{ nb_choc=2; }
    
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

    V0=Etat0[1];
    P0=Etat0[4];

    V1=Etat1[1];
    P1=Etat1[4];
    u1=Etat1[12];
    D1=Etat1[13];

    Vf=EtatFinal[1];
    Pf=EtatFinal[4];
    uf=EtatFinal[12];
    Df=EtatFinal[13];
  
    // simple choc
    if(cas_test==101 || cas_test==106 || cas_test==107 || cas_test==108 || cas_test==109 || cas_test==110 || cas_test==111 || cas_test==112){
      //printf("choc simple (sol_err_exacte/func.c)\n");
      choc1=V0*sqrt((Pf-P0)/(V0-Vf))*Tf - uf*Tf;

      // Etat final
      xaf=-choc1 + 0.2*(2*choc1);
      xaf=choc1 - 0.2*(2*choc1);

      res=err_palier_c(xaf, xbf, nx, 12, Xd, Sol, EtatFinal, ErrEF);
      // u
      res=err_palier_u(ind_u, xaf, xbf, nx, Xd, Xc, U, EtatFinal[12], ErrEF[12]);

    }
    else {
      //printf(" choc double (sol_err_exacte/func.c)\n");
      choc1=V0*sqrt((P1-P0)/(V0-V1))*Tf - uf*Tf;
      choc2=V1*sqrt((Pf-P1)/(V1-Vf))*Tf + u1*Tf - uf*Tf;
      //printf("choc1=%g choc2=%g\n",choc1,choc2 );

      // P1
      delta_choc=choc1-choc2;
      xa1=-choc1 + 0.2*delta_choc;
      xb1=-choc2 - 0.2*delta_choc;
      //res=err_palier(xa1, xb1, nx, Xd, P, P1, &errL1_1, &errL2_1);
      res=err_palier_c(12,xa1, xb1, nx, Xd, Sol, Etat1, ErrEP);
      // u
      res=err_palier_u(ind_u, xa1, xb1, nx, Xd, Xc, U, Etat1[12], ErrEP[12]);

      // Etat final
      xaf=-choc2 + 0.2*(2*choc2);
      xbf=choc2 - 0.2*(2*choc2);
      //res=err_palier(xaf, xbf, nx, Xd, P, Pf, &errL1_f, &errL2_f);
      res=err_palier_c(12,xaf, xbf, nx, Xd, Sol, EtatFinal, ErrEF);
      // u
      
      res=err_palier_u(ind_u, xaf, xbf, nx, Xd, Xc, U, EtatFinal[12], ErrEF[12]);
    }
  }
  else{
    //printf("Erreur cas test=%d (sol_reference.c)\n",cas_test );
    return 1;
  }

  // ##########################################################################
  //   Ecriture de l'erreur
  // ##########################################################################
  
  ////////////////////////////////////////////////////
  // BARPLOT commun a toutes les variables
  FILE *ResErrEP,*ResErrEF;
  if((ResErrEF = fopen("errEF.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  
  if((ResErrEP = fopen("errEP.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  if(cas_test>=101 && cas_test<=112){
    for(int j=0; j<=11; j++){
      fprintf(ResErrEF,"%.15e %.15e \n",ErrEF[j][0],ErrEF[j][1]);
      fprintf(ResErrEP,"%.15e %.15e \n",ErrEP[j][0],ErrEP[j][1]);
      
      //printf(" j=%d | %e %e (sol_err_exacte/func.c)\n",j,ErrEF[j][0],ErrEF[j][1]);
      //printf("        %e %e (sol_err_exacte/func.c)\n",ErrEP[j][0],ErrEP[j][1]);
    }
  }

  // Fichier texte pour toutes les variables sur un graphe different
  // le fichier texte contient une variable par colone et la colonne il y a 
  //    Etat Final Absolue
  //    Etat Final Relative
  //    Etat Premier Choc Absolue
  //    Etat Premier Choc Relative


  ////////////////////////////////////////////////////////////////////////////////
  // BARPLOT avec erreur relative et absolue juxtaposé 
  FILE *ResErr_ligne_abs, *ResErr_ligne_rel;
  if((ResErr_ligne_abs = fopen("err_ligne_abs.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  if((ResErr_ligne_rel = fopen("err_ligne_rel.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }

  if(cas_test>=101 && cas_test<=112){
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne_abs,"%.15e ",ErrEF[j][0]);
      fprintf(ResErr_ligne_rel,"%.15e ",ErrEF[j][1]);
    }
    fprintf(ResErr_ligne_abs,"\n");
    fprintf(ResErr_ligne_rel,"\n");
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne_abs,"%.15e ",ErrEP[j][0]);
      fprintf(ResErr_ligne_rel,"%.15e ",ErrEP[j][1]);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // PLOT classique 
  FILE *ResErr_ligne_EF_abs, *ResErr_ligne_EF_rel;
  FILE *ResErr_ligne_EP_abs, *ResErr_ligne_EP_rel;
  if((ResErr_ligne_EF_abs = fopen("err_ligne_EF_abs.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de err_ligne_EF_abs.txt (func.c/sol_err_exacte.c)\n");
    return 1;
  }
  if((ResErr_ligne_EF_rel = fopen("err_ligne_EF_rel.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de err_ligne_EF_rel.txt (func.c/sol_err_exacte.c)\n");
    return 1;
  }
  if((ResErr_ligne_EP_abs = fopen("err_ligne_EP_abs.txt", "w+")) == NULL){ 
    printf("Erreur d'ouverture de err_ligne_EP_abs.txt (func.c/sol_err_exacte.c)\n");
    return 1;
  }
  if((ResErr_ligne_EP_rel = fopen("err_ligne_EP_rel.txt", "w+")) == NULL){
    printf("Erreur d'ouverture de err_ligne_EP_rel.txt (func.c/sol_err_exacte.c)\n");
    return 1;
  }

  if(cas_test>=101 && cas_test<=112){
    fprintf(ResErr_ligne_EF_abs,"%d ",nx);
    fprintf(ResErr_ligne_EF_rel,"%d ",nx);
    fprintf(ResErr_ligne_EP_abs,"%d ",nx);
    fprintf(ResErr_ligne_EP_rel,"%d ",nx);
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne_EF_abs,"%.15e ",ErrEF[j][0]);
      fprintf(ResErr_ligne_EF_rel,"%.15e ",ErrEF[j][1]);
      fprintf(ResErr_ligne_EP_abs,"%.15e ",ErrEP[j][0]);
      fprintf(ResErr_ligne_EP_rel,"%.15e ",ErrEP[j][1]);
    }
    fprintf(ResErr_ligne_EF_abs,"\n");
    fprintf(ResErr_ligne_EF_rel,"\n");
    fprintf(ResErr_ligne_EP_abs,"\n");
    fprintf(ResErr_ligne_EP_rel,"\n");
  }
  
  
  return 0;
}



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

  //printf(" etat from func.c\n");


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
  Etat0[12]=u0;

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
