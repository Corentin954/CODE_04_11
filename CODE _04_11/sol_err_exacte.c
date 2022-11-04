#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "head_multi.h"
#include "EOS_etain.h"
#include "sol_head.h"

//  NON - UTILISe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int err_palier_c(int n, double xa, double xb, int nx, double* Xd, double** Sol, double* Etat_ref, double** Err);
int err_palier_u(int ind_u, double xa, double xb, int nx, double* Xc, double* Xd, double* U, double U_ref, double* err);


int sol_err_exacte(int cas_test, int nx, int ind_u){

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
  if((Res1 = fopen("output_c.txt", "w+")) == NULL){printf("Erreur ouverture output_c (sol_err_exacte)\n");  return 1;}
  if((Res2 = fopen("output_d.txt", "w+")) == NULL){printf("Erreur ouverture output_d (sol_err_exacte)\n");  return 1;}
  if((Res3 = fopen("output_lambda.txt", "w+")) == NULL){printf("Erreur ouverture output_lambda (sol_err_exacte)\n");  return 1;}
  if((Res4 = fopen("output_Xd.txt", "w+")) == NULL){printf("Erreur ouverture output_lambda (sol_err_exacte)\n");  return 1;}
  
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
  if(cas_test==0||cas_test==1){ // cas test de Sod ou LeBlanc
  	printf(" prise d'erreur pas implémenté\n"); return 1;
  }
  else if(cas_test==3){ //Onde acoustique
  	printf(" prise d'erreur pas implémenté\n"); return 1;
  }
  else if(cas_test>=101 && cas_test<=112){ // Choc étain
  	printf(" Etain cas_test=%d\n",cas_test );
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

	  uf=EtatFinal[12];
	  Df=EtatFinal[13];
	
		// simple choc
		if(cas_test==101 || cas_test==106 || cas_test==107 || cas_test==108 || cas_test==109 || cas_test==110 || cas_test==111 || cas_test==112){
      
		  choc1=V0*sqrt((Pf-P0)/(V0-Vf))*Tf - uf*Tf;

		  // Etat final
      xaf=-choc1 + 0.2*(2*choc1);
      xaf=choc1 - 0.2*(2*choc1);

      res=err_palier_c(xaf, xbf, nx, 12, Xd, Sol, EtatFinal, ErrEF);
      // u
      res=err_palier_u(ind_u, xaf, xbf, nx, Xd, Xc, U, EtatFinal[12], ErrEF[12]);

    }
		else {
		  choc1=V0*sqrt((P1-P0)/(V0-V1))*Tf - uf*Tf;
		  choc2=V1*sqrt((Pf-P1)/(V1-Vf))*Tf + u1*Tf - uf*Tf;
		  //printf("choc1=%g choc2=%g\n",choc1,choc2 );

      // P1
		  delta_choc=choc1-choc2;
		  xa1=-choc1 + 0.2*delta_choc;
		  xb1=-choc2 - 0.2*delta_choc;
      //res=err_palier(xa1, xb1, nx, Xd, P, P1, &errL1_1, &errL2_1);
      res=err_palier_c(xa1, xb1, 12, nx, Xd, Sol, Etat1, ErrEP);
      // u
      res=err_palier_u(ind_u, xa1, xb1, nx, Xd, Xc, U, Etat1[12], ErrEP[12]);

      // Etat final
      xaf=-choc2 + 0.2*(2*choc2);
      xbf=choc2 - 0.2*(2*choc2);
      //res=err_palier(xaf, xbf, nx, Xd, P, Pf, &errL1_f, &errL2_f);
      res=err_palier_c(xaf, xbf, 12, nx, Xd, Sol, EtatFinal, ErrEF);
      // u
      res=err_palier_u(ind_u, xaf, xbf, nx, Xd, Xc, U, EtatFinal[12], ErrEF[12]);

		}
  }
  else{
  	//printf("Erreur cas test=%d (sol_reference.c)\n",cas_test );
  	return 1;
  }
  
  // ##########################################################################
  // Ecriture de l'erreur

  // Fichier texte pour toutes les variables sur le meme graphe
  FILE *ResErrEP,*ResErrEF;
  if((ResErrEF = fopen("fichiers/GNUPLOT/ErrEF.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  
  if((ResErrEP = fopen("fichiers/GNUPLOT/ErrEP.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }

  if(cas_test>=101 && cas_test<=112){
	  for(int j=0; j<=12; j++){
      fprintf(ResErrEF,"%.15e %.15e \n",ErrEF[j][0],ErrEF[j][1]);
      fprintf(ResErrEP,"%.15e %.15e \n",ErrEP[j][0],ErrEP[j][1]);
    }
  }


  // Fichier texte pour toutes les variables sur un graphe different
  // le fichier texte contient une variable par colone et la colonne il y a 
  //    Etat Premier Choc Absolue
  // 		Etat Premier Choc Relative
  // 		Etat Final Absolue
  // 		Etat Final Relative

  FILE *ResErr_ligne;
  if((ResErr = fopen("fichiers/GNUPLOT/err_ligne.txt", "w+")) == NULL){ 
  	printf("Erreur d'ouverture de sol_reference.txt (sol_reference/sol_err_exacte.c)\n");
    return 1;
  }
  

  if(cas_test>=101 && cas_test<=112){
  	for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne,"%.15e ",ErrEF[j][0]);
    }
    fprintf(ResErr_ligne,"\n");
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne,"%.15e ",ErrEF[j][1]);
    }
    fprintf(ResErr_ligne,"\n");
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne,"%.15e ",ErrEP[j][0]);
    }
    fprintf(ResErr_ligne,"\n");
    for(int j=0; j<=12; j++){
      fprintf(ResErr_ligne,"%.15e ",ErrEP[j][1]);
    }
    fprintf(ResErr_ligne,"\n");
  }
  
	return 0;
}


int err_palier_u(int ind_u, double xa, double xb, int nx, double* Xc, double* Xd, double* U, double U_ref, double* err){
  int conda=0, condb=0, ia, ib;
  double S;

  // 
	if(ind_u==0){ // U décentré
	  for(int i=0; i<nx; i++){
	  	if(xa<Xc[i]){ia=i;} else{conda=1;}
	  	if(xb<Xc[i]){ib=i;} else{condb=1;}
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
	  	if(xa<Xd[i]){ia=i;} else{conda=1;}
	  	if(xb<Xd[i]){ib=i;} else{condb=1;}
	  	if(conda*condb){break;}
	  }
	  ia+=1;
	  
	  S=(Xd[ia]-xa)*U[ia-1];
	  S+=(xb-Xd[ib])*U[ib];
	  for(int i=ia; i<=ib-1; i++){
	    S+=U[i]*(Xd[i+1]-Xd[i]);
	  }
	  S=S/(xb-xa);
	}
  
  err[0]=fabs(U_ref-S);
  err[1]=fabs(U_ref-S)/fabs(U_ref);

  return 0;
}


int err_palier_c(int n, double xa, double xb, int nx, double* Xd, double** Sol, double* Etat_ref, double** Err){
  int conda=0, condb=0, ia, ib;
  double S;

  // 
  for(int i=0; i<nx+1; i++){
  	if(xa<Xd[i]){ia=i;} else{conda=1;}
  	if(xb<Xd[i]){ib=i;} else{condb=1;}
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

