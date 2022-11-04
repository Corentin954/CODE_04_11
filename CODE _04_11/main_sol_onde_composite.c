#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "head_multi.h"



int main(int argc, char* argv[]){
  int err;
  double u_init;

  u_init=strtod(argv[1],NULL);
  int N=atoi(argv[2]);

  printf("  u_init=%g, N=%d\n",u_init,N );

  double* tabU=malloc(N*sizeof(double));
  double* tabV=malloc(N*sizeof(double));
  double* tabE=malloc(N*sizeof(double));
  double* tabP=malloc(N*sizeof(double));
  double* tabT=malloc(N*sizeof(double));
  double* tabS=malloc(N*sizeof(double));
  double* tabC=malloc(N*sizeof(double));
  double* tabDf=malloc(N*sizeof(double));
  double* tab_lambda_BETA=malloc(N*sizeof(double));
  double* tab_lambda_GAMMA=malloc(N*sizeof(double));
  double* tab_lambda_LIQ=malloc(N*sizeof(double));


  int onde_compression_isentropique(int N, double u_init, double* tabU, double* tabV, double* tabE, double* tabP, double* tabT, double* tabS, double* tabC, double* tabDf, double* tab_lambda_BETA, double* tab_lambda_GAMMA, double* tab_lambda_LIQ);

  double u1= 320.758508;

  
  //  ****************************************
  // Pour l'onde de compression isentropiques (OCI)
  FILE *fP, *fT, *fXA, *fXB, *fXC, *fC;
  if((fP = fopen("fichiers/Poci.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fT = fopen("fichiers/Toci.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fC = fopen("fichiers/Coci.txt", "w+")) == NULL){printf("erreur ouverture fichier C test\n");return 1;}
  if((fXA = fopen("fichiers/XAoci.txt", "w+")) == NULL){printf("erreur ouverture fichier XA test\n");return 1;}
  if((fXB = fopen("fichiers/XBoci.txt", "w+")) == NULL){printf("erreur ouverture fichier XB test\n");return 1;}
  if((fXC = fopen("fichiers/XCoci.txt", "w+")) == NULL){printf("erreur ouverture fichier XC test\n");return 1;}
  FILE *fV, *fE, *fS, *fU, *fDf;
  if((fV = fopen("fichiers/Voci.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fE = fopen("fichiers/Eoci.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  if((fS = fopen("fichiers/Soci.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  if((fU = fopen("fichiers/Uoci.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  if((fDf = fopen("fichiers/Dfoci.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  // *****************************************
  
  err= onde_compression_isentropique(N, u_init, tabU, tabV, tabE, tabP, tabT, tabS, tabC, tabDf, tab_lambda_BETA, tab_lambda_GAMMA, tab_lambda_LIQ);
  if(err){return err;}
  for(int i=0;  i<N ; i++){
    fprintf(fV, "%.15lf \n",tabV[i]);  fprintf(fE, "%.15lf \n",tabE[i] );
    fprintf(fP, "%.15lf \n",tabP[i]);  fprintf(fT, "%.15lf \n",tabT[i] );
    fprintf(fS, "%.15lf \n",tabS[i]);  fprintf(fC, "%.15lf \n",tabC[i] );
    fprintf(fDf, "%.15lf \n",tabDf[i]);  fprintf(fU, "%.15lf \n",tabU[i] );
    fprintf(fXA, "%.15lf \n",tab_lambda_BETA[i]);  fprintf(fXB, "%.15lf \n",tab_lambda_GAMMA[i]);  fprintf(fXC, "%.15lf \n",tab_lambda_LIQ[i]);
  }



}
