#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrice.h"
#include "EOS_etain.h"


// *******   Newton pour calculer le point triple  ********
/* Valeur d'entree :
    Va,Ea, Vb,Eb, Vc,Ec

    corespodances :  1-a   2-b   3-c

  Valeur sortie
   VA,EA, VB,EB, VC,EC
*/
void Newton_triple(double* A, double* InvA, double* Delta, double* dX, double Va, double Ea, double Vb, double Eb, double Vc, double Ec, double* VA, double* EA, double* VB, double* EB, double* VC, double* EC){
  int dim=6;
  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  double Pc, Tc, Sc, Gc, dPvc, dPec, dTvc, dTec;
  // phase 1
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(Va,1); dPea=dPdE(1);
  dTva=dTdV(Va,1); dTea=dTdE(1);
  //phase 2
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(Vb,2); dPeb=dPdE(2);
  dTvb=dTdV(Vb,2); dTeb=dTdE(2);
  // phase 3
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;
  dPvc=dPdV(Vc,3); dPec=dPdE(3);
  dTvc=dTdV(Vc,3); dTec=dTdE(3);

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
  printf("  * InvA = *\n");
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
   
  double correc=1.;

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
void triple(int Nmax, double epsilon, double* Va0, double* Ea0, double* Vb0, double* Eb0, double* Vc0, double* Ec0){
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


  // phase 1
  Ea=*Ea0; Va=*Va0;
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase 2
  Eb=*Eb0; Vb=*Vb0;
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase 3
  Ec=*Ec0; Vc=*Vc0;
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;


  double critere=1.;
  /*
  printf("init" );
  printf("  P_beta=%.15g P_gamma=%.15g P_liq=%.15g (GPa) ",fP(Va,Ea, 1)*1e-9,fP(Vb,Eb, 2)*1e-9,fP(Vc,Ec, 3)*1e-9 );
  printf("| T_beta=%.15g  T_gamma=%.15g T_liq=%.15g (K)\n",fT(Va,Ea, 1),fT(Vb,Eb, 2),fT(Vc,Ec, 3) );
  */

  while(critere>=epsilon && n<Nmax){
    Newton_triple(A, InvA, Delta, dX, Va, Ea, Vb, Eb, Vc, Ec, &VA, &EA, &VB, &EB, &VC, &EC);

    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    Vc=VC; Ec=EC;
    /*
    printf("n=%d",n );
    printf("  P_beta=%.15g P_gamma=%.15g P_liq=%.15g (GPa) ",fP(Va,Ea, 1)*1e-9,fP(Vb,Eb, 2)*1e-9,fP(Vc,Ec, 3)*1e-9 );
    printf("| T_beta=%.15g  T_gamma=%.15g T_liq=%.15g (K)",fT(Va,Ea, 1),fT(Vb,Eb, 2),fT(Vc,Ec, 3) );
    */
    //critere=fabs(Ga-Gb)+fabs(Ga-Gc);
    critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb) + fabs(dX[4]/Vc) + fabs(dX[5]/Ec);
    //printf(" | critere= %.3g\n",critere );

    n++;
  }

  *Va0=Va; *Ea0=Ea; 
  *Vb0=Vb; *Eb0=Eb; 
  *Vc0=Vc; *Ec0=Ec;

  free(A); free(InvA);
  free(Delta);
  free(dX);
}


/*  ****  Calcul du Point TRIPLE  ****
   utilisé dans  main_calfront.c
*/
void point_triple(double* Pt, double* Tt, double* Vat, double* Eat, double* Vbt, double* Ebt, double* Vct, double* Ect){
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

  triple(Nmax, epsilon, &Va0, &Ea0, &Vb0, &Eb0, &Vc0, &Ec0);
  

  // phase 1
  Ea=Ea0; Va=Va0;
  Pa=fP(Va,Ea, 1); Ta=fT(Va,Ea, 1);
  Sa=fS(Va,Ea, 1); Ga=Ea+Pa*Va-Ta*Sa;
  //phase 2
  Eb=Eb0; Vb=Vb0;
  Pb=fP(Vb,Eb, 2); Tb=fT(Vb,Eb, 2);
  Sb=fS(Vb,Eb, 2); Gb=Eb+Pb*Vb-Tb*Sb;
  // phase 3
  Ec=Ec0; Vc=Vc0;
  Pc=fP(Vc,Ec, 3); Tc=fT(Vc,Ec, 3);
  Sc=fS(Vc,Ec, 3); Gc=Ec+Pc*Vc-Tc*Sc;

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

// *******  Newton sur les lignes de chgmt de phase  ********

/* 
  Utile dans la fonction  ligne_double

  Valeur d'entree :
   Va,Ea, Vb,Eb, Vc,Ec

   corespodances :  1-a   2-b   3-c

  Valeur sortie
   VA,EA, VB,EB, VC,EC

*/
void Newton_double(double* A, double* InvA, double* Delta, double* dX, double Pcst, int phaseA, int phaseB, double Va, double Ea, double Vb, double Eb, double* VA, double* EA, double* VB, double* EB){
  int dim=4;

  // Grandeurs thermodynamiques
  double Pa, Ta, Sa, Ga, dPva, dPea, dTva, dTea;
  double Pb, Tb, Sb, Gb, dPvb, dPeb, dTvb, dTeb;
  // phase A
  Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
  Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
  dPva=dPdV(Va,phaseA); dPea=dPdE(phaseA);
  dTva=dTdV(Va,phaseA); dTea=dTdE(phaseA);
  // phase B
  Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
  Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;
  dPvb=dPdV(Vb,phaseB); dPeb=dPdE(phaseB);
  dTvb=dTdV(Vb,phaseB); dTeb=dTdE(phaseB);


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
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %.15lf\n", det_a);}
  
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
void algo_double(int Nmax, double epsilon, int phaseA, int phaseB, double P, double Va0, double Ea0, double Vb0, double Eb0, double* pVa, double* pEa, double* pVb, double* pEb){
  int n=1;
  int dim=4;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta = (double *) malloc(sizeof(double) * (dim));
  double* dX = (double *) malloc(sizeof(double) * (dim));

  double VA, EA, VB, EB;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;

  Ea=Ea0; Va=Va0;
  Eb=Eb0; Vb=Vb0;


  double critere=1.;
  n=1;
  while(critere>=epsilon && n<Nmax){
    Newton_double(A, InvA, Delta, dX, P, phaseA, phaseB, Va, Ea, Vb, Eb, &VA, &EA, &VB, &EB);
    Va=VA; Ea=EA;
    Vb=VB; Eb=EB;
    critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb);

    Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
    Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
    Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
    Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;

    //printf(" n=%d T=%.15g  critere= %.3g ",n,fT(Va,Ea,phaseA),critere );
    //printf(" Ga-Gb=%.3g  Pa-Pb=%.3g  Ta-Tb=%.3g\n",Ga-Gb,Pa-Pb,Ta-Tb );
    n++;
  }
  *pVa=Va;
  *pEa=Ea;
  *pVb=Vb;
  *pEb=Eb;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
}

//  *******  Fonction qui calcule les points (Va,Ea) et (Vb,Ec) d'un point d'une ligne de chgmt de phases
// On rentre une pression et les deux phases que l'on souhaite
//   utilisé dans  main_calfront.c

void point_double(int phaseA, int phaseB, double P, double* pT, double* pVa, double* pEa, double* pVb, double* pEb){
  int Nmax=100;
  double epsilon=1e-10;

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


  algo_double(Nmax, epsilon, phaseA, phaseB, P, Va0, Ea0, Vb0, Eb0, pVa, pEa, pVb, pEb);

  
  // phase a
  Ea=*pEa; Va=*pVa;
  Pa=fP(Va,Ea, phaseA); Ta=fT(Va,Ea, phaseA);
  Sa=fS(Va,Ea, phaseA); Ga=Ea+Pa*Va-Ta*Sa;
  // phase b
  Eb=*pEb; Vb=*pVb;
  Pb=fP(Vb,Eb, phaseB); Tb=fT(Vb,Eb, phaseB);
  Sb=fS(Vb,Eb, phaseB); Gb=Eb+Pb*Vb-Tb*Sb;

  //sortie des donnes du point triple
  *pVa=Va; *pEa=Ea; 
  *pVb=Vb; *pEb=Eb; 
  *pT=(Ta+Tb)/2.0;


  
  printf("\nResultats double P=%g:\n",P);
  printf(" Ga= %g, Gb= %g,  Ga-Gb= %g\n",Ga, Gb, Ga-Gb );
  printf(" Pa= %g, Pb= %g,  Pa-Pb= %g  (GPa)\n",Pa*1e-9, Pb*1e-9, (Pa-Pb)*1e-9 );
  printf(" Ta= %g, Tb= %g,  Ta-Tb= %g (K)\n",Ta, Tb, Ta-Tb );

  printf("Va=%.5g Ea=%.5g\n",Va,Ea);
  printf("Vb=%.5g Eb=%.5g\n",Vb,Eb);
}

// Fonction qui renvoie le nombre de points sur les frontières en focntions d'un seul entier
void fnb_pts(int nb_points, int* pNTG1, int* pNPB1, int* pNPB2, int* pNTD2 , int* pNTD3, int* pNPH3, int* pNTG3, 
            int* pN12 ,int* pN13, int* pN23, int* pNA, int* pNB, int* pNC, int* pNAB, int* pNAC, int* pNBC ){
  

  int nb_points_front=ceil(nb_points); // /3    La dicretisation des frontières a moins d'importance NOP
  // A
  *pNTG1=nb_points_front;
  *pNPB1=2*nb_points;
  // B
  *pNPB2=12*nb_points_front;
  *pNTD2=3*nb_points_front;
  // C
  *pNTD3=2*nb_points_front;
  *pNPH3=15*nb_points_front;
  *pNTG3=8*nb_points_front;
  
  *pN12=5*nb_points; *pN13=3*nb_points; *pN23=40*nb_points;

  *pNA=*pN13+(*pNTG1-1)+(*pNPB1-1)+(*pN12-2);
  *pNB=*pN12+(*pNPB2-1)+(*pNTD2-1)+(*pN23-2);
  *pNC=*pN23+(*pNTD3-1)+(*pNPH3-1)+(*pNTG3-1)+(*pN13-2);
  *pNAB=2*(*pN12);
  *pNAC=2*(*pN13);
  *pNBC=2*(*pN23);
}



//  *******  Fonction qui calcule les points (V,E) des lignes de chgmt de phases
/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - [Pdeb, Pfin] : segment de presion sur lequel calcule la zone mixte
  - Np : nomnbre de points de discretisatoin du segment
  - Va0, Ea0, Vb0, Eb0 : Valeur initiale de (Va,Ea) et (Vb,Eb)
   (typiquement les coord du point triple pour chaque phase)

Arguments de sortie :
  - (segP, segT) : courbe de T(P) du diagramme de phase (frontières)
  - (segVa, segEa) : courbe pour la phase A
  - (segVb, segEb) : courbe pour la phase B

*/
void ligne_double(int Nmax, double epsilon, int phaseA, int phaseB, double Pdeb, double Pfin, int Np, double Va0, double Ea0, double Vb0, double Eb0, double* segP, double* segT, double* segVa, double* segEa, double* segVb, double* segEb){
  int n=1;
  int dim=4;
  double Pcst;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta = (double *) malloc(sizeof(double) * (dim));
  double* dX = (double *) malloc(sizeof(double) * (dim));

  double hp=(Pfin-Pdeb)/(Np-1);

  double VA, EA, VB, EB;

  double Pa, Ta, Sa, Ga, Ea, Va;
  double Pb, Tb, Sb, Gb, Eb, Vb;

  Ea=Ea0; Va=Va0;
  Eb=Eb0; Vb=Vb0;


  double critere;
  
  for(int i=0; i<Np; i++){
    //printf("\n****  i= %d ******\n",i );
    segP[i]=Pdeb+i*hp;
    Pcst=segP[i];
    //printf("  Pcst= %.15lf, hp= %.15lf, Pdeb= %.15lf\n", Pcst, hp, Pdeb );
    critere=1.;
    n=1;
    while(critere>=epsilon && n<Nmax){
      Newton_double(A, InvA, Delta, dX, Pcst, phaseA, phaseB, Va, Ea, Vb, Eb, &VA, &EA, &VB, &EB);
      Va=VA; Ea=EA;
      Vb=VB; Eb=EB;
      critere=fabs(dX[0]/Va) + fabs(dX[1]/Ea) + fabs(dX[2]/Vb) + fabs(dX[3]/Eb);
      //printf("    critere= %g\n",critere );
      n++;
    }
    segT[i]=fT(Va,Ea,phaseA);
    segVa[i]=Va;
    segEa[i]=Ea;
    segVb[i]=Vb;
    segEb[i]=Eb;
  }

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
}



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
int fVE(int phase, double P, double T, double* pV, double* pE){
  int Nmax=1e2;
  double epsilon=1e-15;
  int n=0;
  double f, fprime;
  double V,E, dV;
  
  double K0, N0, gamma0, Crv, theta0, T0, P0, rho0, v0, E0, Sr, ur;
  coeff(phase, &K0, &N0, &gamma0, &Crv, &theta0, &T0, &P0, &rho0, &v0, &E0, &Sr, &ur);

  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT); 


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
  else{printf("erreur phase= %d (fVE)\n",phase );}

  
  double critere=1.;
  
  // Newton sur E
  while(critere>epsilon && n<Nmax){
    f=fPs(V, phase) + gamma0*rho0*Crv*(T-ur*THETA(V, phase))-P;
    fprime=Ps_prime(V,phase) - Crv*ur*gamma0*gamma0*rho0*rho0*exp(gamma0*fx(V,phase));
    
    dV=f/fprime;
    V-=dV;
    critere=fabs(dV/V);
    //printf("  n= %d, V= %.15lf, E= %.15lf  critere=%g\n",n,V,fEs(V,phase) + (P-fPs(V,phase))/(rho0*gamma0),critere );
    n++;
  }
  E = fEs(V,phase) + (P-fPs(V,phase))/(rho0*gamma0);

  *pV=V;
  *pE=E;

  return 0;
}



// fonction pour determiner si un couple (V,E) est dans le triangle du poit triple
/*  -int *pind : pointeur sur l'indice (valeur retour)  
*/
int ftri(double V, double E){
  // Point triple
  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT); 



  // TRIANGLE
  double penteBA, penteBC, penteAC;

  penteBA=(EaT-EbT)/(VaT-VbT);
  penteBC=(EcT-EbT)/(VcT-VbT);
  penteAC=(EcT-EaT)/(VcT-VaT);

  double EBA=EbT + penteBA*(V-VbT);
  double EAC=EaT + penteAC*(V-VaT);
  double EBC=EbT + penteBC*(V-VbT);

  // triangle
  if(E>EBA && E<EBC && E>EAC){
    return 1;
  } 
  else{ return 0;}
}



// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N12.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsA(int N13, int N12, int NTG, int NPB, double PrefG, double PrefD, double TrefH, double TrefB, double* segVa13, double* segEa13, double* segVa12, double* segEa12, double** S){
  
  // le remplisage de S doit se faire à l'étape avant

  int err=0;
  int phase=1;
  int N=N13+(NTG-1)+(NPB-1)+(N12-2);

  double P,T;
  double V,E;

  /*
  Trefbas=T12[N12-1];
  TrefG=T13[N13-1];
  PrefG=P13[N13-1];
  PrefD=P12[N12-1];
  */

  double hTG=(TrefB-TrefH)/(NTG-1);
  double hPB=(PrefD-PrefG)/(NPB-1);

  
  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a droite 
  for(int j=0; j<=N13-1; j++){ // N13 points
    /*
    T=T13[j];
    P=P13[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    */
    V=segVa13[j];
    E=segEa13[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N13;

  // Parcours de la frontière droite
  P=PrefG;
  for(int j=1; j<=NTG-1; j++){ // N13 points
    T=TrefH+j*hTG;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTG-1;
  
  // Parcours de la frontière haute
  T=TrefB;
  for(int j=1; j<=NPB-2; j++){ // N13 points
    P=PrefG+j*hPB;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NPB-2;

  // Parcours de C12 de T=Trefbas jusqu'au point triple 
  for(int j=0; j<=N12-2; j++){ // N13 points
    /*
    T=T12[N12-1-j];
    P=P12[N12-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsA)\n");return err;}
    */
    V=segVa12[N12-1-j];
    E=segEa12[N12-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];

  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N23.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsB(int N12, int N23, int NPB, int NTD, double PrefG, double PrefD, double TrefB, double TrefH, double* segVb12, double* segEb12, double* segVb23, double* segEb23, double** S){
  
  // le remplisage de S doit se faire à l'étape avant


  int N=N12+(NPB-1)+(NTD-2)+(N23-1);
  int err=0;
  int phase=2;

  double P,T;
  double V,E;
  
  /*
  TrefB=T12[N12-1];
  TrefH=T23[N23-1];
  PrefG=P12[N12-1];
  PrefD=P23[N23-1];
  */

  double hTD=(TrefH-TrefB)/(NTD-1);
  double hPB=(PrefD-PrefG)/(NPB-1);

  
  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C12 duu point triple jusqu'a droite 
  for(int j=0; j<=N12-1; j++){ // N12 points
    /*
    T=T12[j];
    P=P12[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    */
    V=segVb12[j];
    E=segEb12[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N12;

  // Parcours de la frontière basse
  T=TrefB;
  for(int j=1; j<=NPB-1; j++){ 
    P=PrefG+j*hPB;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NPB-1;
  
  // Parcours de la frontière droite
  P=PrefD;
  for(int j=1; j<=NTD-2; j++){ 
    T=TrefB+j*hTD;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTD-2;


  //  Parcours de C23 de la droite jusqu'au point triple 
  for(int j=0; j<=N23-2; j++){ 
    /*
    T=T23[N23-1-j];
    P=P23[N23-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsB)\n");return err;}
    */
    V=segVb23[N23-1-j];
    E=segEb23[N23-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];

  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13, N23.
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple
  - double* P13, T13 : tableaux de N13 éléments qui contients les points de la courbe C13
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers C 
*/
int sommetsC(int N13, int N23, int NTD, int NPH, int NTG, double PrefGB, double PrefGH, double PrefDH, double PrefD, double TrefG, double TrefD, double Trefhaut, double* segVc13, double* segEc13, double* segVc23, double* segEc23 ,double** S){
  
  int N=N23+(NTD-1)+(NPH-1)+(NTG-2)+(N13-1) ;

  int err=0;
  int phase=3;

  double P,T;
  double V,E;
  
  /*
  Trefhaut=2500;
  TrefD=T23[N23-1];
  TrefG=T13[N13-1];
  PrefG=P13[N13-1];
  PrefD=P23[N23-1];
  */

  double hTD=(Trefhaut-TrefD)/(NTD-1);
  double hPD=(PrefDH-PrefD)/(NTD-1);
  
  double hPH=(PrefGH-PrefDH)/(NPH-1);
  // la frontière gauche n'est pas horizontale
  double hTG=(TrefG-Trefhaut)/(NTG-1);
  double hPG=(PrefGB-PrefGH)/(NTG-1);


  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a droite 
  for(int j=0; j<=N23-1; j++){ // N13 points
    /*
    T=T23[j];
    P=P23[j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    */
    V=segVc23[j];
    E=segEc23[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N23;

  // Parcours de la frontière droite
  
  for(int j=1; j<=NTD-1; j++){ 
    T=TrefD+j*hTD;
    P=PrefD+j*hPD;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }

  nref+=NTD-1;
  
  // Parcours de la frontière haute
  T=Trefhaut;
  for(int j=1; j<=NPH-1; j++){ 
    P=PrefDH+j*hPH;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }
  nref+=NPH-1;

  /*
      Calculer les points et les mettre sur la frontières
  */



  // Parcours de la frontière gauche
  for(int j=1; j<=NTG-2; j++){ // N13 points
    T=Trefhaut+j*hTG;
    P=PrefGH+j*hPG;
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    S[nref+j-1][0]=V;
    S[nref+j-1][1]=E;
  }
  nref+=NTG-2;

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N13-2; j++){ // N13 points
    /*
    T=T13[N13-1-j];
    P=P13[N13-1-j];
    err= fVE(phase, P, T, &V, &E);
    if(err){printf("Erreur de fVE (sommetsC)\n");return err;}
    */
    V=segVc13[N13-1-j];
    E=segEc13[N13-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }

  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple


  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsAB(int N12, double* segVa12, double* segEa12, double* segVb12, double* segEb12, double** S){
  
  int N=2*N12;

  int err=0;

  double V,E;
  

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N12-1; j++){ 
    V=segVb12[j];
    E=segEb12[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N12;
  

  //  Parcours de C12 de la gauche jusqu'au point triple 
  for(int j=0; j<=N12-1; j++){ // N13 points
    V=segVa12[N12-1-j];
    E=segEa12[N12-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple


  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsAC(int N13, double* segVa13, double* segEa13, double* segVc13, double* segEc13, double** S){

  int N=2*N13 ;

  int err=0;

  double V,E;

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C13 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N13-1; j++){ 
    V=segVc13[j];
    E=segEc13[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N13;
  

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N13-1; j++){ // N13 points
    V=segVa13[N13-1-j];
    E=segEa13[N13-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


// fonction qui remplie le tableau S des somments  
/* Arguments :
  - int N13.
  - int ND, NG : nombre de points sur les fornitères droites et gauches
  - double* P12, T12 : tableaux de N12 éléments qui contients les points de la courbe C12
                       partant du point triple

  Arguments de sorties: void
  - double** S : les sommets du polygone tendant vers la zone de mélange 12
*/
int sommetsBC(int N23, double* segVb23, double* segEb23, double* segVc23, double* segEc23, double** S){
  
  int N=2*N23;

  int err=0;

  double V,E;

  int nref=0;

  // **********   ECRITURE DE LA MATRICE DES SOMMETS   *************
  //  Parcours de C23 duu point triple jusqu'a gauche (dans le plans VE)
  for(int j=0; j<=N23-1; j++){ 
    V=segVc23[j];
    E=segEc23[j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }
  
  nref+=N23;
  

  //  Parcours de C13 de la gauche jusqu'au point triple 
  for(int j=0; j<=N23-1; j++){ // N13 points
    V=segVb23[N23-1-j];
    E=segEb23[N23-1-j];
    S[nref+j][0]=V;
    S[nref+j][1]=E;
  }


  S[N][0]=S[0][0];
  S[N][1]=S[0][1];
  
  return 0;
}


/*
   Renvoie 1 si (x,y) appartient au poly def par S
   et 0 sinon
*/
int intPOLY(double x, double y, int N, double** S){

  double m, pente;
  int cpt=0;
  
  double xi, yi;
  double xi1, yi1;
  double xI, yI;
  double xext, yext;
  double Yext, Xext;

  Yext=S[0][1];
  for(int i=1; i<N; i++){
    if(Yext<S[i][1]){
      Yext=S[i][1];
    }
  }
  Yext+=abs(floor(Yext/2)); // on s'assure que (., Yext) est hors du poly S
  
  xext=x;
  yext=Yext;
  
  for(int i=0; i<N; i++){
    xi=S[i][0]; xi1=S[i+1][0];
    yi=S[i][1]; yi1=S[i+1][1];
    if(xi!=xi1){
      
      pente=(yi1-yi)/(xi1-xi);
      xI=x;
      yI=pente*(x-xi)+yi;
      if( (xI-xi)*(xI-xi1)<0 &&  (yI-y)*(yI-yext)<0){  // (yI-yi)*(yI-yi1)<0
        cpt++;
        //printf("  #######   cpt +1   ###########\n");
      }
      /*
      printf("   --i= %d\n",i );
      printf("        V= %lf  ,   E= %lf\n",x*1e6,y );
      printf("        xI= %lf  ,   yI= %lf\n",xI*1e6,yI );
      printf("        x_ext= %lf  ,   y_ext= %lf\n",xext*1e6,yext );
      printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi)*1e6, (xI-xi1)*1e6 );
      printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi)*1e6, (yI-yi1) );
      printf("       (yI-y)= %lf , (yI-yext)= %lf\n",(yI-y)*1e6, (yI-yext) );
      printf("       Vi= %lf Ei= %lf   ||   Vi1= %lf Ei1= %lf\n",xi*1e6,yi ,xi1*1e6,yi1);
      */
    }
  }
  //printf("  * cpt= %d\n",cpt );
  if(cpt%2==1){return 1;} // (x,y) est à l'intérieur  entier impair
  else{return 0;}         // (x,y) extérieur

}



// Test si le segment [(xa,ya),(xb,yb)] s'intersecte avec un des segment du polygone S
/*
   Renvoie  1 si intersection
   et       0 si pas d'intersection 
           -1 en cas d'erreur 

*/
int inter_seg(double** REC, int N, double** S, int* pInd){
  
  double penteP, pente;
  double xi,yi,xi1,yi1;
  double xI,yI; // coordonées de l'intersection
  
  double xa, ya, xb, yb;

  for(int ind_rec=0; ind_rec<4; ind_rec++){
    //printf(" ind_rec= %d\n",ind_rec);

    xa=REC[ind_rec][0];    ya=REC[ind_rec][1];
    xb=REC[ind_rec+1][0];  yb=REC[ind_rec+1][1];

    if(xa==xb && ya==yb){printf("deux points du rectangle sont egaux (inter_seg/calfront.c)\n"); return -1;}

    //printf("    xa= %lf  ,   ya= %lf\n",xa,ya );
    //printf("    xb= %lf  ,   yb= %lf\n",xb,yb );

    if(xa==xb){ // pente du segment P non-défini
      xI=xa; // segement vertical
      for(int i=0; i<N; i++){
        xi=S[i][0]; xi1=S[i+1][0];
        yi=S[i][1]; yi1=S[i+1][1];
        
        //printf("   --i= %d\n",i );
        //printf("       xi= %lf  ,   yi= %lf\n",xi,yi );
        //printf("       xi1= %lf  ,   yi1= %lf\n",xi1,yi1 );
        
        if(xi1==xi){printf(" segment tout deux verticaux : xi1= %g, xi= %g, yi1= %g, yi= %g (inter_seg/calfront.c)\n",xi1,xi,yi1,yi); return 1;} // mettre une commande pour passer au i suivant 
        else{
          pente=(yi1-yi)/(xi1-xi);
          yI=pente*(xI-xi)+yi;  
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          if( (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 ){ // on sait que ya!=yb
            if(yi==yi1){*pInd=i; return 1;}
            else{
              if((yI-yi)*(yI-yi1)<0){*pInd=i; return 1;} 
            }
          }
        }  
      }
    }
    else{ // pente défini car xa=!xb
      penteP=(yb-ya)/(xb-xa);
      for(int i=0; i<N; i++){
        xi=S[i][0]; xi1=S[i+1][0];
        yi=S[i][1]; yi1=S[i+1][1];
        //printf("   --i= %d\n",i );
        //printf("        xi= %lf  ,  yi= %lf\n",xi,yi );
        //printf("        xi1= %lf ,  yi1= %lf\n",xi1,yi1 );
        if(xi1==xi){ // pas de pente pour le segment S
          xI=xi;
          yI=penteP*(xI-xa)+ya;
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          if( (xI-xa)*(xI-xb)<0  &&  (yI-yi)*(yI-yi1)<0 ){  // on sait que yi!=yi1
            if(ya==yb){*pInd=i; return 1;}
            else{
              if((yI-ya)*(yI-yb)<0){*pInd=i; return 1;}
            }
          }
        }
        else{ // les deux segment ont une pente  ouf !!
          pente=(yi1-yi)/(xi1-xi);
          xI=(yi-ya+penteP*xa-pente*xi)/(penteP-pente);
          yI=pente*(xI-xi)+yi; 
          //printf("       xI= %lf  ,   yI= %lf\n",xI,yI );
          //printf("       (xI-xi)= %lf , (xI-xi1)= %lf\n",(xI-xi), (xI-xi1) );
          //printf("       (yI-yi)= %lf , (yI-yi1)= %lf\n",(yI-yi), (yI-yi1) );
          //printf("       (xI-xa)= %lf , (xI-xb)= %lf\n",(xI-xa), (xI-xb) );
          //printf("       (yI-ya)= %lf , (yI-yb)= %lf\n",(yI-ya), (yI-yb) );
          //if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){return 1;}
          
          if(penteP==pente){} // segment parralel //
          else if(yi1==yi){
            if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 ){*pInd=i; return 1;}
          }
          else if(ya==yb){
            if( (xI-xa)*(xI-xb)<0 &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){*pInd=i; return 1;}
          }
          else{
            if( (xI-xa)*(xI-xb)<0 &&  (yI-ya)*(yI-yb)<0  &&  (xI-xi)*(xI-xi1)<0 &&  (yI-yi)*(yI-yi1)<0 ){*pInd=i; return 1;}
          }
        }
      }
    }
  }


  return 0;  // pas d'intsection
}


/*  
    Creation des fichiers de maillage
*/
int sommets_polygone(int NA, int NB, int NC, int NTG1, int NPB1, int NPB2, int NTD2, int NTD3, int NPH3, int NTG3, int N12, int N13, int N23,
                     double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC,
                     double* segP12, double* segT12, double* segP13, double* segT13, double* segP23, double* segT23,
                     double* segVa12, double* segEa12, double* segVb12, double* segEb12,
                     double* segVa13, double* segEa13, double* segVc13, double* segEc13,
                     double* segVb23, double* segEb23, double* segVc23, double* segEc23){

  int err;
  int Nmax=1e4;
  double epsilon=5e-5;

  double Vat,Vbt,Vct;
  double Eat,Ebt,Ect;
  double Ptriple, Ttriple;
  point_triple(&Ptriple, &Ttriple, &Vat, &Eat, &Vbt, &Ebt, &Vct, &Ect);
  double Pdeb=Ptriple;

  // Limitation de la ligne de chgmt de phase beta/gamma (A/B)
  int phaseA12=1;
  int phaseB12=2;
  //double Pfin12=14.62153e9; //intersecte T=0   ancien jeu de coeff
  double Pfin12=15.3e9;  // point d'intersection entre T=0 et la ligne de chgmt de phase beta/gamma
  
  // Limitation de la ligne de chgmt de phase beta/liquide (A/C)
  int phaseA13=1;
  int phaseB13=3;
  //double Pfin13=-7e9;   // 0
  double Pfin13=-5e9;   // 0
  

  // Limitation de la ligne de chgmt de phase gamma/liquide (B/C)
  int phaseA23=2;
  int phaseB23=3;
  //double Pfin23=90e9;
  double Pfin23=40e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;

  
  // Valeur d'initialisation
  Va012=Vat;  Ea012=Eat;
  Vb012=Vbt;  Eb012=Ebt;
  
  Va013=Vat;  Ea013=Eat;
  Vb013=Vct;  Eb013=Ect;
  
  Va023=Vbt;  Ea023=Ebt;
  Vb023=Vct;  Eb023=Ect;


  // frontiere 12
  //printf(" ligne 12  N12= %d\n",N12);
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, N12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);
  
  // frontiere 13
  //printf(" ligne 13   N13= %d\n",N13);
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, N13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  //printf(" ligne 23   N23= %d\n",N23);
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, N23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  
  
  //  Sommets de la phase BETA (A ds le code)
  //int NA=N13+(NTG1-1)+(NPB1-1)+(N12-2);
  //printf(" -- NA= %d\n",NA );
  double TrefBA=segT12[N12-1];
  double TrefHA=segT13[N13-1];
  double PrefGA=segP13[N13-1];
  double PrefDA=segP12[N12-1];
  err=sommetsA(N13, N12, NTG1, NPB1, PrefGA, PrefDA, TrefHA, TrefBA, segVa13, segEa13, segVa12, segEa12, SA);
  if(err){return err;}
  
  //   Sommets de la phase GAMMA (B ds le code)
  //int NB=N12+(NPB2-1)+(NTD2-1)+(N23-2);
  //printf(" -- NB= %d\n",NB );
  double TrefBB=segT12[N12-1];
  double TrefHB=segT23[N23-1];
  double PrefGB=segP12[N12-1];
  double PrefDB=segP23[N23-1];
  err=sommetsB(N12, N23, NPB2, NTD2, PrefGB, PrefDB, TrefBB, TrefHB, segVb12, segEb12, segVb23, segEb23, SB);
  if(err){return err;}
  
  //   Sommets de la phase LIQUIDE (C ds le code)
  //int NC=N23+(NTD3-1)+(NPH3-1)+(NTG3-1)+(N13-2);
  //printf(" -- NC= %d\n",NC );
  /*
  double TrefhautC=25e3;//25e3;
  double PrefGHC=71e9;  
  double PrefDHC=170e9;
  */
  
  // 2 points sur la droite affine qui est la frontière de la zone de non-definition de la phase liquide (sur la gauche du diagramme)
  // permet de deduire la pression minimal en fonction d'une temperature voulue  (min tout en etant sur du domaine de def de la phase liquide)
  double y1=segT13[N13-1], x1=segP13[N13-1];
  double y2=5e3, x2=10e9;
  // coeff de la droite affine
  double a=(y1-y2)/(x1-x2);
  double b=(x2*y1 - x1*y2)/(x2-x1);
  double TrefhautC=3.5e3;//25e3;
  double PrefGHC=(TrefhautC-b)/a;  
  printf("PrefGHC=%g\n",PrefGHC);
  //PrefGHC=70e9;  
  double PrefDHC=segP23[N23-1];
  double TrefDC=segT23[N23-1];
  double TrefGC=segT13[N13-1];
  double PrefGBC=segP13[N13-1];
  double PrefDC=segP23[N23-1];
  err=sommetsC(N13, N23, NTD3, NPH3, NTG3, PrefGBC,PrefGHC,PrefDHC, PrefDC, TrefGC, TrefDC, TrefhautC, segVc13, segEc13, segVc23, segEc23, SC);
  if(err){return err;}
  

  //printf("Sommets AB\n");
  err=sommetsAB(N12, segVa12, segEa12, segVb12, segEb12, SAB);
  if(err){return err;}
  //printf("Sommets AC\n");
  err=sommetsAC(N13, segVa13, segEa13, segVc13, segEc13, SAC);
  if(err){return err;}
  //printf("Sommets BC\n");
  err=sommetsBC(N23, segVb23, segEb23, segVc23, segEc23, SBC);
  if(err){return err;}
  

  
  return 0;
}



/* Valeur d'entree :
    - Va,Ea, Vb,Eb
    - Vmel, Emel

    - A,InvA,dX : pour limiter les allocations 

  Valeur sortie
    - Vx,Ex, Vy,Ey, 
    - x : fraction massique de la phase Y

   double* A et InvA permet une seule allocation 
*/
int Newton_mixte(double* A, double* InvA, double* Delta, int phaseX, int phaseY, double Vmel, double Emel, double Vx, double Ex, double Vy, double Ey, double x, double* pVX, double* pEX, double* pVY, double* pEY, double* pX, double* dX){
  int dim=5;

  // Grandeurs thermodynamiques
  double Px, Tx, Sx, Gx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, Sy, Gy, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  
  
  //printf("Px= %.10lf, Tx= %.10lf, Sx= %.10lf, Gx= %.10lf\n",Px,Tx,Sx,Gx );
  //printf("Py= %.10lf, Ty= %.10lf, Sy= %.10lf, Gy= %.10lf\n",Py,Ty,Sy,Gy );

  
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
   
  double correc=1.;

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
int VE_mixte(int Nmax, double epsilon, int phaseX, int phaseY, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, 
             double* pVx, double* pEx, double* pVy, double* pEy, double* px, int* n_iter, double* pcritere){
  int err, n=1, dim=5;
  double errV, errE;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* dX=malloc(dim*sizeof(double));
  double* Delta=malloc(dim*sizeof(double));
  
  double Px, Tx, Sx, Gx;
  double Py, Ty, Sy, Gy;

  double Vx, Ex, Vy, Ey, x;
  double VX, EX, VY, EY, X;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  x=1./2;

  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  Sx=fS(Vx,Ex, phaseX); Gx=Ex+Px*Vx-Tx*Sx;
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  Sy=fS(Vy,Ey, phaseY); Gy=Ey+Py*Vy-Ty*Sy;


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){

    err=Newton_mixte(A, InvA, Delta, phaseX, phaseY, Vmel, Emel, Vx, Ex, Vy, Ey, x, &VX, &EX, &VY, &EY, &X, dX);

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    x=X;
    errV=((1-x)*Vx + x*Vy-Vmel)/Vmel;
    errE=((1-x)*Ex + x*Ey-Emel)/Emel;

    //critere=fabs(dX[0]/Vx) +fabs(dX[1]/Ex) +fabs(dX[2]/Vy) +fabs(dX[3]/Ey)+fabs(dX[4]/x);
    //critere=fmax(fabs(errV),fabs(errE
    critere=fabs(dX[4]);

    n++;
  }
  //printf(" --n= %d, critere= %g\n",n,critere );
  //for(int i=0; i<dim; i++){  printf(" dX[%d]= %g  ",i,dX[i] );  } printf("\n");
  
  if( isnan(critere) ){printf("Erreur critere=%g , V=%g; E=%g; zone: %d%d (VE_mixte)\n ",critere,Vmel,Emel,phaseX, phaseY ); return 1;}
  if(n==Nmax){printf("Newton_mixte non convergent : -n= %d (Nmax=%d) , critere= %g V=%g, E=%g (VE_mixte)\n",n,Nmax,critere, Vmel, Emel );return 1;}
  *n_iter=n;

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *px=x;
   
   *pcritere=critere;
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



// Aire du parallélépipède porté par v1 et v2 =======================================
double aire(double* v1, double* v2){
  return v1[0]*v2[1] - v1[1]*v2[0];
}



//  coorbary  =======================================================================
/*
Renvoie les 3 coordonées barycentrique de (x,y) pour un triangle du maillage
Arguments d'entree :
 - double *xy        : coordonnées d'un point de la discretisation
 - double **sommets  : coordonnées des trois points consistuant un triangle du maillage

Arguments de sortie :
 - double **pcoorbar : pointeur sur coorbar, tableau contenant les trois coordonnées barycentrique du pts
*/
int coorbary(double *xy, double **sommets, double **pcoorbar){
  double epsilon=1e-7;
  double ajax[2], akai[2], ajai[2], akaj[2];
  double *coorbar=*pcoorbar;
  double a1a2a3;
  // (i,j,k)=(1,2,3) à une permutation circulaire pres
  // ai, aj, ak les trois sommets du triangle
  // ajx= vecteur(ak,xy) aji=vecteur(aj,ai)
  int i,j,k;
  for(i=0; i<2; i++){
    j=(i+1)%3; // permutation circulaire de (i,j,k)
    k=(i+2)%3;
    for(int d=0; d<2; d++){
      ajax[d]=xy[d]-sommets[j][d];
      akaj[d]=sommets[j][d]-sommets[k][d];
      akai[d]=sommets[i][d]-sommets[k][d];
      ajai[d]=sommets[i][d]-sommets[j][d];
    }
    a1a2a3=aire(ajai,akai);
    if(fabs(a1a2a3)<epsilon){ return 1;}
    coorbar[i]=aire(ajax,akaj)/a1a2a3;
  }
  coorbar[2]=1-coorbar[0]-coorbar[1];
  return 0;
}


/*
  Fonction qui renvoie les deux indices(corespondant à des droites) entre lesquelles se trouve le couple (V,E)
  --> le but étant de trouver la pression et la température en ce point
*/
int interp_zone_mixte(double V, double E, int phaseALPHA, int phaseBETA, int N, double* Valpha, double* Vbeta, double* Ealpha, double* Ebeta, double* pP, double* pT, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
    double P,T;
    int res;
    double epsilon=1e-5;

    // reperage du quadrilatere
    double** S=alloctabd(5,2);
    int ind=0;
    int testNO=1;
    for(int i=0; i<=N-2; i++){
      S[0][0]=Valpha[i];   S[0][1]=Ealpha[i];
      S[1][0]=Valpha[i+1]; S[1][1]=Ealpha[i+1];
      S[2][0]=Vbeta[i+1];  S[2][1]=Ebeta[i+1];
      S[3][0]=Vbeta[i];    S[3][1]=Ebeta[i];
      S[4][0]=Valpha[i];   S[4][1]=Ealpha[i];

      if(intPOLY(V, E, 4, S)){ind=i; testNO=0;break;} // ou i=N
    }

    if(testNO){
      printf("\n");
      printf("le point est dans aucun des quadrilatere\n");
      printf(" ind= %d\n",ind );
      printf("V= %.10lf , E= %.10lf\n",V*1e6,E); 
      return 1;
    }

    *pVmx=Valpha[ind];  *pEmx=Ealpha[ind];
    *pVmy=Vbeta[ind];   *pEmy=Ebeta[ind];

    
    double pente1, pente2;

    // Calcul du poit milieu
    // diag de [Malpha_i, Mbeta_i+1]
    pente1=(Ebeta[ind+1]-Ealpha[ind])/(Vbeta[ind+1]-Valpha[ind]);
    //Ediag1=Ealpha[i] + pente1*(V-Valpha[i]);

    // diag de [Malpha_i+1, Mbeta_i]
    pente2=(Ebeta[ind]-Ealpha[ind+1])/(Vbeta[ind]-Valpha[ind+1]);
    //Ediag2=Ealpha[i+1] + pente2*(V-Valpha[i+1]);
    
    // On cherche V* tq     Ealpha[i+1] + pente2*(V*-Valpha[i+1]) = Ealpha[i] + pente1*(V*-Valpha[i])
    double VI=(Ealpha[ind]-Ealpha[ind+1] + pente2*Valpha[ind+1] - pente1*Valpha[ind])/(pente2-pente1);
    double EI=Ealpha[ind+1] + pente2*(VI-Valpha[ind+1]);

    //printf("  VI= %.10lf,  EI= %f\n",VI,EI );

    // Coordonées d'un triangle
    double** coorTri=alloctabd(3,2);

    // Calcul des 3 coordonnées barycentriques pour les 4 triangles
    double* coorBar=malloc(3*sizeof(double));  // une ligne = un triangle

    double* VE=malloc(2*sizeof(double));  // une ligne = un triangle
    VE[0]=V; VE[1]=E;


    // Calcul des pressions et températures aux sommets des triangles
    double* Ptri=malloc(3*sizeof(double));  // une ligne = un triangle
    double* Ttri=malloc(3*sizeof(double));  // une ligne = un triangle
    
    int test0, test1, test2;    

     
    //printf("  $$ 1er triangle  $$\n");
    // 1er triangle
    coorTri[0][0]=Valpha[ind];    coorTri[0][1]=Ealpha[ind];
    coorTri[1][0]=Valpha[ind+1];  coorTri[1][1]=Ealpha[ind+1];
    coorTri[2][0]=Vbeta[ind+1];   coorTri[2][1]=Ebeta[ind+1];

    Ptri[0]=fP(coorTri[0][0], coorTri[0][1], phaseALPHA);  
    Ptri[1]=fP(coorTri[1][0], coorTri[1][1], phaseALPHA);  
    Ptri[2]=fP(coorTri[2][0], coorTri[2][1], phaseBETA);  

    Ttri[0]=fT(coorTri[0][0], coorTri[0][1], phaseALPHA);  
    Ttri[1]=fT(coorTri[1][0], coorTri[1][1], phaseALPHA);  
    Ttri[2]=fT(coorTri[2][0], coorTri[2][1], phaseBETA);  


    res=coorbary(VE, coorTri, &coorBar);
    
    /*
    printf(" Triangle 1 :\n");
    for(int i=0; i<3; i++){
      printf(" coorBar[%d]= %lf",i,coorBar[i] );
    }
    printf("\n");
    */

    test0=( 0-epsilon<=coorBar[0] && coorBar[0]<=1+epsilon );
    test1=( 0-epsilon<=coorBar[1] && coorBar[1]<=1+epsilon );
    test2=( 0-epsilon<=coorBar[2] && coorBar[2]<=1+epsilon );
    if(test0 && test1 && test2){  // on est dans ce triangle
      P=0; T=0;
      for(int i=0; i<3; i++){
        P+=coorBar[i]*Ptri[i];  
        T+=coorBar[i]*Ttri[i];
      }
      *pP=P;
      *pT=T;
      return 0;
    }
    else{

      //printf("  $$ 2nd triangle  $$\n");
      // 2nd triangle
      coorTri[0][0]=Valpha[ind];    coorTri[0][1]=Ealpha[ind];
      coorTri[1][0]=Vbeta[ind];     coorTri[1][1]=Ebeta[ind];
      coorTri[2][0]=Vbeta[ind+1];   coorTri[2][1]=Ebeta[ind+1];

      Ptri[0]=fP(coorTri[0][0], coorTri[0][1], phaseALPHA);  
      Ptri[1]=fP(coorTri[1][0], coorTri[1][1], phaseBETA);  
      Ptri[2]=fP(coorTri[2][0], coorTri[2][1], phaseBETA);  

      Ttri[0]=fT(coorTri[0][0], coorTri[0][1], phaseALPHA);  
      Ttri[1]=fT(coorTri[1][0], coorTri[1][1], phaseBETA);  
      Ttri[2]=fT(coorTri[2][0], coorTri[2][1], phaseBETA);  

     
      res=coorbary(VE, coorTri, &coorBar);
      
      /*
      printf(" Triangle 2 :\n");
      for(int i=0; i<3; i++){
        printf(" coorBar[%d]= %lf",i,coorBar[i] );
      }
      printf("\n");
      */

      test0=( 0-epsilon<=coorBar[0] && coorBar[0]<=1+epsilon );
      test1=( 0-epsilon<=coorBar[1] && coorBar[1]<=1+epsilon );
      test2=( 0-epsilon<=coorBar[2] && coorBar[2]<=1+epsilon );
      if(test0 && test1 && test2){  // on est dans ce triangle
        P=0; T=0;
        for(int i=0; i<3; i++){
          P+=coorBar[i]*Ptri[i];  
          T+=coorBar[i]*Ttri[i];  
        }
        *pP=P;
        *pT=T;
        return 0;
      }
      else{
        printf("\nLe triangle n'a pas été trouvé (interp_zone_mixte)\n");
        
        printf("  ind= %d       N-1= %d\n",ind,N-1 );
        printf("V= %.10lf, E= %lf\n",V,E );
        printf("  4 sommets du quadrilatere :\n");
        printf("Valpha[%d]= %.10lf , Ealpha[%d]= %lf\n", ind, Valpha[ind], ind, Ealpha[ind]);
        printf("Valpha[%d]= %.10lf , Ealpha[%d]= %lf\n", ind+1, Valpha[ind+1], ind+1, Ealpha[ind+1]);
        printf("Vbeta[%d]= %.10lf , Ebeta[%d]= %lf\n", ind+1, Vbeta[ind+1], ind+1, Ebeta[ind+1]);
        printf("Vbeta[%d]= %.10lf , Ebeta[%d]= %lf\n", ind, Vbeta[ind], ind, Ebeta[ind]);
    
        for(int i=0; i<3; i++){
          printf("    coorBar[%d]= %lf\n",i,coorBar[i] );
        }
        return 1;  
      }
    }

  free(S);
  free(Ptri); free(Ttri);
  freetab(coorTri); 
  free(coorBar);
  free(VE);
}



/*
   Calcule dans quel quadrilatere d'une zone mixte donné se trouve le couple (V,E)
*/
int graine_zone_mixte(double V, double E, int N, double* Valpha, double* Vbeta, double* Ealpha, double* Ebeta, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
    int res;

    // reparage du quadrilatere
    double** S=alloctabd(5,2);
    int ind=0;
    int testNO=1;
    for(int i=0; i<=N-2; i++){
      S[0][0]=Valpha[i];   S[0][1]=Ealpha[i];
      S[1][0]=Valpha[i+1]; S[1][1]=Ealpha[i+1];
      S[2][0]=Vbeta[i+1];  S[2][1]=Ebeta[i+1];
      S[3][0]=Vbeta[i];    S[3][1]=Ebeta[i];
      S[4][0]=Valpha[i];   S[4][1]=Ealpha[i];

      if(intPOLY(V, E, 4, S)){ind=i; testNO=0;break;} // ou i=N
    }

    if(testNO){
      printf("\n");
      printf("le point est dans aucun des quadrilatere (graine_zone_mixte)\n");
      printf(" ind= %d\n",ind );
      printf("V= %.10lf , E= %.10lf\n",V*1e6,E); 
      return 1;
    }

    *pVmx=Valpha[ind];  *pEmx=Ealpha[ind];
    *pVmy=Vbeta[ind];   *pEmy=Ebeta[ind];

    


  free(S);
}



// Calcul de la célérite du son
/*  Pp=P(V+dV,E+dE(+dV) )  
    Pm=P(V-dV,E+dE(-dV) )  

    Argument de sortie :
     - double* pc : pointeur sur la celerité du son
*/
int c_celerite_son(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                   double V0x, double E0x, double V0y, double E0y, double *pc, double *pG){

  int err, n_iter;
  double critere;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,x;

  int phaseA,phaseB;

  if(zone==-1){printf("Erreur out of domain : zone= %d (c_celerite_son)\n",zone); return 1;  }
  else if(zone==1 || zone==2 || zone==3){  
    Pp=fP(V+dV,E+dE,zone); 
    Pm=fP(V-dV,E-dE,zone); 
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseB=zone%10;
    phaseA=(zone-phaseB)/10;
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  if(err){return 1;}
    Pp=fP(Vx,Ex,phaseA);

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  if(err){return 1;}
    Pm=fP(Vx,Ex,phaseA);  
  }
  else if(zone==0){
    *pc=0;
    return 0;
  }

  //printf("zone=%d   Pm=%g  Pp=%g\n",zone,Pm,Pp );

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV;  

  double d2PdV2=(Pm+Pp-2*P)/(dV*dV);
  *pG=-(V/2)*d2PdV2/dPdV;
  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g\n",V*1e6,E );  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}



// Fonction pour trier les zone par ordre croisant
int fdiff(const void* A, const void* B){
  //int x=*(*A);
  //int y=*(*B);

  int const *pa = A;
  int const *pb = B;

  return *pb-*pa;
} 
// Fonction pour avoir une benchmark de (P,T) sur (V,E)
/*

*/
int ecr_benchmark_PT(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax){

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisation juste en lisant les fichiers
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fprintf(fparams, "%d \n",nb_points );
  fprintf(fparams, "%d %d\n",Nv,Ne );
  fprintf(fparams, "%.15lf %.15lf\n",Vmin,Vmax );
  fprintf(fparams, "%.15lf %.15lf\n",Emin,Emax );
  fprintf(fparams, "\n\n Organisation du fichier\n");
  fprintf(fparams, "nb_points\n");
  fprintf(fparams, "Nv_dis  Ne_dis\n");
  fprintf(fparams, "Vmin_dis Vmax_dis\n");
  fprintf(fparams, "Emin_dis Emax_dis\n");


  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;

  int err,res;
  double E,V;

  double Vmx, Emx;
  double Vmy, Emy;


  // Variables pour discretisation de la maille
  double V_sub, E_sub;
  int N_rec=10;
  double he_sub=he/(N_rec-1);
  double hv_sub=hv/(N_rec-1);
  double** VE_front_rec=alloctabd(4*(N_rec-1),2);
  int count;
  

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  printf("NA=%d   NB=%d   NC=%d\n",NA,NB,NC );
  printf("N12=%d   N13=%d   N23=%d\n",N12,N13,N23 );

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);
  
   
  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));
  double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));
  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));
  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));
  double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));
  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));
  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));
  double* segEc23=malloc(N23*sizeof(double));


  // CALCUL DES SOMMETS  
  printf(" Sommets_polygone \n");
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}
  
  // Calcul des coord du triangle triple
  double** STRI_12=alloctabd(2,2);
  double** STRI_13=alloctabd(2,2);
  double** STRI_23=alloctabd(2,2);

  double Va_temp,Vb_temp,Vc_temp;
  double Ea_temp,Eb_temp,Ec_temp;
  double PT,TT;
  point_triple(&PT, &TT, &Va_temp, &Ea_temp, &Vb_temp, &Eb_temp, &Vc_temp, &Ec_temp);
  STRI_12[0][0]=Va_temp;   STRI_12[0][1]=Ea_temp;
  STRI_12[1][0]=Vb_temp;   STRI_12[1][1]=Eb_temp;
  
  STRI_13[0][0]=Va_temp;   STRI_13[0][1]=Ea_temp;
  STRI_13[1][0]=Vc_temp;   STRI_13[1][1]=Ec_temp;
  
  STRI_23[0][0]=Vb_temp;   STRI_23[0][1]=Eb_temp;
  STRI_23[1][0]=Vc_temp;   STRI_23[1][1]=Ec_temp;

  FILE *ftriple;
  if((ftriple = fopen("fichiers/triple.txt", "w+")) == NULL){printf("erreur ouverture fichier triple.txt \n"); return 1;}
  fprintf(ftriple, "%g %g\n",Va_temp, Ea_temp );
  fprintf(ftriple, "%g %g\n",Vb_temp, Eb_temp );
  fprintf(ftriple, "%g %g\n",Vc_temp, Ec_temp );
  fprintf(ftriple, "%g %g\n",PT,TT );
  /*
  fprintf(ftriple, "\n\n Organisation du fichier\n");
  fprintf(ftriple, "VbetaT EbetaT\n");
  fprintf(ftriple, "VgammaT EgammaT\n");
  fprintf(ftriple, "VliqT EliqT\n");
  fprintf(ftriple, "Ptriple Ttriple\n");
  */
  fclose(ftriple);

  // Test affichage STRI
  /*
  printf(" STRI_12= [ %g  %g \n            %g  %g ]\n",STRI_12[0][0],STRI_12[0][1],STRI_12[1][0],STRI_12[1][1] );
  printf(" STRI_13= [ %g  %g \n            %g  %g ]\n",STRI_13[0][0],STRI_13[0][1],STRI_13[1][0],STRI_13[1][1] );
  printf(" STRI_23= [ %g  %g \n            %g  %g ]\n",STRI_23[0][0],STRI_23[0][1],STRI_23[1][0],STRI_23[1][1] );
  */
  
  /*
  STRI_12[0][0]=129.905671875844291e-6;   STRI_12[0][1]=72.168190127265447e3;
  STRI_12[1][0]=127.666112185063994e-6;   STRI_12[1][1]=99.480256119236515e3;
  
  STRI_13[0][0]=129.905671875844291e-6;   STRI_13[0][1]=72.168190127265447e3;
  STRI_13[1][0]=132.538478805442338e-6;   STRI_13[1][1]=136.154346228414312e3;
  
  STRI_23[0][0]=127.666112185063994e-6;   STRI_23[0][1]=99.480256119236515e3;
  STRI_23[1][0]=132.538478805442338e-6;   STRI_23[1][1]=136.154346228414312e3;
  */

  // Coordonnées de la ligne Sa12  (sépartion de A et de la zone mixte 12 ou A/B)
  //  et de Sb12
  double** Sa12=alloctabd(N12,2);  
  double** Sb12=alloctabd(N12,2); 
  //printf("\n\nN12=%d\n", N12);
  for(int i=0; i<N12; i++){
    Sa12[i][0]=segVa12[i]; 
    Sa12[i][1]=segEa12[i];

    Sb12[i][0]=segVb12[i]; 
    Sb12[i][1]=segEb12[i]; 

    //printf("Sa12=[%g, %g]\n",Sa12[i][0],Sa12[i][1]); 
    //printf("Sb12=[%g, %g]\n",Sb12[i][0],Sb12[i][1]); 
  }

  // Coordonnées de la ligne Sa13  (sépartion de A et de la zone mixte 13 ou A/C)
  //  et de Sb12
  double** Sa13=alloctabd(N13,2);  
  double** Sc13=alloctabd(N13,2); 
  //printf("\n\nN13=%d\n",N13 );
  for(int i=0; i<N13; i++){
    Sa13[i][0]=segVa13[i];
    Sa13[i][1]=segEa13[i];

    Sc13[i][0]=segVc13[i]; 
    Sc13[i][1]=segEc13[i];

    //printf("Sa13=[%g, %g]\n",Sa13[i][0],Sa13[i][1]); 
    //printf("Sc13=[%g, %g]\n",Sc13[i][0],Sc13[i][1]); 
  }

  // Coordonnées de la ligne Sa13  (sépartion de A et de la zone mixte 13 ou A/C)
  //  et de Sb12
  double** Sb23=alloctabd(N23,2);  
  double** Sc23=alloctabd(N23,2); 
  //printf("\n\nN23=%d\n",N23 );
  for(int i=0; i<N23; i++){
    Sb23[i][0]=segVb23[i]; 
    Sb23[i][1]=segEb23[i];

    Sc23[i][0]=segVc23[i]; 
    Sc23[i][1]=segEc23[i];  

    //printf("Sb23=[%g, %g]\n",Sb23[i][0],Sb23[i][1]); 
    //printf("Sc23=[%g, %g]\n",Sc23[i][0],Sc23[i][1]); 
  }

  // frontières extérieures des phases
  int NextA=NA-N12-N13+3;
  double** SextA=alloctabd(NextA,2); 
  //printf("\n\nNextA=%d\n",NextA ); 
  for(int i=0; i<NextA; i++){
    SextA[i][0]=SA[N13-1+i][0]; 
    SextA[i][1]=SA[N13-1+i][1];
    //printf("SextA=[%g, %g]\n",SextA[i][0],SextA[i][1]);   
  }
  /*
  printf(" VextA : \n[ ");
  for(int i=0; i<NextA; i++){    printf("%g ",SextA[i][0]);    } printf(" ];\n");
  printf(" EextA : \n[ ");
  for(int i=0; i<NextA; i++){    printf("%g ",SextA[i][1]);    } printf(" ];\n");
  */

  int NextB=NB-N12-N23+3;
  double** SextB=alloctabd(NextB,2);  
  //printf("\n\nNextB=%d\n",NextB );
  for(int i=0; i<NextB; i++){
    SextB[i][0]=SB[N12-1+i][0]; 
    SextB[i][1]=SB[N12-1+i][1];
    //printf("SextB=[%g, %g]\n",SextB[i][0],SextB[i][1]);
  }

  int NextC=NC-N23-N13+3;
  double** SextC=alloctabd(NextC,2);  
  //printf("\n\nNextC=%d\n",NextC );
  for(int i=0; i<NextC; i++){
    SextC[i][0]=SC[N23-1+i][0]; 
    SextC[i][1]=SC[N23-1+i][1]; 
    //printf("SextC=[%g, %g]\n",SextC[i][0],SextC[i][1]); 
  }

  // frontières extérieures des zones biphasiques
  double** SextAB=alloctabd(2,2);  
  SextAB[0][0]=Sa12[N12-1][0];
  SextAB[0][1]=Sa12[N12-1][1];
  
  SextAB[1][0]=Sb12[N12-1][0];  
  SextAB[1][1]=Sb12[N12-1][1]; 
 
  double** SextAC=alloctabd(2,2);  
  SextAC[0][0]=Sa13[N13-1][0];
  SextAC[0][1]=Sa13[N13-1][1];
  
  SextAC[1][0]=Sc13[N13-1][0];  
  SextAC[1][1]=Sc13[N13-1][1]; 

  double** SextBC=alloctabd(2,2);  
  SextBC[0][0]=Sb23[N23-1][0];
  SextBC[0][1]=Sb23[N23-1][1];
  
  SextBC[1][0]=Sc23[N23-1][0];  
  SextBC[1][1]=Sc23[N23-1][1]; 
  
  
  // Pour tester les valeurs (coord des points du maillage)
  FILE *fileV, *fileE;
  if((fileV = fopen("fichiers/Vbench.txt", "w+")) == NULL){printf("erreur ouverture fichier V bench \n"); return 1;}
  if((fileE = fopen("fichiers/Ebench.txt", "w+")) == NULL){printf("erreur ouverture fichier E bench \n"); return 1;}
  for(int j=0; j<Ne; j++){
    E=Emin+j*he;
    fprintf(fileE, "%.15lf ",E);
  }
  for(int j=0; j<Nv; j++){
    V=Vmin+j*hv;
    fprintf(fileV, "%.15lf ",V);
  }
  fclose(fileV);  fclose(fileE);
  /////////////
  


  //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif

  //  fichier texte diagramme (V,E)
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.10g %.10g\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.10g %.10g\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.10g %.10g\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10g %.10g\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10g %.10g\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10g %.10g\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
  

  // Fichier texte diagramme (P,T)
  FILE *fSApt, *fSBpt, *fSCpt, *fSABpt, *fSACpt, *fSBCpt;
  if((fSApt = fopen("fichiers/pSApt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSBpt = fopen("fichiers/pSBpt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSCpt = fopen("fichiers/pSCpt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSABpt = fopen("fichiers/pSABpt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSACpt = fopen("fichiers/pSACpt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBCpt = fopen("fichiers/pSBCpt.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSApt,"%.10g %.10g\n", fP(SA[i][0],SA[i][1],1),fT(SA[i][0],SA[i][1],1) );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSBpt,"%.10g %.10g\n", fP(SB[i][0],SB[i][1],2),fT(SB[i][0],SB[i][1],2) );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSCpt,"%.10g %.10g\n", fP(SC[i][0],SC[i][1],3),fT(SC[i][0],SC[i][1],3) );  }

  for(int i=0; i<N12; i++){  fprintf(fSABpt,"%.10g %.10g\n",fP(Sa12[i][0],Sa12[i][1],1),fT(Sa12[i][0],Sa12[i][1],1) );  }
  for(int i=0; i<N13; i++){  fprintf(fSACpt,"%.10g %.10g\n",fP(Sa13[i][0],Sa13[i][1],1),fT(Sa13[i][0],Sa13[i][1],1) );  }
  for(int i=0; i<N23; i++){  fprintf(fSBCpt,"%.10g %.10g\n",fP(Sc23[i][0],Sc23[i][1],3),fT(Sc23[i][0],Sc23[i][1],3) );  }

  fclose(fSApt);   fclose(fSBpt);   fclose(fSCpt);  fclose(fSABpt);  fclose(fSACpt);  fclose(fSBCpt);
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////



  double* tabVx=malloc((Nv-1)*(Ne-1)*sizeof(double));  double* tabEx=malloc((Nv-1)*(Ne-1)*sizeof(double));
  double* tabVy=malloc((Nv-1)*(Ne-1)*sizeof(double));  double* tabEy=malloc((Nv-1)*(Ne-1)*sizeof(double));
  
  int* tabNbZONE=malloc((Nv-1)*(Ne-1)*sizeof(int)); // nombre de zones présentent sur le triangle
  int** tabZONE=alloctab(5,(Nv-1)*(Ne-1));  // mémorise les zones présentent dans le rectangle
  int** tabINTER=alloctab(9,(Nv-1)*(Ne-1)); // mémorise les lignes avec lesquelles il y a une intersection
  int** tabInd_Mixte=alloctab(9,(Nv-1)*(Ne-1)); // indice du segment de la ligne de chgmt de phase qui s'intersecte avec le rectangle
  

  // une ligne un carré puis sur les colonnes, les zones qui le composent on rempli par -99 si rien dans les zones 
  // (ou 0 peut-être plus simple à implementer mais il faudra changer l'indice du tri) 

  int** numZONE=alloctab(234,2);
  // zone pure/zone mixte
  numZONE[121][0]=1;  numZONE[121][1]=12;
  numZONE[122][0]=2;  numZONE[122][1]=12;
  numZONE[131][0]=1;  numZONE[131][1]=13;
  numZONE[133][0]=3;  numZONE[133][1]=13;
  numZONE[232][0]=2;  numZONE[232][1]=23;
  numZONE[233][0]=3;  numZONE[233][1]=23;
  // triangle/zone mixte
  numZONE[120][0]=0;  numZONE[120][1]=12;
  numZONE[130][0]=0;  numZONE[130][1]=13;
  numZONE[230][0]=0;  numZONE[230][1]=23;
  // zone pure/hors domaine
  numZONE[11][0]=1;   numZONE[11][1]=-1;
  numZONE[12][0]=2;   numZONE[12][1]=-1;
  numZONE[13][0]=3;   numZONE[13][1]=-1;
  // zone mixte/hors domaine
  numZONE[112][0]=-1; numZONE[112][1]=12;
  numZONE[113][0]=-1; numZONE[113][1]=13;
  numZONE[123][0]=-1; numZONE[123][1]=23;


  double** REC=alloctabd(5,2); // rectangle du cadrillage
  
  int ind_mixte;
  int nb_zone;
  int nb_inter=0;
  int ind_inter;
  int test;

  int k=0;
  for(int i=0; i<Nv-1; i++){
    printf("  --i= %d\n",i );
    V=Vmin+i*hv + hv/2; // pour prendre la valeur au centre du rectangle
    for(int j=0; j<Ne-1; j++){
      E=Emin+j*he + he/2;

      //printf("  -- i= %d, j= %d",i,j );
      //printf("  V1=%g; V2=%g; E1=%g; E2=%g; \n",Vmin+i*hv,Vmin+(i+1)*hv,Emin+j*he,Emin+(j+1)*he );
      
      nb_inter=0;

      //   Calcul du chevanchement des zones
      REC[0][0]=Vmin+i*hv;       REC[0][1]=Emin+j*he;
      REC[1][0]=Vmin+(i+1)*hv;   REC[1][1]=Emin+j*he;
      REC[2][0]=Vmin+(i+1)*hv;   REC[2][1]=Emin+(j+1)*he;
      REC[3][0]=Vmin+i*hv;       REC[3][1]=Emin+(j+1)*he;
      REC[4][0]=REC[0][0];       REC[4][1]=REC[0][1];


      res=inter_seg(REC, N12-1, Sa12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=121; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N12-1, Sb12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=122; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}

      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N13-1, Sa13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=131; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N13-1, Sc13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=133; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N23-1, Sb23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=232; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, N23-1, Sc23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=233; tabInd_Mixte[nb_inter][k]=ind_mixte; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      
      // Pour les frontières triangles, on mets ind_mixte=0 pour avoir les points triples comme graine
      res=inter_seg(REC, 1, STRI_12, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=120; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, STRI_13, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=130; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, STRI_23, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=230; ind_mixte=0; nb_inter++;}
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextA-1, SextA, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=11; ind_mixte=0; nb_inter++; }
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextB-1, SextB, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=12; ind_mixte=0; nb_inter++; }
      //else if(res==0){tabZONE[1][k]=-99}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, NextC-1, SextC, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=13; ind_mixte=0; nb_inter++;}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      
      res=inter_seg(REC, 1, SextAB, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=112; ind_mixte=N12-1; nb_inter++; }
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, SextAC, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=113; ind_mixte=N13-1; nb_inter++;}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}

      res=inter_seg(REC, 1, SextBC, &ind_mixte);
      if(res==1){tabINTER[nb_inter][k]=123; ind_mixte=N23-1; nb_inter++;}
      else if(res==-1){printf("erreur de inter_seg (ecr_benchmark_PT)\n"); return res;}
      

      /*
      printf("   Nb intersection= %d -- [",nb_inter );
      for(int i=0; i<nb_inter; i++){  printf("  %d ",tabINTER[i][k] ); }
      printf("]\n\n");
      */
      
      // Cas de rectangle sans chevauchement
      if(nb_inter==0){
        nb_zone=1;
        tabNbZONE[k]=1;
        // Attribution d'une zone aux points milieux des rectangles
        if(intPOLY(V,E,NC,SC)){   tabZONE[0][k]=3;  }
        else if(intPOLY(V,E,NA,SA)){  tabZONE[0][k]=1;  }
        else if(intPOLY(V,E,NB,SB)){  tabZONE[0][k]=2;  }
        else if(ftri(V,E)){  tabZONE[0][k]=0;  }
        else if(intPOLY(V,E,NAB,SAB)){
          err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=12;
        }
        else if(intPOLY(V,E,NAC,SAC)){
          err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=13;
        }
        else if(intPOLY(V,E,NBC,SBC)){
          err=graine_zone_mixte(V, E, N23, segVb23, segVc23, segEb23, segEc23, &Vmx, &Emx, &Vmy, &Emy);
          if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
          tabVx[k]=Vmx; tabEx[k]=Emx;
          tabVy[k]=Vmy; tabEy[k]=Emy;
          tabZONE[0][k]=23;
        }
        else{  // hors domaine --> on ne se soucis des pb de chevauchement des rectangles sur les zones
          tabZONE[0][k]=-1;
        }
      }
      // ## 1 intersection donc seulement deux zones ##
      else if(nb_inter==1){
        nb_zone=2;
        tabNbZONE[k]=2;
        tabZONE[0][k]=numZONE[ tabINTER[0][k] ][0];
        tabZONE[1][k]=numZONE[ tabINTER[0][k] ][1];
        
        //debug
        if(tabZONE[0][k]==3 && tabZONE[1][k]==1){ 
          printf("tabZONE=[%d %d %d], tabINTER=%d  numZONE=[%d %d]\n",tabZONE[0][k],tabZONE[1][k],tabZONE[2][k],tabINTER[0][k],numZONE[ tabINTER[0][k] ][0],numZONE[ tabINTER[0][k] ][1] );
          printf("  Vdis=%g Edis=%g\n",V,E );
        }
        //
        
        // On met les graines pour le Newton mixte
        if(tabZONE[1][k]==12){
          tabVx[k]=segVa12[ind_mixte];  tabEx[k]=segEa12[ind_mixte];
          tabVy[k]=segVb12[ind_mixte];  tabEy[k]=segEb12[ind_mixte];
        }
        else if(tabZONE[1][k]==13){
          tabVx[k]=segVa13[ind_mixte];  tabEx[k]=segEa13[ind_mixte];
          tabVy[k]=segVc13[ind_mixte];  tabEy[k]=segEc13[ind_mixte];
        }
        else if(tabZONE[1][k]==23){
          tabVx[k]=segVb23[ind_mixte];  tabEx[k]=segEb23[ind_mixte];
          tabVy[k]=segVc23[ind_mixte];  tabEy[k]=segEc23[ind_mixte];
        }
        else if(tabZONE[1][k]==-1){  }
        else{  printf("Erreur dans tabZONE[1]= %d (ecr_benchmark_PT)\n",tabZONE[1][k]); return 1; }
      }
      // ## Plusieurs intersections ##
      else if(nb_inter==2 || nb_inter==3  || nb_inter==4 || nb_inter==5){
        
        //printf("  !!!  zone multi-multi-phase  nb_inter= %d !!!!\n",nb_inter);
        //for(int i=0; i<nb_inter; i++){
        //  printf("  %d : %d\n",i,tabINTER[i][k] );
        //}
 
        nb_zone=0;
        // On cree un tableau de la discretisation des frontières du rectangle
        //printf(" Rectangle :\n");
        //printf(" Vmin= %.10lf  Vmax= %.10lf , Emin= %.10lf, Emax= %.10lf\n",(Vmin+i*hv)*1e6,(Vmin+(i+1)*hv)*1e6,Emin+j*he,Emin+(j+1)*he );
        count=0;
        for(int iv=0; iv<2; iv++){
          for(int ie=0; ie<N_rec; ie++){
            VE_front_rec[count][0]= Vmin + (i+iv)*hv;
            VE_front_rec[count][1]= Emin + j*he + ie*he_sub;
            count++;
          }
        }
        for(int ie=0; ie<2; ie++){
          for(int iv=1; iv<N_rec-1; iv++){
            VE_front_rec[count][0]= Vmin + i*hv + iv*hv_sub;
            VE_front_rec[count][1]= Emin + (j+ie)*he;
            count++;
          }
        }
        // On parcours les frontières du rectangle
        for(int ir=0; ir<4*(N_rec-1); ir++){
          //printf("  ir/(4*N_rec-4)= %d / %d\n", ir, 4*(N_rec-1) );
          
          V_sub=VE_front_rec[ir][0];
          E_sub=VE_front_rec[ir][1];

          if(intPOLY(V_sub,E_sub,NC,SC)){ 
            test=1;
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==3){test=0;}
            }
            if(test==1){ tabZONE[nb_zone][k]=3; nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NA,SA)){
            test=1; 
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==1){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=1;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NB,SB)){
            test=1;            
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==2){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=2;  nb_zone++; }
          }
          else if(ftri(V_sub, E_sub)){
            test=1;
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==0){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=0;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NAB,SAB)){
            test=1;
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==12){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=12;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NAC,SAC)){
            test=1;
            for(int i=0; i<nb_zone; i++){
              if(tabZONE[i][k]==13){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=13;  nb_zone++; }
          }
          else if(intPOLY(V_sub,E_sub,NBC,SBC)){
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==23){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=23;  nb_zone++; }
          }
          else{
            test=1;
            for(int i=0;i<nb_zone;i++){
              if(tabZONE[i][k]==-1){test=0;}
            }
            if(test==1){tabZONE[nb_zone][k]=-1;  nb_zone++; }
          }
        }
        tabNbZONE[k]=nb_zone;
      }
      else{printf("Erreur dans la valeur de nb_inter= %d (ecr_benchmark_PT)\n",nb_inter); return 1;}
      if (tabNbZONE[k]==0){printf("Erreur aucune zone atribuer tabNbZONE[%d]= %d  V=%.10lf , E=%.10lf\n",k,tabNbZONE[k],V,E); return 1; }
      


      //  tri de la fonction enregistrement de zone
      int ind_i, ind_j, temp;
      for(int ind_i=1; ind_i<nb_zone; ++ind_i){        
        while( tabZONE[ind_i-1][k]<tabZONE[ind_i][k] ){
          ind_j=ind_i;
          while( ind_j>0  && tabZONE[ind_j-1][k]<tabZONE[ind_j][k]  ){
            temp=tabZONE[ind_j][k];
            tabZONE[ind_j][k]=tabZONE[ind_j-1][k];
            tabZONE[ind_j-1][k]=temp;
            ind_j--;
          }
        }
      }
      

      k++;
    }
  }




  // fichier binaire   couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  FILE *fileVmy, *fileEmy;
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}

  // Ecriture des tableaux
  fwrite(tabVx,sizeof(double),(Nv-1)*(Ne-1),fileVmx);
  fwrite(tabEx,sizeof(double),(Nv-1)*(Ne-1),fileEmx);
  fwrite(tabVy,sizeof(double),(Nv-1)*(Ne-1),fileVmy);
  fwrite(tabEy,sizeof(double),(Nv-1)*(Ne-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // fichier binaire ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabZONE[0],sizeof(int),(Nv-1)*(Ne-1),fileZONE_0);
  fwrite(tabZONE[1],sizeof(int),(Nv-1)*(Ne-1),fileZONE_1);
  fwrite(tabZONE[2],sizeof(int),(Nv-1)*(Ne-1),fileZONE_2);
  fwrite(tabZONE[3],sizeof(int),(Nv-1)*(Ne-1),fileZONE_3);
  fwrite(tabZONE[4],sizeof(int),(Nv-1)*(Ne-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}

  // Ecriture des tableaux
  fwrite(tabNbZONE,sizeof(int),(Nv-1)*(Ne-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  for(int i=0; i<NA+1; i++){  SA_V[i]=SA[i][0];  SA_E[i]=SA[i][1];  }
  for(int i=0; i<NB+1; i++){  SB_V[i]=SB[i][0];  SB_E[i]=SB[i][1];  }
  for(int i=0; i<NC+1; i++){  SC_V[i]=SC[i][0];  SC_E[i]=SC[i][1];  }

  for(int i=0; i<NAB+1; i++){  SAB_V[i]=SAB[i][0];  SAB_E[i]=SAB[i][1];  }
  for(int i=0; i<NAC+1; i++){  SAC_V[i]=SAC[i][0];  SAC_E[i]=SAC[i][1];  }
  for(int i=0; i<NBC+1; i++){  SBC_V[i]=SBC[i][0];  SBC_E[i]=SBC[i][1];  }

  // Ecriture des tableaux
  fwrite(SA_V,sizeof(double),NA+1,fileA_V);     fwrite(SA_E,sizeof(double),NA+1,fileA_E);
  fwrite(SB_V,sizeof(double),NB+1,fileB_V);     fwrite(SB_E,sizeof(double),NB+1,fileB_E);
  fwrite(SC_V,sizeof(double),NC+1,fileC_V);     fwrite(SC_E,sizeof(double),NC+1,fileC_E);
  fwrite(SAB_V,sizeof(double),NAB+1,fileAB_V);  fwrite(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fwrite(SAC_V,sizeof(double),NAC+1,fileAC_V);  fwrite(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fwrite(SBC_V,sizeof(double),NBC+1,fileBC_V);  fwrite(SBC_E,sizeof(double),NBC+1,fileBC_E);
  

  // debug   //////////////
  /*
  int zone,zone2;
  for(int i=0; i<(Nv-1)*(Ne-1); i++){
    nb_zone=tabNbZONE[i];
    if(nb_zone==2){
      zone=tabZONE[0][i];
      zone2=tabZONE[1][i];

      // On teste d'abord si la premiere zone est mixte
      if(zone==12 || zone==13 || zone==23){ }
      else if(zone==0 || zone==1 || zone==2 || zone==3){ // alors la premiere zone est pure
        if(zone2!=-1){printf("i=%d : Erreur tabZONE=[%d, %d, %d, %d, %d]  |zone =%d, zone2=%d, nb_zone=%d| (lec_benchmark_zone)\n",i,tabZONE[0][i],tabZONE[1][i],tabZONE[2][i],tabZONE[3][i],tabZONE[4][i],zone,zone2,nb_zone );
                      //printf("plot([%g]*1e6 ,[ %g ],'marker','+', 'linewidth',5) \n",V,E );
        }
      }
    }
  }
  */
  //////////



  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);



  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);
  freetab(REC);

  free(tabNbZONE);
  freetab(tabINTER);
  freetab(tabInd_Mixte);
  freetab(numZONE);
  freetab(VE_front_rec);


  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  freetab(SA);   freetab(SB);   freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  freetab(Sa12); freetab(Sb12);
  freetab(Sa13); freetab(Sc13);
  freetab(Sb23); freetab(Sc23);

  freetab(STRI_12);  freetab(STRI_13);  freetab(STRI_23);

  fprintf(fparams, "\nEcriture bien terminé\n");
  fclose(fparams);


  return 0;
}



// Renvoie la zone et les fractoins massiques de chaque phases
/*  Arguments d'entrées
     - Xfm : tableau des 3 fractions massiques  
*/
int lec_benchmark_zone(int affichage, int EOS, double epsilon, int Nmax, double V, double E, 
                       int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                       double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE, 
                       int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax, 
                       double* Xfm, int* pzone, double* pV0x, double *pE0x, double* pV0y, double* pE0y, double *pP, double*pT, int* n_iter, double* pcritere){

  //int Nmax=1e4;
  //double epsilon=1e-10;
  
  // On initialise les valeurs de pression et temperature pour s'assurer qu'on les écrase dans la suite de la fct
  //*pP=-50e9; *pT=-3e3;


  double eps_lambda=1e-13; // tolerance sur la valeur de lambda

  int err;
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;

  if(V<Vmin || Vmax<V){printf("V hors tabulation : Vmin= %g  Vmax= %g  et V= %g (lec_benchmark_zone)\n",Vmin,Vmax,V);   return 1;}  
  if(E<Emin || Emax<E){printf("E hors tabulation : Emin= %g  Emax= %g  et E= %g (lec_benchmark_zone)\n",Emin,Emax,E);   return 1;}

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  // Point triple
  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);


  Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
  
  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);
  int i=(V-Vmin)/hv;
  int j=(E-Emin)/he;


  int pos_fic=i*(Ne-1) + j;

  Vmx=tabVx[pos_fic];  Emx=tabEx[pos_fic];
  Vmy=tabVy[pos_fic];  Emy=tabEy[pos_fic];

  
  // pour le calcul des fractions massiques dans la zone trinagle --> calcul des coord barycentriques dans le triangle
  double* VE=malloc(2*sizeof(double));
  VE[0]=V; VE[1]=E;
  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;
  double* coorBar=malloc(3*sizeof(double));


  // On lit le nombre de zones sur le rectangle dans la tabulation
  int nb_zone=tab_Nb_ZONE[pos_fic];

  int zone, zone2;
  int phaseA, phaseB;

  // On regarde si lorsque qu'il y a 3 zones ou plus, on a -1 dans les zones
  // Dans ce cas on met comme graines segV[N]
  // sinon c'est qu'on est au abord du point triple et dans ce cas on met comme graines les sommets du point triple
  int isout=0;
  if(nb_zone>=3){  for(int i=0; i<nb_zone; i++){  if(tabZONE[i][pos_fic]==-1){isout=1;} } }
  // on se sert du isout dans la boucle ou nb_zone>=3

  //printf("   nb_zone= %d\n",nb_zone );

  //  ##  1 seule zone ##
  if(nb_zone==1){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );

    if(zone==-1){
      if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1;  *pzone=-1; *pP=-100e9; *pT=-1e3; return 0;}
      else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
      else{printf("Erreur dans la valeur affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
    }
    else if(zone==3 || zone==2 || zone==1){  Xfm[zone-1]=1.;  *pzone=zone; *pP=fP(V,E,zone);  *pT=fT(V,E,zone);  return 0;  }
    else if(zone==23 || zone==13 || zone==12){
      phaseB=zone%10;  phaseA=(zone-phaseB)/10;
      err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}
      Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
      *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);  *pT=fT(*pV0x,*pE0x,phaseA);
      return 0;
    }
    else if(zone==0){ 
      err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
      *pzone=0; 
      *pP=PT; *pT=TT; 
      return 0; 
      
    }
  }

  //    TEST D'UN NOUVEAU CODAGE DU REPERAGE DE ZONE DANS LE CAS D'UN RECTANGLE AVEC 2 ZONES
  //    !!!    SACHANT QUE LES ZONES SONT TRIER PAR ORDRE D2CROISSANT DE FACON 0 AVOIR -1 TJRS EN DERNIER  !!!

  //  ##  2 zones   ##   la première zone est forcément la zone pur (ou triple)
  else if(nb_zone==2){
    zone=tabZONE[0][pos_fic];
    zone2=tabZONE[1][pos_fic];
    //printf("  zone=[%d %d]\n",zone,zone2 );

    // On teste d'abord si la premiere zone est mixte
    if(zone==12 || zone==13 || zone==23){ 
      phaseB=zone%10;  phaseA=(zone-phaseB)/10;
      Xfm[0]=0.;  Xfm[1]=0.;  Xfm[2]=0.;
      
      if(zone2==0){ 
        // On est au frontière du point triple donc on teste directement l'appartenance au triangle
        if( ftri(V, E) ){ err=coorbary(VE, TRI, &Xfm); if(err){return 1;}  *pzone=0; *pP=PT; *pT=TT; return 0; }
        else{ // on est dans la zone biphasique --> on lance le Newton
          err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}

          Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
          *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);  *pT=fT(*pV0x,*pE0x,phaseA);

          return 0;
        }
      }
      else if(zone2==-1){ // on frontière externe --> on procede de manière classique : algortihme des piquets
        if(zone==12){
          if( intPOLY(V,E,NAB,SAB) ){
            err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}

            Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
            *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);  *pT=fT(*pV0x,*pE0x,phaseA);
            return 0;
          }
          else{
            if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
            else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
            else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
          }
        }
        else if(zone==13){
          if( intPOLY(V,E,NAC,SAC) ){
            err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}

            Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
            *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);  *pT=fT(*pV0x,*pE0x,phaseA);
            return 0;
          }
          else{
            if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
            else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
            else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
          }
        }
        else if(zone==23){
          if( intPOLY(V,E,NBC,SBC) ){
            err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}

            Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
            *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);  *pT=fT(*pV0x,*pE0x,phaseA);
            return 0;
          }
          else{
            if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1; *pP=-100e9; *pT=-1e3;  return 0;}
            else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
            else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
          }
        }
        else{ printf("Erreur dans les zones de tabulation tabZONE=[%d %d %d ..]\n (lec_benchmark_zone)",zone,zone2,tabZONE[2][pos_fic] ); return 1; }
      }
      else if(zone2==1 || zone2==2 || zone2==3){ // on lance systématiquement le newton --> test sur les fractions massiques
        
        if(EOS==0){ // Hybride + absence de détection des zones
          err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}
          if( x<-eps_lambda || x> (1 + eps_lambda) ){ // on est en zone pure
            Xfm[zone2-1]=1.; *pzone=zone2; *pP=fP(V,E,zone2); *pT=fT(V,E,zone2); return 0; 
          }
          else{ // zone biphasique
            Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
            *pzone=zone;  *pP=fP(*pV0x,*pE0x,phaseA);   *pT=fT(*pV0x,*pE0x,phaseA);
            return 0;
          }
        }
        else if(EOS==1|| EOS==2){ // Hybride + détection de zones
          if(zone==12){
            if(intPOLY(V,E,NAB,SAB)){
              err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}
              Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
              *pzone=zone;
              return 0;
            }
            else{  Xfm[zone2-1]=1.;  *pzone=zone2; *pP=fP(V,E,zone2); *pT=fT(V,E,zone2); return 0;}
          }
          else if(zone==13){
            if(intPOLY(V,E,NAC,SAC)){
              err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}
              Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
              *pzone=zone;
              return 0;
            }
            else{  Xfm[zone2-1]=1.;  *pzone=zone2; *pP=fP(V,E,zone2); *pT=fT(V,E,zone2); return 0;}
          }
          else if(zone==23){
            if(intPOLY(V,E,NBC,SBC)){
              err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);  if(err){return 1;}
              Xfm[phaseA-1]=1.-x;  Xfm[phaseB-1]=x;
              *pzone=zone;
              return 0;
            }
            else{  Xfm[zone2-1]=1.;  *pzone=zone2; *pP=fP(V,E,zone2); *pT=fT(V,E,zone2); return 0;}
          }
        }
        else{printf("Erreur EOS=%d (lec_benchmark_zone)\n",EOS ); return 1;}
      }
      else{ printf("Erreur dans les zones de tabulation tabZONE=[%d %d %d ..]\n (lec_benchmark_zone)",zone,zone2,tabZONE[2][pos_fic] ); return 1; }
    }
    else if(zone==0 || zone==1 || zone==2 || zone==3){ // alors la premiere zone est pure
      if(zone2!=-1){printf("Erreur tabZONE=[%d, %d, %d, %d, %d] tabZONE[1] devrait être égale à -1  |zone =%d, zone2=%d, nb_zone=%d| V=%g E=%g (lec_benchmark_zone)\n",tabZONE[0][pos_fic],tabZONE[1][pos_fic],tabZONE[2][pos_fic],tabZONE[3][pos_fic],tabZONE[4][pos_fic],zone,zone2,nb_zone,V,E );
                    //printf("plot([%g]*1e6 ,[ %g ],'marker','+', 'linewidth',5) \n",V,E );
      }
      if(zone==3){
        if( intPOLY(V,E,NC,SC) ){  Xfm[zone-1]=1.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone);  return 0; }
        else{ // hors domaine
          if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
          else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
          else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;}
        }
      }
      else if(zone==2){
        if( intPOLY(V,E,NB,SB) ){  Xfm[zone-1]=1.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone);  return 0; }
        else{
          if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
          else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
          else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
        }
      }
      else if(zone==1){
        if( intPOLY(V,E,NA,SA) ){  Xfm[zone-1]=1.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone);  return 0; }
        else{
          if(affichage==1){ Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
          else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
          else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
        }
      }
      else if(zone==0){
        if( ftri(V, E)){ err=coorbary(VE, TRI, &Xfm); if(err){return 1;}  *pzone=0; *pP=PT; *pT=TT; return 0; }
        else{
          if(affichage==1){ Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
          else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
          else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
        }
      }
    }
  }

  //  ##  3 zones ou plus  ##
  else{
    for(int ind=0; ind<nb_zone-1; ind++){

      zone=tabZONE[ind][pos_fic];
      //printf("   - zone= %d\n",zone );
      if(zone==3 && intPOLY(V,E,NC,SC)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone); return 0;  }
      else if(zone==2 && intPOLY(V,E,NB,SB)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone); return 0;  }
      else if(zone==1 && intPOLY(V,E,NA,SA)){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone); return 0;  }
      else if(zone==23 && intPOLY(V,E,NBC,SBC)){
        if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[N23-1][1];}
        else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
        err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
        if(err){return 1;}
        Xfm[0]=0.;  Xfm[1]=1.-x;  Xfm[2]=x;   *pzone=23;
        *pP=fP(*pV0x,*pE0x,2);  *pT=fT(*pV0x,*pE0x,2);

        return 0;
      } 
      else if(zone==13 && intPOLY(V,E,NAC,SAC)){
        if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
        else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
        err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
        if(err){return 1;}
        Xfm[0]=1.-x;  Xfm[1]=0;  Xfm[2]=x;   *pzone=13;
        *pP=fP(*pV0x,*pE0x,1);  *pT=fT(*pV0x,*pE0x,1);

        return 0;
      } 
      else if(zone==12 && intPOLY(V,E,NAB,SAB)){
        if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
        else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
        err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
        if(err){return 1;}
        Xfm[0]=1.-x;  Xfm[1]=x;  Xfm[2]=0.;  *pzone=12;
        *pP=fP(*pV0x,*pE0x,1);  *pT=fT(*pV0x,*pE0x,1);

        return 0;
      }
      else if(zone==0 && ftri(V, E)){  err=coorbary(VE, TRI, &Xfm); if(err){return 1;}   *pzone=0; *pP=PT; *pT=TT;  return 0;  }  
    }
    
    zone=tabZONE[nb_zone-1][pos_fic];
    if(zone==3 || zone==2 || zone==1){  Xfm[zone-1]=1.;  Xfm[(zone-1+1)%3]=0.;  Xfm[(zone-1+2)%3]=0.;  *pzone=zone; *pP=fP(V,E,zone); *pT=fT(V,E,zone); return 0;  }
    else if(zone==23){
      if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
      else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
      if(err){return 1;}
      Xfm[0]=0.;  Xfm[1]=1.0-x;  Xfm[2]=x;  *pzone=23;
      *pP=fP(*pV0x,*pE0x,2);  *pT=fT(*pV0x,*pE0x,2);

      return 0;
    }
    else if(zone==13){
      if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
      if(err){return 1;}
      Xfm[0]=1.-x;  Xfm[1]=0.;  Xfm[2]=x;  *pzone=13;
      *pP=fP(*pV0x,*pE0x,1);  *pT=fT(*pV0x,*pE0x,1);

      return 0;
    }  
    else if(zone==12){
      if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, n_iter, pcritere);
      if(err){return 1;}
      Xfm[0]=1.-x;  Xfm[1]=x;  Xfm[2]=0.;   *pzone=12;
      *pP=fP(*pV0x,*pE0x,1);  *pT=fT(*pV0x,*pE0x,1);

      return 0;
    }
    else if(zone==0){  err=coorbary(VE, TRI, &Xfm); if(err){return 1;}   *pzone=0; *pP=PT; *pT=TT; return 0;  }
    else if(zone==-1){  
      if(affichage==1){  Xfm[0]=-1; Xfm[1]=-1; Xfm[2]=-1; *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %g , E= %g (lec_benchmark_zone)\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d (lec_benchmark_zone)\n",affichage); return 1;} 
    }
    else{printf("Erreur dans tabZONE[%d][%d]= %d (lec_benchmark_zone)\n\n",nb_zone-1,pos_fic,tabZONE[nb_zone-1][pos_fic]); return 1;}
  }
  // Si on arrive la c'est que la bonne zone n'a pas été trouvé
  printf("Erreur choix zone (lec_benchmark_zone)\n");
  printf("  nb_zone= %d",nb_zone );
  printf("  tabZONE[0][%d]= %d",pos_fic,tabZONE[0][pos_fic]);
  for(int ind=1; ind<nb_zone; ind++){
    printf("  tabZONE[%d][%d]= %d",ind,pos_fic,tabZONE[ind][pos_fic]); 
  }
  printf("\n V= %.10lf , E= %.10lf\n",V,E);
  return 1; 
}




// Script pour tester le Benchmark
int use_benchmark_PT(int EOS, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err;
  double E,V;
  double P,T;
  double x;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  double* Xfm=malloc(3*sizeof(double));

  int n_iter;
  double critere;


  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  /* pour pouvoir utiliser la discretisations juste en lisant les fichiers
  Organisation du fichier :
    nb_points
    Nv_dis , Ne_dis
    Vmin_dis , Vmax_dis
    Emin_dis , Emax_dis
  */
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3;
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC);

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixtes
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  
  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);


  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);

 
  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);


 //////////////////////////////////////////////////////////////////////////////


  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXtest, *fCtest, *fGtest;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXtest = fopen("fichiers/Xtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier C test\n");return 1;}
  if((fGtest = fopen("fichiers/Gtest.txt", "w+")) == NULL){printf("erreur ouverture fichier G test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
    E=Emin_test+j*he_test;
    fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
    V=Vmin_test+j*hv_test;
    fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);


  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.15lf \n",Vmin_dis+(i+0.5)*hv_dis );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.15lf \n",Emin_dis+(i+0.5)*he_dis );
  }
  fclose(file_E_dis);
  //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  
  double G;
  int Nmax=20;
  double epsilon=1e-8;

  double hv=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he=(Emax_dis-Emin_dis)/(Ne_dis-1);
  int i_tab,j_tab;
  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      E=Emin_test+j*he_test;

      //printf("V=%g; E=%g;\n",V,E );

      // Tabulation pure (ne produit que l'effet mais n'est pas optimisé)
      // on se contente de modifire le couple (V,E) d'entreé pour qu'il coresponde au centre de l'élement de tabulation 
      // puis calcul les granduers thermo a precision machine
      if(EOS==2){
        i_tab=(V-Vmin_dis)/hv;
        j_tab=(E-Emin_dis)/he; 
        V = Vmin_dis + (i_tab+0.5)*hv;
        E = Emin_dis + (j_tab+0.5)*he;
        /*
        printf("tab pure : Vnew=%g Enew=%g\n",V,E );
        printf("rectangle :\n");
        printf("Vl=%g  Vr=%g\n", Vmin_dis + (i_tab)*hv, Vmin_dis + (i_tab+1)*hv);
        printf("El=%g  Er=%g\n", Emin_dis + (j_tab)*he, Emin_dis + (j_tab+1)*he);
        */
      }
      else{
        // rien on garde le couple V,E d'entreé
      }

      // Fonction renvoie la zone, la fraction massique et les graines pour les Newton_mixte
      err= lec_benchmark_zone(1, EOS, epsilon, Nmax, V, E, nb_points, 
                              SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                              Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, 
                              Xfm, &zone, &V0x, &E0x, &V0y, &E0y, &P, &T, &n_iter, &critere);
      if(err){printf("lec_benchmark_zone depuis use_benchmark_PT\n");  return err;}
      
      
      if(zone==-1){c=-1e3; G=-1e1;}
      else{
        err=c_celerite_son(Nmax, epsilon,zone, P, V, E, dV, V0x, E0x, V0y, E0y, &c, &G);
        if(err){return err;}
      }


      fprintf(fPtest,"%.5g\n",P);
      fprintf(fTtest,"%.5g\n",T);
      fprintf(fXtest,"%.5g %.5g %.5g\n",Xfm[0],Xfm[1],Xfm[2]);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.5g\n",c);      
      fprintf(fGtest,"%.5g\n",G); // dérivée fondamentale      
    }
  }
  fclose(fPtest);  fclose(fTtest);  fclose(fXtest); fclose(fCtest); fclose(fGtest);
  

  // Liberation de la mémoire
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}



//  Remplie les tableaux suivants pour les utiliser dans le main 
/*
*/
int lec_tabulation(int nb_points, int Nv_dis, int Ne_dis, int** tabZONE, int* tab_Nb_ZONE, double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC){

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);
  
  //printf(" ok 3\n");
  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);
  
  //printf(" ok 4\n");
  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  

  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);

  
  //printf(" ok fin lec\n");
  return 0;
}

// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec tabulation + calcul #############################



// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec calcul uniquement (sans tabulation) #############################

int graine_zone_mixte_modif(double V, double E, int N, double** Sp, double* pVmx, double* pEmx, double* pVmy, double* pEmy){
  int res;

  int Ntot=2*N-1;
  
  // reparage du quadrilatere
  double** S=alloctabd(5,2);
  int ind=0;
  int testNO=1;
  for(int i=0; i<=N-2; i++){
    S[0][0]=Sp[Ntot-i][0];     S[0][1]=Sp[Ntot-i][1];
    S[1][0]=Sp[Ntot-(i+1)][0]; S[1][1]=Sp[Ntot-(i+1)][1];
    S[2][0]=Sp[i+1][0];        S[2][1]=Sp[i+1][1];
    S[3][0]=Sp[i][0];          S[3][1]=Sp[i][1];
    S[4][0]=Sp[Ntot-i][0];     S[4][1]=Sp[Ntot-i][1];

    /*
    if( fabs(Sp[i][0] - Sp[i+1][0])>fabs(Sp[i][0] - V) || fabs(Sp[N+i][0] - Sp[N+i+1][0])>fabs(Sp[i][0] - V) || fabs(Sp[i][1] - Sp[i+1][1])>fabs(Sp[i][1] - E) || fabs(Sp[N+i][1] - Sp[N+i+1][1])>fabs(Sp[i][1] - E)   ){  
      printf(" V= %g; E= %g; \n", V,E);
      printf(" Sp[N+%d][0]= %g;  Sp[N+%d][1]= %g; \n",i,Sp[N+i][0],i,Sp[N+i][1]);
      printf(" Sp[N+%d+1][0]= %g;  Sp[N+%d+1][1]= %g; \n",i,Sp[N+i+1][0],i,Sp[N+i+1][1]);
      printf(" Sp[%d+1][0]= %g;  Sp[%d+1][1]= %g; \n",i,Sp[i+1][0],i,Sp[i+1][1]);
      printf(" Sp[%d][0]= %g;  Sp[%d][1]= %g; \n",i,Sp[i][0],i,Sp[i][1]);
    }
    */

    if(intPOLY(V, E, 4, S)){ind=i; testNO=0;   break;  } // ou i=N
  
  }

  if(testNO){
    printf("\n");
    printf("le point est dans aucun des quadrilateres (graine_zone_mixte)\n");
    printf(" ind= %d\n",ind );
    printf("V= %g , E= %g \n",V,E); 
    return 1;
  }

  *pVmx=Sp[N+ind][0]; *pEmx=Sp[N+ind][1];
  *pVmy=Sp[ind][0];   *pEmy=Sp[ind][1];

  
  free(S);
}


// #######################  FONCTIONS POUR LE CALCUL DE P,T et C avec tabulation uniquement #############################

int Newton_E_init(double* A, double* InvA, double* Delta, double* dX, int phaseX, int phaseY, double P, double T, double Vx, double Ex, double Vy, double Ey, double* pVX, double* pEX, double* pVY, double* pEY){
  int dim=4;

  // Grandeurs thermodynamiques
  double Px, Tx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  
  
  //printf("Px= %g, Tx= %g\n",Px,Tx );
  //printf("Py= %g, Ty= %g\n",Py,Ty );
  
  
  Delta[0]=P-Px;   Delta[1]=P-Py;  
  Delta[0]=T-Tx;   Delta[1]=T-Ty;  
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;  A[dim*0+1]=dPex;  A[dim*0+2]=0.;    A[dim*0+3]=0.; 
  
  A[dim*1+0]=0.;    A[dim*1+1]=0.;    A[dim*1+2]=dPvy;  A[dim*1+3]=dPey;

  A[dim*2+0]=dTvx;  A[dim*2+1]=dTex;  A[dim*2+2]=0.;    A[dim*2+3]=0.;     
  
  A[dim*3+0]=0.;    A[dim*3+1]=0.;    A[dim*3+2]=dTvy;  A[dim*3+3]=dTey;    
  
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
  //printf("detA= %g\n",det_a );
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
  
  
  // Mutiplcation matricielle de A et InvA  TEST
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
  


  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */

  return 0;
}


int E_mixte_init(int Nmax, double epsilon, int phaseX, int phaseY, double* lambda, double P, double T, double V0x, double E0x, double V0y, double E0y,  double* pVx, double* pEx, double* pVy, double* pEy){
  int n=1;
  int dim=4;
  int err;
   
  // printf(" epsilon= %g\n",epsilon);
  /*
  printf("  phaseX= %d , phaseY= %d\n",phaseX,phaseY );
  printf("  Vmel= %g , Emel= %g\n",Vmel, Emel );
  */


  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx;
  double Py, Ty;

  double Vx, Ex, Vy, Ey;
  double VX, EX, VY, EY;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    err=Newton_E_init(A, InvA, Delta, dX, phaseX, phaseY, P, T, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY);
    if(err){return 1;}

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;

    //critere=fabs(Gx-Gy)+fabs(Px-Py)+fabs(Ty-Tx);
    critere=fabs(dX[0]/Vx)+fabs(dX[1]/Ex)+fabs(dX[2]/Vy)+fabs(dX[3]/Ey);
    //critere=fabs(Delta[0])+fabs(Delta[1])+fabs(Delta[2])+fabs(Delta[3]);
    //printf("  n= %d    critere= %g\n",n,critere );
    n++;
  }
  //printf(" --n= %d, critere= %g (VE_mixte_cin_phase)\n",n,critere );

  if(n==Nmax){printf("Newton_mixte_cin_phase non convergent : -n= %d , critere= %g (E_mixte_init)\n",n,critere ); return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
  
  return 0;
}


/// Petite fonction pour le schema Godunov en energie totale
int epsilon_etain(double V, double P, double T, double* pE, double* lambda, int nb_points,
                  double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC ){
  
  int err, Nmax=3e2;
  double epsilon=1e-9;

  double Vf,Ef;
  int zone, phaseX, phaseY;
  // nb_points est l'argument d'entree de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  
  double Vx, Ex, Vy, Ey;
  double V0x, E0x, V0y, E0y;
  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0;
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  // Valeur d'initialisation
  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);

  //printf(" phase= %d\n",phase );
  if(lambda[0]==1 && lambda[1]==0 && lambda[2]==0){zone=1; phaseX=1; Vx=V;}
  else if(lambda[0]==0 && lambda[1]==1 && lambda[2]==0){zone=2; phaseX=2; Vx=V; }
  else if(lambda[0]==0 && lambda[1]==0 && lambda[2]==1){zone=3; phaseX=3; Vx=V;}
  
  else if(lambda[0]>0 && lambda[1]>0 && lambda[2]==0){
    zone=12; phaseX=1; phaseY=2;
    // VALEUR INITIALE
    V0x=Va0; E0x=Ea0; V0y=Vb0; E0y=Eb0;
    err=E_mixte_init(Nmax, epsilon, phaseX, phaseY, lambda, P, T, V0x, E0x, V0y, E0y,  &Vx, &Ex, &Vy, &Ey);
    if(err){return 1;}
  }
  else if(lambda[0]>0 && lambda[1]==0 && lambda[2]>0){
    zone=13; phaseX=1; phaseY=3;
    // VALEUR INITIALE
    V0x=Va0; E0x=Ea0; V0y=Vc0; E0y=Ec0;
    err=E_mixte_init(Nmax, epsilon, phaseX, phaseY, lambda, P, T, V0x, E0x, V0y, E0y,  &Vx, &Ex, &Vy, &Ey);
    if(err){return 1;}
  }
  else if(lambda[0]==0 && lambda[1]>0 && lambda[2]>0){
    zone=23; phaseX=2; phaseY=3;
    // VALEUR INITIALE
    V0x=Vb0; E0x=Eb0; V0y=Vc0; E0y=Ec0;
    err=E_mixte_init(Nmax, epsilon, phaseX, phaseY, lambda, P, T, V0x, E0x, V0y, E0y,  &Vx, &Ex, &Vy, &Ey);
    if(err){return 1;}
  }
  else if(lambda[0]>0 && lambda[1]>0 && lambda[2]>0){
    zone=0; 
    *pE=lambda[0]*EaT + lambda[1]*EbT +lambda[2]*EcT;
    return 0;
  } // on connait E0a, etc..
  

  *pE=fE_VP(Vx,P,phaseX);
  return 0;
}



// ********                                                                  ************
///////////*************  CINEMATIQUE DE CHANGEMENT DE PHASE *************////////////////
// ********                                                                  ************


// *******  Newton pour la cinétique de chgmt de phase  ********

/*  Valeur d'entree :
      Vmel, Emel : 
      lambda : fraction massique
      Va,Ea, Vb, Eb : Vraiables


   Valeur sortie
     Vx,Ex, Vy,Ey, 

   double* A et InvA permet une seule allocation 

*/
void Newton_cin_phase(double* A, double* InvA, double* Delta, int phaseX, int phaseY, double Vmel, double Emel, double lambda, double Vx, double Ex, double Vy, double Ey, double* pVX, double* pEX, double* pVY, double* pEY, double* dX){
  int dim=4;

  // Grandeurs thermodynamiques
  double Px, Tx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, dPvy, dPey, dTvy, dTey;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  
  Delta[0]=Py-Px;   Delta[1]=Ty-Tx;  
  Delta[2]=Vmel - (1.-lambda)*Vx - lambda*Vy;
  Delta[3]=Emel - (1.-lambda)*Ex - lambda*Ey;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;       A[dim*0+1]=dPex;       A[dim*0+2]=-dPvy;   A[dim*0+3]=-dPey; 
  
  A[dim*1+0]=dTvx;       A[dim*1+1]=dTex;       A[dim*1+2]=-dTvy;   A[dim*1+3]=-dTey;
  
  A[dim*2+0]=1.-lambda;  A[dim*2+1]=0.;         A[dim*2+2]=lambda;  A[dim*2+3]=0.;     
  
  A[dim*3+0]=0.;         A[dim*3+1]=1.-lambda;  A[dim*3+2]=0.;      A[dim*3+3]=lambda;    
  
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
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-8){printf("A non inversible : detA= %g\n", det_a);}
  
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
  
  
  // Mutiplcation matricielle de A et InvA  TEST
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
  


  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];

  /*
  printf("  Vx= %.10lf, Ex= %.10lf   Vy= %.10lf, Ey= %.10lf    x= %.10lf\n",Vx,Ex,Vy,Ey,x );
  printf("  (Vmel-Vx)/(Vy-Vx)= %.10lf  ,    (Emel-Ex)/(Ey-Ex)= %.10lf\n", (Vmel-Vx)/(Vy-Vx), (Emel-Ex)/(Ey-Ex));
  printf("  (1-x)*Vx + x*Vy - Vmel= %.10lf  ,    (1-x)*Ex + x*Ey - Emel= %.10lf\n", (1-x)*Vx + x*Vy - Vmel, (1-x)*Ex + x*Ey - Emel);
  */

}

/*
Arguments d'entree :
  - Nmax : nombre max d'iteration
  - epsilon : tolérance sur le critère  (ici fabs(G1-G2))
  - Vmel, Emel : Point dont on sait qu'il est dans la zone de mélange (alpha/beta) (X/Y)
  - phaseX, phaseY : les deux phases dans la zone de mélange

Arguments de sortie :
  - Vx, Ex : valeur de la phase X
  - Vy, Ey : valeur de la phase Y
  -   x    : fraction massique de Y
*/
int VE_mixte_cin_phase(int Nmax, double epsilon, int phaseX, int phaseY, double lambda, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy){
  int n=1;
  int dim=4;
   
  // printf(" epsilon= %g\n",epsilon);
  /*
  printf("  phaseX= %d , phaseY= %d\n",phaseX,phaseY );
  printf("  Vmel= %g , Emel= %g\n",Vmel, Emel );
  */

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx;
  double Py, Ty;

  double Vx, Ex, Vy, Ey;
  double VX, EX, VY, EY;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 

  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);


  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    Newton_cin_phase(A, InvA, Delta, phaseX, phaseY, Vmel, Emel, lambda, Vx, Ex, Vy, Ey, &VX, &EX, &VY, &EY, dX);

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;

    critere=fabs(dX[0]/Vx)+fabs(dX[1]/Ex)+fabs(dX[2]/Vy)+fabs(dX[3]/Ey);
    
    n++;
  }
  //printf(" --n= %d, critere= %g (VE_mixte_cin_phase)\n",n,critere );

  if(n==Nmax){printf("Newton_mixte_cin_phase non convergent : -n= %d , critere= %g (VE_mixte_cin_phase)\n",n,critere ); return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;

  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
  
  return 0;
}



int Newton_cin_3_phase_all_dom(double* A, double* InvA, double* Delta, double Vmel, double Emel, double* lambda, double Vx, double Ex, double Vy, double Ey, double Vz, double Ez, double* pVX, double* pEX, double* pVY, double* pEY, double* pVZ, double* pEZ, double* dX){
  int dim=6;

  double lambdaX,lambdaY,lambdaZ;
  lambdaX=lambda[0];
  lambdaY=lambda[1];
  lambdaZ=lambda[2];

  int phaseX=1, phaseY=2, phaseZ=3;

  // Grandeurs thermodynamiques
  double Px, Tx, dPvx, dPex, dTvx, dTex;
  double Py, Ty, dPvy, dPey, dTvy, dTey;
  double Pz, Tz, dPvz, dPez, dTvz, dTez;
  // phase X
  Px=fP(Vx,Ex, phaseX); Tx=fT(Vx,Ex, phaseX);
  dPvx=dPdV(Vx,phaseX); dPex=dPdE(phaseX);
  dTvx=dTdV(Vx,phaseX); dTex=dTdE(phaseX);
  // phase Y
  Py=fP(Vy,Ey, phaseY); Ty=fT(Vy,Ey, phaseY);
  dPvy=dPdV(Vy,phaseY); dPey=dPdE(phaseY);
  dTvy=dTdV(Vy,phaseY); dTey=dTdE(phaseY);
  // phase Z
  Pz=fP(Vz,Ez, phaseZ); Tz=fT(Vz,Ez, phaseZ);
  dPvz=dPdV(Vz,phaseZ); dPez=dPdE(phaseZ);
  dTvz=dTdV(Vz,phaseZ); dTez=dTdE(phaseZ);
  
  Delta[0]=Py-Px;   Delta[1]=Pz-Px;
  Delta[2]=Ty-Tx;   Delta[3]=Tz-Tx;
  Delta[4]=Vmel - lambdaX*Vx - lambdaY*Vy - lambdaZ*Vz;
  Delta[5]=Emel - lambdaX*Ex - lambdaY*Ey - lambdaZ*Ez;
   
  /*
  for(int j=0; j<dim; j++){
    printf("Delta[%d]= %g  ",j, Delta[j] );
  }
  printf("\n");
  */
  
  
  // on met le second menbre dans la dernière colonne
  //matrice A
  A[dim*0+0]=dPvx;     A[dim*0+1]=dPex;     A[dim*0+2]=-dPvy;    A[dim*0+3]=-dPey;    A[dim*0+4]=0.;       A[dim*0+5]=0.;  

  A[dim*1+0]=dPvx;     A[dim*1+1]=dPex;     A[dim*1+2]=0.;       A[dim*1+3]=0.;       A[dim*1+4]=-dPvz;    A[dim*1+5]=-dPez; 

  A[dim*2+0]=dTvx;     A[dim*2+1]=dTex;     A[dim*2+2]=-dTvy;    A[dim*2+3]=-dTey;    A[dim*2+4]=0.;       A[dim*2+5]=0.; 

  A[dim*3+0]=dTvx;     A[dim*3+1]=dTex;     A[dim*3+2]=0.;       A[dim*3+3]=0.;       A[dim*3+4]=-dTvz;    A[dim*3+5]=-dTez; 

  A[dim*4+0]=lambdaX;  A[dim*4+1]=0.;       A[dim*4+2]=lambdaY;  A[dim*4+3]=0.;       A[dim*4+4]=lambdaZ;  A[dim*4+5]=0.; 

  A[dim*5+0]=0.;       A[dim*5+1]=lambdaX;  A[dim*5+2]=0.;       A[dim*5+3]=lambdaY;  A[dim*5+4]=0.;       A[dim*5+5]=lambdaZ;
  
 
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
  //printf("detA= %g\n",det_a );
  if ( fabs(det_a)<1e-9){printf("A non inversible : detA= %g , V=%g; E=%g;\n", det_a,Vmel, Emel); return 1;}
  
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
  
  // Mutiplcation matricielle de A et InvA  TEST
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
    //printf("  dX[%d]= %g\n",i,dX[i] );
  }
  //printf("\n");
  

  double correc=1.;

  *pVX=Vx+correc*dX[0];
  *pEX=Ex+correc*dX[1];
  *pVY=Vy+correc*dX[2];
  *pEY=Ey+correc*dX[3];
  *pVZ=Vz+correc*dX[4];
  *pEZ=Ez+correc*dX[5];


  return 0;
}



int VE_mixte_cin_3_phase_all_dom(int Nmax, double epsilon, double* lambda, double V0x, double E0x, double V0y, double E0y, double V0z, double E0z, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy, double* pVz, double* pEz){
  int err;
  int n=1;
  int dim=6;
  double errV,errE;

  double* A = (double *) malloc(sizeof(double) * (dim*dim));
  double* InvA = (double *) malloc(sizeof(double) * (dim*dim));
  double* Delta=malloc(dim*sizeof(double));
  double* dX=malloc(dim*sizeof(double));
  
  double Px, Tx;
  double Py, Ty;
  double Pz, Tz;

  double Vx, Ex, Vy, Ey, Vz, Ez;
  double VX, EX, VY, EY, VZ, EZ;

  Vx=V0x; Ex=E0x; 
  Vy=V0y; Ey=E0y; 
  Vz=V0z; Ez=E0z; 

  // phase X
  Px=fP(Vx,Ex, 1); Tx=fT(Vx,Ex, 1);
  // phase Y
  Py=fP(Vy,Ey, 2); Ty=fT(Vy,Ey, 2);
  // phase Z
  Pz=fP(Vz,Ez, 3); Tz=fT(Vz,Ez, 3);

  double critere;
  
  critere=1.;
  //printf("  critere= %.14lf\n",critere );
    
  n=0;
  while(critere>=epsilon && n<Nmax){
    
    err= Newton_cin_3_phase_all_dom(A, InvA, Delta, Vmel, Emel, lambda, Vx, Ex, Vy, Ey, Vz, Ez, &VX, &EX, &VY, &EY, &VZ, &EZ, dX);
    if(err){return 1;}

    Vx=VX; Ex=EX;
    Vy=VY; Ey=EY;
    Vz=VZ; Ez=EZ;

    errV=(Vmel - lambda[0]*Vx - lambda[1]*Vy - lambda[2]*Vz)/Vmel;
    errE=(Emel - lambda[0]*Ex - lambda[1]*Ey - lambda[2]*Ez)/Emel;

    critere=fabs(dX[0]/Vx)+fabs(dX[1]/Ex)+fabs(dX[2]/Vy)+fabs(dX[3]/Ey)+fabs(dX[4]/Vz)+fabs(dX[5]/Ez);
    //critere=fmax(fabs(errV),fabs(errE));
    
    n++;
  }
  //printf(" --n= %d, critere= %g (VE_mixte_cin_3_phase_all_dom)\n",n,critere );

  if(n==Nmax){printf("Newton_mixte_cin_phase non convergent : -n= %d (Nmax=%d), critere= %g  V=%g E=%g (VE_mixte_cin_3_phase_all_dom)\n",n,Nmax,critere,Vmel,Emel ); return 1;}
  if(isnan(critere)){printf("Newton_mixte_cin_phase non convergent : -n= %d , critere= %g  V=%g E=%g (VE_mixte_cin_3_phase_all_dom)\n",n,critere,Vmel,Emel ); return 1;}

  *pVx=Vx; *pEx=Ex;
  *pVy=Vy; *pEy=Ey;
  *pVz=Vz; *pEz=Ez;

  //printf("  ok fin boucle i\n");
  free(A); 
  free(InvA);
  free(Delta);
  free(dX);
  
  return 0;
}



// Calcul de la célérite du son
/*  Pp=P(V+dV,E+dE(+dV) )  
    Pm=P(V-dV,E+dE(-dV) )  

    Argument de sortie :kugh
     - double* pc : pointeur sur la celerité du son
*/
int c_celerite_son_cin_phase(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, 
                             double V0x, double E0x, double V0y, double E0y, double* lambda, double *pc){

  int err;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,x;
  double lambdaB;

  int phaseA,phaseB;


  if(zone==-1){printf("Erreur out of domain : zone= %d (c_celerite_son_cin_phase)\n",zone); return 1;  }
  else if(zone==1 || zone==2 || zone==3){  
    Pp=fP(V+dV,E+dE,zone); 
    Pm=fP(V-dV,E-dE,zone); 
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseB=zone%10;
    phaseA=(zone-phaseB)/10;
    lambdaB=lambda[phaseB-1];

    err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
    //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pp=fP(Vx,Ex,phaseA);

    err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
    //err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x);  if(err){return 1;}
    Pm=fP(Vx,Ex,phaseA);
  }
  else if(zone==0){
    *pc=0;
    return 0;
  }

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 

  
  if(c2<0){
    printf("Erreur c2= %g negatif \n",c2);
    printf("  V= %g , E= %g , zone= %d \n",V,E,zone);  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaB= %g , c^2= %g \n",lambdaB, c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}


int c_celerite_son_cin_3_phase_all_dom(int Nmax, double epsilon, double P, double V, double E, double dV, double* lambda, double *pc){

  int err;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;


  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pp=fP(Vx,Ex,1);

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pm=fP(Vx,Ex,1);

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  //double dPdV=(dPpdV + dPmdV)/2;
  double dPdV=(Pp - Pm)/(2*dV);
  double c2=-V*V*dPdV; 

  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g \n",V,E);  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}

int c_G_3_phase_all_dom(int Nmax, double epsilon, double Va0, double Ea0, double Vb0, double Eb0, double Vc0, double Ec0, double P, double V, double E, double dV, double* lambda, double *pc, double *pG){
  int err;

  double Pb,Pa;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  //Pb=(fP(Vx,Ex,1)+fP(Vx,Ex,2)+fP(Vx,Ex,3))/3;
  Pb=fP(Vx,Ex,1);
 
  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  //Pa=(fP(Vx,Ex,1)+fP(Vx,Ex,2)+fP(Vx,Ex,3))/3;
  Pa=fP(Vx,Ex,1);

  double dPdV_S=(Pb-Pa)/(2*dV);
  double d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);
  
  double c2=-V*V*dPdV_S;
  *pc=sqrt(c2);
  *pG=(-V/2)*d2PdV2_S/dPdV_S;
  
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g ",V,E);  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pb,Pa );
    printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
    return 1;
  }
  return 0;
}

/*  Calcul de la célrité du son à l'équilibre thermo (par taux de variation)

*/
int c_G_eq_thermo(int Nmax, double epsilon, int zone, double P, double V, double E, double dV, double* lambda, double *pc, double *pG){
  int err;
  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;
  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0;
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  if (zone==1 || zone==2 || zone==3){ // zone monophasique
    Pp=fP(V+dV,E+dE,zone);
    Pm=fP(V-dV,E-dE,zone);
  }
  else if (zone==12 || zone==13 || zone==23){ // zone biphasique
    err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
    //err= VE_mixte_cin_phase(Nmax, epsilon, int phaseX, int phaseY, double lambda, double V0x, double E0x, double V0y, double E0y, double Vmel, double Emel, double* pVx, double* pEx, double* pVy, double* pEy){

    if(err){printf("erreur c_celerite_son_init\n"); return 1;}
    Pp=fP(Vx,Ex,1);

    err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
    if(err){printf("erreur c_celerite_son_init\n"); return 1;}
    Pm=fP(Vx,Ex,1);
  }
  else if(zone==0){ // zone triphasique
    
  }

  

  
  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 


  return 0;
}

int c_celerite_son_init(double P, double V, double E, double* lambda, double *pc){
  int err;

  int Nmax=1e3;
  double epsilon=1e-5;
  
  double dV=V*1e-4;
  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;
  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0;
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur c_celerite_son_init\n"); return 1;}
  Pp=fP(Vx,Ex,1);

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur c_celerite_son_init\n"); return 1;}
  Pm=fP(Vx,Ex,1);

  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 



  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g \n",V,E);  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}

int c_G_eq(int Nmax, double epsilon, double V, double E, int zone, double P, double T, 
                      double V0x, double E0x, double V0y, double E0y, double* pc, double* pG){
  int err;
  double critere=0;
  int n_iter=0;

  double Pb, Pa;

  double dPdV_S, d2PdV2_S;

  double Pe,Pv, Pvv;

  double dV=1e-8;
  double dE=-P*dV;

  double Vx, Ex, Vy, Ey;
  int phaseA, phaseB;
  double x=0;

  // Calcul des dérivées partielles
  if(zone==1||zone==2||zone==3){
    
    Pe=dPdE(zone);
    Pv=dPdV(V, zone);
    Pvv=d2PdV2(V, zone);

    dPdV_S=Pv-P*Pe;
    d2PdV2_S=Pvv - Pv*Pe + P*Pe*Pe;
    
    /*
    Pb=fP(V+dV,E+dE,zone);
    Pa=fP(V-dV,E-dE,zone);
    dPdV_S=(Pb-Pa)/(2*dV);
    d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);
    */

    *pc=sqrt(-V*V*dPdV_S);
    *pG=(-V/2)*d2PdV2_S/dPdV_S;
  }
  else if(zone==12||zone==13||zone==23){
    phaseB=zone%10;  phaseA=(zone-phaseB)/10;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);
    if(err){printf("(c_G_eq)\n"); return err;}
    Pb=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);
    if(err){printf("(c_G_eq) \n"); return err;}
    Pa=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;

    dPdV_S=(Pb-Pa)/(2*dV);
    d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);

    *pc=sqrt(-V*V*dPdV_S);
    *pG=(-V/2)*d2PdV2_S/dPdV_S;
  }
  else{printf(" Erreur zone=%d (c_G_eq)\n",zone ); return 1;}

  return 0;
}


int G_S_cin_3_phase(double epsilon, double P, double V, double E, double* lambda, double *pc, double* pG, double* pS){

  int err;
  int Nmax=50;
  //double epsilon=1e-8;

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  double dV=V*1e-6;

  double Pp,Pm;
  double dE=-P*dV;

  double Vx,Ex,Vy,Ey,Vz,Ez;

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pp=(fP(Vx,Ex,1)+fP(Vy,Ey,2)+fP(Vz,Ez,3))/3;

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  Pm=(fP(Vx,Ex,1)+fP(Vy,Ey,2)+fP(Vz,Ez,3))/3;

  // Pour le calcul des (V,E) de chaque phase pour le calcul d'entropie S
  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V, E, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
  //P=fP(Vx,Ex,1);


  // Vitesse du son
  double dPpdV=(Pp-P)/dV;
  double dPmdV=(Pm-P)/(-dV);
  double dPdV=(dPpdV + dPmdV)/2;
  double c2=-V*V*dPdV; 

  // Derivée fondamentale
  double d2PdV2=(Pm+Pp-2*P)/(dV*dV);
  *pG=-(V/2)*d2PdV2/dPdV;
  
  /*
  printf("Erreur c2= %g \n",c2);
  printf("  V= %g , E= %g \n",V,E);  
  printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
  printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
  printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
  printf("  d2P/dV2= %g \n",d2PdV2 );
  printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
  printf("  *pG=-(V/2)*d2PdV2/dPdV \n");
  printf("  %g =-(%g/2)*%g/%g \n", *pG,V,d2PdV2,dPdV);
  */
  
 
  // Entropie 
  double Sa=fS(Vx, Ex, 1);
  double Sb=fS(Vy, Ey, 2);
  double Sc=fS(Vz, Ez, 3);
  *pS = lambda[0]*Sa + lambda[1]*Sb + lambda[2]*Sc;
   
  if(c2<0){
    printf("Erreur c2= %g negatif\n",c2);
    printf("  V= %g , E= %g \n",V,E);  
    printf("  dV/V= %g , dE/E= %g \n",dV/V,dE/E );
    printf("  P= %g , Pp= %g , Pm= %g\n",P,Pp,Pm );
    printf("  dPp/dV= %g , dPm/dV= %g , dPdV= %g \n",dPpdV,dPmdV,dPdV );
    printf("  lambdaA= %g , lambdaB= %g , lambdaC= %g , c^2= %g \n",lambda[0],lambda[1],lambda[2], c2 );
    //*pc=3.3e2;
    return 1;
  }
  *pc=sqrt(c2);
  return 0;
}


//  Implemenation du schema temporelle des fraction massiques
int lambda_cin_phase(int sch_cin_phase, double* lambda_out, double* lambda_in, double* lambda_eq, double dt, double Ttau ){
  int sgn;
  if(sch_cin_phase==0){ // Schéma implcite d'odre 1 classique,
    lambda_out[0]=((dt/Ttau)*lambda_eq[0]+lambda_in[0])/((dt/Ttau)+1.);
    lambda_out[1]=((dt/Ttau)*lambda_eq[1]+lambda_in[1])/((dt/Ttau)+1.);
    lambda_out[2]=((dt/Ttau)*lambda_eq[2]+lambda_in[2])/((dt/Ttau)+1.);
  }
  //else if(sch==1){ // Schéma implcite d'odre 1 Fickett-Rivard

  return 0;
}


int test_consis_thermo_eq(int Nmax, double epsilon, double V, double E, int zone, 
                          double V0x, double E0x, double V0y, double E0y, double* pres){
  int err;
  double critere=0;
  int n_iter=0;

  double P,T;

  double Pp_e,Pm_e,Tp_e,Tm_e, Tp_v, Tm_v, Pp_v, Pm_v;

  double Pp1_e,Pm1_e,Tp1_e,Tm1_e, Tp1_v, Tm1_v, Pp1_v, Pm1_v;
  double Pp2_e,Pm2_e,Tp2_e,Tm2_e, Tp2_v, Tm2_v, Pp2_v, Pm2_v;

  double Pe,Te,Pv,Tv;
  double Pvv,Tvv;

  double Pe_o2,Te_o2,Pv_o2,Tv_o2;
  double Pe_o4,Te_o4,Pv_o4,Tv_o4;
  
  double Pvv_o2,Tvv_o2,Pvv_o4,Tvv_o4;

  double dE=1e-7*E;
  //double dV=1e-5*V;
  double dV=1e-10;

  double Vx, Ex, Vy, Ey;
  int phaseA, phaseB;
  double x=0;
  double fac;
  double Pa, Pb, c2, G, c2_VE, G_VE;
  double dPdV_S, d2PdV2_S;

  double V0x_scale, V0y_scale, E0x_scale, E0y_scale;

  // Calcul des dérivées partielles
  if(zone==1||zone==2||zone==3){
    P=fP(V,E,zone);  T=fT(V,E,zone);
    
    Pe=dPdE(zone);     Te=dTdE(zone);
    Pv=dPdV(V, zone);  Tv=dTdV(V, zone);
    
    Pvv=d2PdV2(V, zone);
    Tvv=d2TdV2(V, zone);

    //Pe=0; Te=0; Tv=0; Pv=0;
    fac=1.0;
  }
  else if(zone==12||zone==13||zone==23){
    phaseB=zone%10;  phaseA=(zone-phaseB)/10;
    
    /*
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    P=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    T=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E+2*dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pp2_e=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tp2_e=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;
    
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E+dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pp1_e=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tp1_e=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E-2*dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pm2_e=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tm2_e=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;
    
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E-dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pm1_e=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tm1_e=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+2*dV, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pp2_v=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tp2_v=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;
    
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pp1_v=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tp1_v=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;

    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-2*dV, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pm2_v=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tm2_v=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;
    
    err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
    if(err){printf("test_consis_thermo_eq \n"); return err;}
    Pm1_v=(fP(Vx,Ex,phaseA)+fP(Vy,Ey,phaseB))/2;
    Tm1_v=(fT(Vx,Ex,phaseA)+fT(Vy,Ey,phaseB))/2;

    Tv_o2=(Tp1_v - Tm1_v)/(2*dV);
    Pe_o2=(Pp1_e - Pm1_e)/(2*dE);
    Te_o2=(Tp1_e - Tm1_e)/(2*dE);
    Pv_o2=(Pp1_v - Pm1_v)/(2*dV);

    Tvv_o2=(Tp1_v + Tm1_v - 2*T)/(dV*dV);
    Pvv_o2=(Pp1_v + Pm1_v - 2*T)/(dV*dV);
    
    Tv_o4=(-Tp2_v + 8*Tp1_v - 8*Tm1_v + Tm2_v)/(12*dV);
    Pe_o4=(-Pp2_e + 8*Pp1_e - 8*Pm1_e + Pm2_e)/(12*dE);
    Te_o4=(-Tp2_e + 8*Tp1_e - 8*Tm1_e + Tm2_e)/(12*dE);
    Pv_o4=(-Pp2_v + 8*Pp1_v - 8*Pm1_v + Pm2_v)/(12*dV);

    Tvv_o4=(-Tp2_v + 16*Tp1_v -30*T + 16*Tm1_v - Tm2_v)/(12*dV*dV);
    Pvv_o4=(-Pp2_v + 16*Pp1_v -30*P + 16*Pm1_v - Pm2_v)/(12*dV*dV);
    */

    //ordre 2
    Tv=Tv_o2;  Pe=Pe_o2;  Te=Te_o2;  Pv=Pv_o2; Tvv=Tvv_o2; Pvv=Pvv_o2;

    //ordre 4
    //Tv=Tv_o4;  Pe=Pe_o4;  Te=Te_o4;   Pv=Pv_o4; Tvv=Tvv_o4; Pvv=Pvv_o4;

    fac=0.0;
  }
  else{  // zone triple
    // a completer et hors frontières
    //T=0; P=1; Te=0; Tv_o4=0; Pe=0;
    fac=0.0;
  }
  

  // Calcul classique (dans le code) de la celerite du son :
  if(1){
    dE=-P*dV;
    if (zone==1||zone==2||zone==3){
      Pb=fP(V+dV,E+dE,zone);
      Pa=fP(V-dV,E-dE,zone);
      
      dPdV_S=(Pb-Pa)/(2*dV);
      d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);
      
      c2=-V*V*dPdV_S;
      
      G=(-V/2)*d2PdV2_S/dPdV_S;
    }
    else if(zone==12||zone==13||zone==23){
      phaseB=zone%10;  phaseA=(zone-phaseB)/10;

      /*
      err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
      if(err){printf("test_consis_thermo_eq \n"); return err;}
      Pb=(fP(Vx,Ex,phaseA));//+fP(Vy,Ey,phaseB))/2;
      
      err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);  
      if(err){printf("test_consis_thermo_eq \n"); return err;}
      Pa=(fP(Vx,Ex,phaseA));//+fP(Vy,Ey,phaseB))/2;
      */
      
      dPdV_S=(Pb-Pa)/(2*dV);
      d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);

      c2=-V*V*dPdV_S;
      
      G=(-V/2)*d2PdV2_S/dPdV_S;
      G=0;
    }
    else{c2=0;}
  }

  // Calcul avec prise en compte des 3 phases comme dans le code actuel
  double* lambda=malloc(3*sizeof(double));
  int ind, phaseC;
  if(0){
    dV=1e-8;
    dE=-P*dV;
    if(zone==12||zone==13||zone==23){
      phaseB=zone%10;  phaseA=(zone-phaseB)/10; phaseC=-(phaseA+phaseB)+6;
      err=VE_mixte(Nmax, epsilon, phaseA, phaseB, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);
      if(err){printf("test_consis_thermo_eq \n"); return err;}

      //printf("zone=%d [%d,%d,%d]\n",zone,phaseA,phaseB,phaseC );

      lambda[phaseA-1]=1.-x; 
      lambda[phaseB-1]=x;
      lambda[phaseC-1]=0.0; 
    }
    else if(zone==1||zone==2||zone==3){      
      //printf("zone=%d zone/3=%d [%d,%d,%d]\n",zone,(zone)%3,zone,(zone)%3+1,(zone+1)%3+1 );

      lambda[zone-1]=1.;
      ind=(zone)%3;
      lambda[ind]=0.0;
      ind=(zone+1)%3;
      lambda[ind]=0.0;
    }

    double Vz,Ez;
    double Va0, Ea0, Vb0, Eb0, Vc0, Ec0;
    Va0=125e-6;  Ea0=2.5e4;
    Vb0=108e-6;  Eb0=3.4e5;
    Vc0=112e-6;  Ec0=6.8e5;

    err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V+dV, E+dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
    if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
    Pb=fP(Vx,Ex,1);

    err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V-dV, E-dE, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
    if(err){printf("erreur VE_mixte_cin_3_phase_all_dom\n"); return 1;}
    Pa=fP(Vx,Ex,1);

    dPdV_S=(Pb-Pa)/(2*dV);
    d2PdV2_S=(Pb+Pa-2*P)/(dV*dV);

    c2=-V*V*dPdV_S;  
    G=(-V/2)*d2PdV2_S/dPdV_S;

    if (zone==-1||zone==0){
      c2=0;
    }
  }

  // Autre test de consistance thermo :
  //*pres=fac*(T*Pe - P*Te + Tv);
  //*pres=fac*sqrt(-V*V*(Pv - P*Pe));  // test du calcul de c
  //pres[2]=sqrt(c2);
  //pres[3]=fac*sqrt(-V*V*(Pv - P*Pe));

  double dPdV_S_VE=Pv-P*Pe;
  double d2PdV2_S_VE=fac*(Pvv - Pv*Pe + P*Pe*Pe);

  c2_VE=-V*V*dPdV_S_VE;
  G_VE=(-V/2)*d2PdV2_S_VE/dPdV_S_VE;

  //pres[2]=G;
  //pres[3]=G_VE;

  pres[2]=fac*d2PdV2_S*1e-19;
  pres[3]=fac*d2PdV2_S_VE*1e-19;
  //pres[2]=c2;
  //pres[3]=c2_VE;

  //printf("  P=%g, T=%g, zone=%d\n", P,T,zone);
  //printf("  T.Pe=%g, P.Te=%g, Tv=%g test=%g\n", T*Pe,P*Te,Tv,*pres);

  free(lambda);
  return 0;
}


int fPTC_cin_phase(int affichage, int EOS, int CIN_PHASE, int sch_cin_phase, double epsilon, double V, double E, double* pP, double *pc, double* pT, double* pVAR, double* lambda_in, double* lambda_out,
                   int Ntau, double dx, double dt, double Ttau, int* pzone, int* n_iter, double* pcritere, double* p_test_consis_thermo,
                   double dUdx, double Pn, double Vn, double Cn, // cinétique
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC,
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){

  // Ttau est calculé avec c(Vn, En)
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);  }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);

  double old_critere=*pcritere;
  double old_nbiter=*n_iter;

  int Nmax=50, err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y, V0z, E0z;
  double Vx, Ex, Vy, Ey, Vz, Ez;
  double dV=1e-8;

  int zone, phaseA, phaseB;
  double G,S;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple double Va_temp,Vb_temp,Vc_temp;
  double PT,TT;
  double VaT,VbT,VcT;
  double EaT,EbT,EcT;

  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);
  /*  ancien coeff
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  */

  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  double Peq, NUM, DEN;


  // Tabulation pure (ne produit que l'effet mais n'est pas optimisé)
  // on se contente de modifire le couple (V,E) d'entreé pour qu'il coresponde au centre de l'élement de tabulation 
  // puis calcul les granduers thermo a precision machine
  double hv=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he=(Emax_dis-Emin_dis)/(Ne_dis-1);
  int i=(V-Vmin_dis)/hv;
  int j=(E-Emin_dis)/he;
  if(EOS==2){ 
    V = Vmin_dis + (i+0.5)*hv;
    E = Emin_dis + (j+0.5)*he;
  }
  else{
    // rien on garde le couple V,E d'entreé
  }


  // Fonction renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  err= lec_benchmark_zone(affichage, EOS, epsilon, Nmax, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y, pP, pT, n_iter, pcritere);
  if(err){printf("lec_benchmark_zone depuis fPTC_cin_phase\n");  return err;}
  
  *pzone=zone;
  *pcritere=fmax(*pcritere, old_critere);
  *n_iter=fmax(*n_iter, old_nbiter);

  // Calcul de l'entropie à l'équilibre
  if(zone==1||zone==2||zone==3){
    S=fS(V,E,zone);
  }
  else if(zone==12||zone==13||zone==23){
    phaseB=zone%10;  phaseA=(zone-phaseB)/10;
    S=Xfm[phaseA-1]*fS(V0x,E0x,phaseA)+(Xfm[phaseB-1])*fS(V0y,E0y,phaseB);
  }
  else if(zone==0){
    S=Xfm[0]*fS(Va0,Ea0,1);
    S+=Xfm[1]*fS(Vb0,Eb0,2);
    S+=Xfm[2]*fS(Vc0,Ec0,3);
  }

  //printf(" ZONE= %d\n",zone );
  for(int i=0; i<3; ++i){  if( isnan(Xfm[i]) ){ printf("Erreur lambda=[%g, %g, %g]\n",Xfm[0],Xfm[1],Xfm[2] ); return 1;}  }
  
  if(affichage){}
  else{
    for(int i=0; i<3; i++){
      if(Xfm[i]<0 ){
        printf("lambda=[%g, %g, %g] V=%g E=%g (fPTC_cin_phase)\n",Xfm[0],Xfm[1],Xfm[2],V,E );
        Xfm[i]=0;
      }
      if(Xfm[i]>1 ){
        printf("lambda=[%g, %g, %g] V=%g E=%g (fPTC_cin_phase)\n",Xfm[0],Xfm[1],Xfm[2],V,E );
        Xfm[i]=1;
       }
     }
  }
  
  //printf(" V= %g , E= %g , zone= %d, cin_phase= %d \n",V,E,zone,CIN_PHASE );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  
  
  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.; *pc=0.; pVAR[0]=0.; pVAR[1]=S; *p_test_consis_thermo=0; lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;   return 0; } //printf("     hors domaine zone=%d\n",zone );
    else if(affichage==0){printf("Error out of domain : V= %g, E= %g (fPTC_cin_phase)\n",V,E); return 1;}
    else{printf("Erreur valeur affichage= %d (fPTC_cin_phase)",affichage); return 1;} 
  }

  //printf(" ZONE= %d | xeqA= %g , xeqB= %g, xeqC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,Xfm[0],Xfm[1],Xfm[2],lambda_out[0],lambda_out[1],lambda_out[2] );  

  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);
  if(CIN_PHASE==0){  // equilibre thermo
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; }
    if(zone==0){ *pc=0; pVAR[0]=0.; pVAR[1]=S; *p_test_consis_thermo=0; return 0; } 

    //err=c_G_eq(Nmax, epsilon, V, E, zone, *pP, *pT, V0x, E0x, V0y, E0y, pc, &G);
    Vx=125e-6;  Ex=2.5e4;
    Vy=108e-6;  Ey=3.4e5;
    Vz=112e-6;  Ez=6.8e5;
    err=c_G_3_phase_all_dom(Nmax, epsilon, Vx, Ex, Vy, Ey, Vz, Ez, *pP, V, E, dV, lambda_out, pc, &G);  
    if(err){printf("Erreur c_G_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}
    
    pVAR[0]=G;
    pVAR[1]=S;
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    err= lambda_cin_phase(sch_cin_phase, lambda_out, lambda_in, Xfm, dt, Ttau );
    if(err){ printf("Erreur lambda_cin_phase (fPTC_cin_phase)\n"); return 1; }

    err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda_out, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V, E, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
    if(err){printf("erreur VE_mixte_cin_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}
    *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);
    S=lambda_out[0]*fP(Vx,Ex,1) + lambda_out[1]*fP(Vx,Ex,2) + lambda_out[2]*fP(Vx,Ex,3);  
    pVAR[1]=S;

    err=c_G_3_phase_all_dom(Nmax, epsilon, Vx, Ex, Vy, Ey, Vz, Ez, *pP, V, E, dV, lambda_out, pc, &G);  
    if(err){printf("Erreur c_G_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}

    pVAR[0]=G;
  }
  else if(CIN_PHASE==3){ // dtP + rho.c2.div(u) = (-1/tau)(P-Peq)
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; }

    Peq=*pP;
    NUM= (-Cn*Cn/Vn)*dUdx + Peq/Ttau;
    DEN= 1+dt/Ttau;
    *pP= (Pn + dt*NUM)/DEN;

    err=c_celerite_son_cin_3_phase_all_dom(Nmax, epsilon, Peq, V, E, dV, lambda_out, pc);  if(err){printf("Erreur c_celerite_son_cin_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}
  }
  else{printf("Erreur CIN_PHASE= %d (fPTC_cin_phase)\n",CIN_PHASE ); return 1;}


  if(0){
    err= test_consis_thermo_eq(Nmax, epsilon, V, E, zone, V0x, E0x, V0y, E0y, p_test_consis_thermo); if(err){printf("Erreur test_consis_thermo_eq (fPTC_cin_phase)\n"); return err;}
  }

  return 0;
}


int fPTC_cin_phase_test_Newton(int affichage, int EOS, int CIN_PHASE, int sch_cin_phase, double epsilon, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out,
                   int Ntau, double dx, double dt, double Ttau, int* pzone, int* n_iter, double* pcritere, double* p_test_consis_thermo,
                   double dUdx, double Pn, double Vn, double Cn, // cinétique
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC,
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){

  double old_critere=*pcritere;
  double old_nbiter=*n_iter;

  int Nmax=50;
  //double epsilon=1e-13;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y, V0z, E0z;
  double Vx, Ex, Vy, Ey, Vz, Ez;
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double Va_temp,Vb_temp,Vc_temp;
  double Ea_temp,Eb_temp,Ec_temp;
  double PT,TT;
  point_triple(&PT, &TT, &Va_temp, &Ea_temp, &Vb_temp, &Eb_temp, &Vc_temp, &Ec_temp);


  double Va0, Ea0, Vb0, Eb0, Vc0, Ec0; 
  Va0=125e-6;  Ea0=2.5e4;
  Vb0=108e-6;  Eb0=3.4e5;
  Vc0=112e-6;  Ec0=6.8e5;

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  err= lec_benchmark_zone(affichage, EOS, epsilon, Nmax, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                          Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y, pP, pT, n_iter, pcritere);
  if(err){printf("lec_benchmark_zone depuis fPTC_cin_phase\n");  return err;}
  
  p_test_consis_thermo[0]=*pP;
  p_test_consis_thermo[1]=*pT;
  // p_test_consis_thermo[2]=C

  *pzone=zone;

  *pcritere=fmax(*pcritere, old_critere);
  *n_iter=fmax(*n_iter, old_nbiter);

  //printf(" ZONE= %d\n",zone );
  for(int i=0; i<3; ++i){  if( isnan(Xfm[i]) ){ printf("Erreur lambda=[%g, %g, %g]\n",Xfm[0],Xfm[1],Xfm[2] ); return 1;}  }
  
  if(affichage){}
  else{
    for(int i=0; i<3; i++){
     if(Xfm[i]<0 ){
      printf("lambda=[%g, %g, %g] V=%g E=%g (fPTC_cin_phase)\n",Xfm[0],Xfm[1],Xfm[2],V,E );
      Xfm[i]=0;
     }
     if(Xfm[i]>1 ){
      printf("lambda=[%g, %g, %g] V=%g E=%g (fPTC_cin_phase)\n",Xfm[0],Xfm[1],Xfm[2],V,E );
      Xfm[i]=1;
     }
    }
  }
  
  //printf(" V= %g , E= %g , zone= %d, cin_phase= %d \n",V,E,zone,CIN_PHASE );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  
  
  if(zone==-1){
    if(affichage==1){*pP=-0; *pT=0.; *pc=0.; p_test_consis_thermo[0]=0; p_test_consis_thermo[1]=0; p_test_consis_thermo[2]=0; p_test_consis_thermo[3]=0; lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;   return 0; } //printf("     hors domaine zone=%d\n",zone );
    else if(affichage==0){printf("Error out of domain : V= %g, E= %g (fPTC_cin_phase)\n",V,E); return 1;}
    else{printf("Erreur valeur affichage= %d (fPTC_cin_phase)",affichage); return 1;} 
  }

  for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; }

  //printf(" ZONE= %d | xeqA= %g , xeqB= %g, xeqC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,Xfm[0],Xfm[1],Xfm[2],lambda_out[0],lambda_out[1],lambda_out[2] );  

  if(zone==0){ *pc=0; p_test_consis_thermo[2]=0; p_test_consis_thermo[3]=0; return 0; } 

  err= VE_mixte_cin_3_phase_all_dom(Nmax, epsilon, lambda_out, Va0, Ea0, Vb0, Eb0, Vc0, Ec0, V, E, &Vx, &Ex, &Vy, &Ey, &Vz, &Ez);
  if(err){printf("erreur VE_mixte_cin_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}
  *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);
  if( (fP(Vx,Ex,2)-fP(Vx,Ex,3)) >=1e3){
    printf(" Erreur Pbeta=%g, Pgamma=%g, Pliq=%g\n",fP(Vx,Ex,1),fP(Vx,Ex,2),fP(Vx,Ex,3) );
    printf("        Tbeta=%g, Tgamma=%g, Tliq=%g\n",fT(Vx,Ex,1),fT(Vx,Ex,2),fT(Vx,Ex,3) );
    printf("        V=%g, E=%g | (V,E)beta=%g; %g; (V,E)gamma=%g; %g; (V,E)liq=%g; %g;\n",V,E,Vx,Ex,Vy,Ey,Vz,Ez );
    return 1;
  }

  // cérité du son 
  err=c_celerite_son_cin_3_phase_all_dom(Nmax, epsilon, *pP, V, E, dV, lambda_out, pc);  if(err){printf("Erreur c_celerite_son_cin_3_phase_all_dom (fPTC_cin_phase)\n"); return 1;}
  

  if (affichage==1){
    err= test_consis_thermo_eq(Nmax, epsilon, V, E, zone, V0x, E0x, V0y, E0y, p_test_consis_thermo); if(err){printf("Erreur test_consis_thermo_eq (fPTC_cin_phase)\n"); return err;}
  }

  return 0;
}

// fonction ENERGIE INTERNE SPECIFIQUE epsilon(tau, p)
// Algorithme de Newton pour calculer EPSn+1 sans conaitre Pn+1
// Necessaire pour le schema vNR
int fPn1_etain_cin_phase(int tst, int EOS, int cin_phase, int sch_cin_phase, double epsilon, double TAUn, double TAUn1, double Pn, double Qn, double EPSn, double* pEPSn1, double* pPn1, double* lambda_in, int Ntau, double dX, double dt, double Ttau,
                         int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                         double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                         double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){

  int Nmax=1e3;
  double eps=1e-9;
  if(tst/100==1){  eps=1e-12;  }
  else if(tst/100==0){  eps=1e-10;  }

  //double epsilon=1e-7;

  int err, n_iter;
  double res;
  double dTAU=TAUn1-TAUn; 
  double dE;
  double E=EPSn;
  double C;
  double Cunuse, Vunuse, Punuse, dUdx_unuse;
  double test_consis_thermo;

  //printf(" V= %g , E= %g , TAUn1= %g , Pn= %g , xA= %g , xB= %g , xC= %g\n",TAUn,EPSn,TAUn1,Pn,lambda[0],lambda[1],lambda[2] );
  
  // Calcul de la derivée
  double dE_inf=EPSn*1e-5;
  double P,T,x;
  double* pVAR=malloc(2*sizeof(double));
  int zone;
  double Pp,Pm;
  int phaseA,phaseB;
  double dPpdE, dPmdE, dPdE;
  double V0x, E0x, V0y, E0y, Vx, Ex, Vy, Ey;

  double hv=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he=(Emax_dis-Emin_dis)/(Ne_dis-1);
  int i=(TAUn1-Vmin_dis)/hv;
  int j=(EPSn-Emin_dis)/he;

  int pos_fic=i*(Ne_dis-1) + j;
  //printf("  --i_dis= %d, j_dis= %d \n",i,j );

  V0x=tabVx[pos_fic];  E0x=tabEx[pos_fic];
  V0y=tabVy[pos_fic];  E0y=tabEy[pos_fic];

  // Calcul de la derivee
  double lambdaB;
  double* lambda_test=malloc(3*sizeof(double)); // récupére la modif de lambda dans fPTC_cin_phase (inutile)
  lambda_test[0]=lambda_in[0];
  lambda_test[1]=lambda_in[1];
  lambda_test[2]=lambda_in[2];


  double f,f_prime;
  double critere=1.;
  int n=0;
  while(critere>eps && n<Nmax){
    err=fPTC_cin_phase(0, EOS, 0, sch_cin_phase, epsilon, TAUn1, E, &P, &C, &T, pVAR, lambda_test, lambda_test, 
                       Ntau, dX, dt, Ttau, &zone, &n_iter, &critere, &test_consis_thermo, dUdx_unuse, Punuse, Vunuse, Cunuse,
                       nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                       SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
    if(err){return err;}

    if(zone==-1){printf("Erreur out of domain : zone= %d (fPn1_etain_cin_phase)\n",zone); return 1;  }
    else if(zone==1 || zone==2 || zone==3){
      Pp=fP(TAUn1,E+dE_inf,zone);
      Pm=fP(TAUn1,E-dE_inf,zone);
    }
    else if(zone==12 || zone==13 || zone==23){
      j=(E-Emin_dis)/he;
      pos_fic=i*(Ne_dis-1) + j;
      V0x=tabVx[pos_fic];  E0x=tabEx[pos_fic];
      V0y=tabVy[pos_fic];  E0y=tabEy[pos_fic];

      phaseB=zone%10;
      phaseA=(zone-phaseB)/10;
      lambdaB=lambda_test[phaseB-1];

      err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, TAUn1, E+dE_inf, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
      Pp=fP(Vx,Ex,phaseA);

      err=VE_mixte_cin_phase(Nmax, epsilon, phaseA, phaseB, lambdaB, V0x, E0x, V0y, E0y, TAUn1, E-dE_inf, &Vx, &Ex, &Vy, &Ey);  if(err){return 1;}
      Pm=fP(Vx,Ex,phaseA);
    }
    else if(zone==0){
      Pp=P;
      Pm=P;
    }

    dPpdE=(Pp-P)/dE_inf;
    dPmdE=(Pm-P)/(-dE_inf);
    dPdE=(dPpdE + dPmdE)/2;

    f= E-EPSn + ((P+Pn)/2 + Qn)*dTAU;
    f_prime= 1 + dPdE*dTAU/2;
    dE=(1.)*(-f/f_prime);
    if(isnan(dE)){printf("Erreur dE=%g (fPn1_etain_cin_phase)\n",dE ); return 1;}

    E+=dE;

    //critere=fabs(E-EPSn + ((P-Pn)/2 + Qn)*dTAU);
    critere=fabs(dE/E);
    
    n++;
  }
  //printf("  n= %d , dE= %g, critere= %g , |f|= %g , zone= %d \n",n,dE,critere,fabs(f),zone );
  if(n>=Nmax){printf("  n= %d , critere= %g , |f|= %g , zone= %d \n",n,critere,fabs(f),zone ); return 1;}

  *pEPSn1=E;
  *pPn1=P;

  free(lambda_test);

  return 0;
}


// Script pour tester le Benchmark
int use_benchmark_PT_cin_phase(int EOS, double epsilon, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err, n_iter;
  double E,V;
  double P,T;
  double x;
  double dt;
  //double test_consis_thermo;
  double* test_consis_thermo=malloc(4*sizeof(double));
  double* pVAR=malloc(2*sizeof(double));
  double critere;

  int sch_cin_phase;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  double* lambda=malloc(3*sizeof(double));
  double Ttau;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  /* pour pouvoir utiliser la discretisation juste en lisant les fichiers
  Organisation du fichier :
    nb_points
    Nv_dis , Ne_dis
    Vmin_dis , Vmax_dis
    Emin_dis , Emax_dis
  */
  
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);



 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n"); return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n"); return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n"); return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n"); return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n"); return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n"); return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.14g %.14g\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.14g %.14g\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.14g %.14g\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.14g %.14g\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.14g %.14g\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.14g %.14g\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////


  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXAtest, *fXBtest, *fXCtest, *fCtest, *fGtest, *fStest, *ftest_cons_thermo;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXAtest = fopen("fichiers/XAtest.txt", "w+")) == NULL){printf("erreur ouverture fichier XA test\n");return 1;}
  if((fXBtest = fopen("fichiers/XBtest.txt", "w+")) == NULL){printf("erreur ouverture fichier XB test\n");return 1;}
  if((fXCtest = fopen("fichiers/XCtest.txt", "w+")) == NULL){printf("erreur ouverture fichier XC test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier c test\n");return 1;}
  if((fGtest = fopen("fichiers/Gtest.txt", "w+")) == NULL){printf("erreur ouverture fichier G test\n");return 1;}
  if((fStest = fopen("fichiers/Stest.txt", "w+")) == NULL){printf("erreur ouverture fichier S test\n");return 1;}
  if((ftest_cons_thermo = fopen("fichiers/test_cons_thermo.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      fprintf(fEtest, "%.14g ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin_test+j*hv_test;
      fprintf(fVtest, "%.14g ",V);
  }
  fclose(fVtest);  fclose(fEtest);


  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.14g \n",Vmin_dis+i*hv_dis+hv_dis/2 );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.14g \n",Emin_dis+i*he_dis+he_dis/2 );
  }
  fclose(file_E_dis);
  
  // Arguement pour cinétique non utilisé 
  double Cunuse, Vunuse, Punuse, dUdx_unuse;
  double Ntau=1, dx=1;

  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      //printf("      - i= %d, j= %d ",i,j );
      //printf("   V= %g  E= %g\n",V,E );
      //err=lec_benchmark_PT(1,V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
      //                     Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
      
      //fPTC_cin_phase_test_Newton
      err=fPTC_cin_phase(1, EOS, 0, sch_cin_phase, epsilon, V, E, &P, &c, &T, pVAR, lambda, lambda, 
                         Ntau, dx, dt, Ttau, &zone, &n_iter, &critere, test_consis_thermo, dUdx_unuse, Punuse, Vunuse, Cunuse,
                         nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                         SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){printf(" | fPTC_cin_phase (use_benchmark_PT_cin_phase)\n"); return err;}
      
      fprintf(fPtest,"%.14g\n",P);
      fprintf(fTtest,"%.14g\n",T);
      fprintf(fXAtest,"%.14g\n",lambda[0]);
      fprintf(fXBtest,"%.14g\n",lambda[1]);
      fprintf(fXCtest,"%.14g\n",lambda[2]);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.14g\n",c);      
      fprintf(fGtest,"%.14g\n",pVAR[0]);      
      fprintf(fStest,"%.14g\n",pVAR[1]);      
      fprintf(ftest_cons_thermo,"%.14g %.14g %.14g %.14g \n",test_consis_thermo[0],test_consis_thermo[1],test_consis_thermo[2],test_consis_thermo[3]);   
    } 
  }

  fclose(fPtest);  fclose(fTtest);  fclose(fXAtest); fclose(fXBtest); fclose(fXCtest); fclose(fGtest); fclose(fStest); fclose(fCtest); fclose(fZONEtest); fclose(ftest_cons_thermo);
  

  // Liberation de la mémoire
  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}



int debug_Newton(){
  printf("\n\n  DEBUG NEWTON  \n");
  int Nmax=1e2;
  double epsilon=1e-5;
  int err, n_iter;
  double critere;
   
  double x;
  double V,E;
  V=0.0001; E=831111;

  double Vx, Ex, Vy, Ey;
  double V0x, E0x, V0y, E0y;
  V0x= 9.92845e-05; E0x= 806926; V0y= 0.000101004; E0y= 866875;
  //V0y*=(1.0+1e-3);
  //E0y*=(1.0+1e-3);
  printf(" V0x= %g , E0x= %g , V0y= %g , E0y= %g\n",V0x,E0x,V0y,E0y );

  err=VE_mixte(Nmax, epsilon, 2, 3, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey, &x, &n_iter, &critere);
  if(err){return 1;}
  printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey );
  printf(" x= %g ,\n",x );
  
  double lambda3=x;


  //V0x= Vx; E0x= Ex; V0y= Vy; E0y= Ey;
  lambda3*=1.3;
  printf("lambda= %g\n",lambda3 );
  err= VE_mixte_cin_phase(Nmax, epsilon, 2, 3, lambda3, V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
  if(err){return 1;}

  printf("Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey );
  printf("P= %g  T= %g \n",fP(Vx,Ex,2), fT(Vx,Ex,3));
}



// IMPLEMENTATION DES METHODES TABULATIONS OU CALCUL PUR


// Fonction pour avoir une benchmark de (P,T) sur (V,E)
/*

*/
/*
int ecr_benchmark_pur_METH2(int nb_points, int Nv ,int Ne, double Vmin, double Vmax, double Emin, double Emax){

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisations juste en lisant les fichiers
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fprintf(fparams, "%d \n",nb_points );
  fprintf(fparams, "%d %d\n",Nv,Ne );
  fprintf(fparams, "%.15lf %.15lf\n",Vmin,Vmax );
  fprintf(fparams, "%.15lf %.15lf\n",Emin,Emax );
  fprintf(fparams, "\n\n Organisation du fichier\n");
  fprintf(fparams, "nb_points\n");
  fprintf(fparams, "Nv_dis  Ne_dis\n");
  fprintf(fparams, "Vmin_dis Vmax_dis\n");
  fprintf(fparams, "Emin_dis Emax_dis\n");


  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  

  int err,res;
  

  double E,V;

  double Vmx, Emx;
  double Vmy, Emy;


  // Varaibles pour discretisation de la maille
  double V_sub, E_sub;
  int N_rec=10;
  double he_sub=he/(N_rec-1);
  double hv_sub=hv/(N_rec-1);
  double** VE_front_rec=alloctabd(4*(N_rec-1),2);
  int count;
  

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );
  

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);
  
   
  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));
  double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));
  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));
  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));
  double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));
  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));
  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));
  double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));
  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));
  double* segEc23=malloc(N23*sizeof(double));
  


  // CALCUL DES SOMMETS  
  printf(" Sommets_polygone \n");
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}


  
  // Pour tester les valeurs (coord des points du maillage)
  FILE *fileV, *fileE;
  if((fileV = fopen("fichiers/Vbench_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier V bench \n"); return 1;}
  if((fileE = fopen("fichiers/Ebench_METH2.txt", "w+")) == NULL){printf("erreur ouverture fichier E bench \n"); return 1;}
  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fileE, "%.15lf ",E);
  }
  for(int j=0; j<Nv; j++){
      V=Vmin+j*hv;
      fprintf(fileV, "%.15lf ",V);
  }
  fclose(fileV);  fclose(fileE);
  /////////////
  

  double* tabVx=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEx=malloc((Nv-1)*(Ne-1)*sizeof(double));
  double* tabVy=malloc((Nv-1)*(Ne-1)*sizeof(double)); double* tabEy=malloc((Nv-1)*(Ne-1)*sizeof(double));
  
  int* tabZONE=malloc((Nv-1)*(Ne-1)*sizeof(double));  // mémorise la zone dans le rectangle


  int test;

  int k=0;
  for(int i=0; i<Nv-1; i++){
    printf("  --i= %d\n",i );
    V=Vmin+i*hv + hv/2; // pour prendre la valeur au centre du rectangle
    for(int j=0; j<Ne-1; j++){
      //printf("  -- i= %d, j= %d\n",i,j );

      E=Emin+j*he + he/2;  
      

      // Attribution d'une zone aux points milieux des rectangles
      if(intPOLY(V,E,NC,SC)){ 
        tabZONE[k]=3;
      }
      else if(intPOLY(V,E,NB,SB)){
        tabZONE[k]=2;
      }      
      else if(intPOLY(V,E,NA,SA)){
        tabZONE[k]=1;
      }
      else if(ftri(V,E)){
        tabZONE[k]=0;
      }
      else if(intPOLY(V,E,NAB,SAB)){
        err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=12;
      }
      else if(intPOLY(V,E,NAC,SAC)){
        err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=13;
      }
      else if(intPOLY(V,E,NBC,SBC)){
        err=graine_zone_mixte(V, E, N23, segVb23, segVc23, segEb23, segEc23, &Vmx, &Emx, &Vmy, &Emy);
        if(err){printf("err= %d (benchmark_PT)\n",err ); return err;}
        tabVx[k]=Vmx; tabEx[k]=Emx;
        tabVy[k]=Vmy; tabEy[k]=Emy;
        tabZONE[k]=23;
      }
      else{  // hors domaine --> on ne se soucis des pb de chevauchement des rectangles sur les zones
        tabZONE[k]=-1;
      }
    }
  }


  // fichier binaire   couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  FILE *fileVmy, *fileEmy;
  if((fileVmy = fopen("fichiers/Vbench_mixte_y_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabVx,sizeof(double),(Nv-1)*(Ne-1),fileVmx);
  fwrite(tabEx,sizeof(double),(Nv-1)*(Ne-1),fileEmx);
  fwrite(tabVy,sizeof(double),(Nv-1)*(Ne-1),fileVmy);
  fwrite(tabEy,sizeof(double),(Nv-1)*(Ne-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // fichier binaire ZONE
  FILE *fileZONE;
  if((fileZONE = fopen("fichiers/ZONEbench_METH2.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE bench \n"); return 1;}
  
  // Ecriture des tableaux
  fwrite(tabZONE,sizeof(int),(Nv-1)*(Ne-1),fileZONE);

  fclose(fileZONE_0);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE A (ecr_benchmark_PT)\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE B (ecr_benchmark_PT)\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE C (ecr_benchmark_PT)\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AB (ecr_benchmark_PT)\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE AC (ecr_benchmark_PT)\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "wb")) == NULL){printf("erreur ouverture fichier ZONE BC (ecr_benchmark_PT)\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  for(int i=0; i<NA+1; i++){  SA_V[i]=SA[i][0];  SA_E[i]=SA[i][1];  }
  for(int i=0; i<NB+1; i++){  SB_V[i]=SB[i][0];  SB_E[i]=SB[i][1];  }
  for(int i=0; i<NC+1; i++){  SC_V[i]=SC[i][0];  SC_E[i]=SC[i][1];  }

  for(int i=0; i<NAB+1; i++){  SAB_V[i]=SAB[i][0];  SAB_E[i]=SAB[i][1];  }
  for(int i=0; i<NAC+1; i++){  SAC_V[i]=SAC[i][0];  SAC_E[i]=SAC[i][1];  }
  for(int i=0; i<NBC+1; i++){  SBC_V[i]=SBC[i][0];  SBC_E[i]=SBC[i][1];  }

  // Ecriture des tableaux
  fwrite(SA_V,sizeof(double),NA+1,fileA_V);     fwrite(SA_E,sizeof(double),NA+1,fileA_E);
  fwrite(SB_V,sizeof(double),NB+1,fileB_V);     fwrite(SB_E,sizeof(double),NB+1,fileB_E);
  fwrite(SC_V,sizeof(double),NC+1,fileC_V);     fwrite(SC_E,sizeof(double),NC+1,fileC_E);
  fwrite(SAB_V,sizeof(double),NAB+1,fileAB_V);  fwrite(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fwrite(SAC_V,sizeof(double),NAC+1,fileAC_V);  fwrite(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fwrite(SBC_V,sizeof(double),NBC+1,fileBC_V);  fwrite(SBC_E,sizeof(double),NBC+1,fileBC_E);


  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);



  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);
  freetab(REC);

  free(tabNbZONE);
  freetab(tabINTER);
  freetab(tabInd_Mixte);
  freetab(numZONE);
  freetab(VE_front_rec);


  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  freetab(SA);   freetab(SB);   freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  freetab(Sa12); freetab(Sb12);
  freetab(Sa13); freetab(Sc13);
  freetab(Sb23); freetab(Sc23);

  freetab(STRI_12);  freetab(STRI_13);  freetab(STRI_23);

  fprintf(fparams, "\nEcriture bien terminé\n");
  fclose(fparams);


  return 0;

}
*/

/*  Arguments d'entrée :
     - V,E
     - nb_points : nombre de points de discretisation des frontières
     - SA,SB,SC etc.. contiennent les discretisations des frontières
    Arguments de sorties :
     - Xfm,pzone : pointeurs sur les fraction massique et la zone de (V,E)
     - pV0x, pE0x, etc...  pointeur sur la valeur initiale choisi pour lancer une zone mixte (On peut y avoir acces à la fonction avant grace à tabVx)
*/
int cal_zone(int affichage, double V, double E, int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                double* Xfm, int* pzone, double* pV0x, double* pE0x, double* pV0y, double* pE0y){
 
  int Nmax=20;
  double epsilon=1e-10;

  int err, n_iter; 
  double critere;
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;


  // pour le calcul des fractions massiques dans la zone trinagle --> calcul des coord barycentriques dans le triangle
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  double* VE=malloc(2*sizeof(double));
  VE[0]=V; VE[1]=E;
  double** TRI=alloctabd(3,2);
  TRI[0][0]=VaT; TRI[0][1]=EaT;
  TRI[1][0]=VbT; TRI[1][1]=EbT;
  TRI[2][0]=VcT; TRI[2][1]=EcT;


  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );

  
  for(int i=0; i<3; i++){ Xfm[i]=0; }

  //printf(" nb_zone= %d\n",nb_zone );
  if(intPOLY(V,E,NC,SC)){   Xfm[2]=1.0;  *pzone=3;  return 0;  }
  else if(intPOLY(V,E,NB,SB)){ Xfm[1]=1.0;  *pzone=2;  return 0;  }
  else if(intPOLY(V,E,NA,SA)){ Xfm[0]=1.0;  *pzone=1;  return 0;  }
  else if(intPOLY(V,E,NBC,SBC)){
    err= graine_zone_mixte_modif(V, E, N23, SBC, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d BC (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, &n_iter, &critere);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[2]=x;  Xfm[1]=1.0-x;  *pzone=23;

    return 0;
  }
  else if(intPOLY(V,E,NAC,SAC)){
    //err=graine_zone_mixte(V, E, N13, segVa13, segVc13, segEa13, segEc13, &Vmx, &Emx, &Vmy, &Emy);
    err= graine_zone_mixte_modif(V, E, N13, SAC, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d AC (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, &n_iter, &critere);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[2]=x;  Xfm[0]=1.0-x;  *pzone=13;

    return 0;
  }
  else if(intPOLY(V,E,NAB,SAB)){
    //err=graine_zone_mixte(V, E, N12, segVa12, segVb12, segEa12, segEb12, &Vmx, &Emx, &Vmy, &Emy);
    err= graine_zone_mixte_modif(V, E, N12, SAB, &Vmx, &Emx, &Vmy, &Emy);
    if(err){printf("err= %d AB (graine_zone_mixte_modif)\n",err ); return err;}
    err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, pV0x, pE0x, pV0y, pE0y, &x, &n_iter, &critere);
    if(err){printf("err= %d (VE_mixte)\n",err ); return 1;}
    Xfm[1]=x;  Xfm[0]=1.0-x;  *pzone=12;

    return 0;
  }
  else if(ftri(V,E)){
    err=coorbary(VE, TRI, &Xfm); if(err){return 1;} 
    *pzone=-1;  return 0;
  }
  else{
   if(affichage==1){Xfm[0]=-1.0;  Xfm[1]=-1.0;  Xfm[2]=-1.0;  *pzone=-1;  return 0;}
   else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
   else{printf("Erreur dans la valeur de affichage= %d",affichage); return 1;} 
  }

  return 1; 
}



// Renvoie la valeur de P,T et C pour un couple V,E en prenant en compte la cinétique de chgmt de phase
/*  Argments d'entrees :
      - CIN_PAHSE :  if CIN_PHASE==0 : pas de cinétique de phase
                     if CIN_PHASE==1 : prise en compte de la cinétique de phase

   Metre un param entier tels que 
    meth==0 : calcul + tabulation
    meth==1 : calcul pure
    meth==2 : tabulation pure
    meth==3 : maillage ??

    Pour le moment uniquement avec cal_tab
*/
/*
int fPTC_METH_cin_phase(int methode, int affichage, int CIN_PHASE, double V, double E, double* pP, double *pc, double* pT, double* lambda_in, double* lambda_out, int Ntau, double dx, double dt, double Ttau, int* pzone,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE ){
  
  //      Ttau est calculé avec c(Vn, En)   
  
  if(CIN_PHASE==2){  Ttau=Ntau*dx/(*pc);   }
  //double Ttau=1e-9;   //Ntau*dx/(*pc);

  int Nmax=3e2;
  double epsilon=2e-9;
  int err;
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;  
  double Vx, Ex, Vy, Ey;  
  double dV=1e-8;

  int zone;

  double* Xfm=malloc(3*sizeof(double)); // tableau des fractions massiques
  int phaseX,phaseY;

  // Point triple
  double PT= 4510421781.609182357788085;
  double TT= 624.63982538492518;

  //err= PT_tab_cal(V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
  //                Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, pP, pT, px, &zone, &V0x, &E0x, &V0y, &E0y);

  // Focntion renvoie la zone, la fraction massique et les graines pour les Newton_mixte
  
  if(methode==0){ // methode mixte
    err= lec_benchmark_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
                            Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  }
  else if(methode==1){ // calcul pur
    err= cal_zone(affichage, V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, Xfm, &zone, &V0x, &E0x, &V0y, &E0y);
  }
  else if(methode==1){ // tabulation pure
  if(err){return err;}
  
  *pzone=zone;


  //printf(" ZONE= %d\n",zone );
  
  //for(int i=0; i<3; i++){ if(Xfm[i]<0 ){Xfm[i]=0;}  if(Xfm[i]>1 ){Xfm[i]=1;} }
  
  //printf("  V= %g , E= %g , zone= %d\n",V,E,zone );
  //printf(" V= %g, E= %g, Xfm[0]= %g, Xfm[2]= %g, Xfm[2]= %g \n",V,E,Xfm[0],Xfm[1],Xfm[2] );  

  if(zone==-1){
    if(affichage==1){*pP=-10e9; *pT=0.;  lambda_out[0]=-1; lambda_out[1]=-1; lambda_out[2]=-1;  return 0;}
    else if(affichage==0){printf("Error out of domain : V= %g , E= %g\n",V,E); return 1;}
    else{printf("Erreur affichage= %d",affichage); return 1;} 
  }
  if(CIN_PHASE==0){
    for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i];  }
    // if mehode==2 (tabulation pure) return 0   car on a deja toutes les valeurs à l'équilibre 
    // en revanche si on active la cinétique de chgmt de phase, on doit recalculer P,T et C si on est en zone de cgmt de phase 
  }
  else if(CIN_PHASE==1 || CIN_PHASE==2){
    // Calcul de x^{n+1} 
    if(zone==1 || zone==2 || zone==3){ for(int i=0; i<3; i++){  lambda_out[i]=Xfm[i]; } }
    else if(zone==0){
      //printf("DlambdaA= %g , DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1),(Ttau*(lambda[1])+Xfm[1])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]-lambda_out[0];

    }
    else if(zone==12 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=1.-lambda_out[0]; // (autre phase)
      lambda_out[2]=0.;
      //printf(" ok 12  lambda[0]= %g , Xfm[0]= %g\n",lambda[0],Xfm[0]);
    }
    else if( zone==13 ){
      //printf("DlambdaA= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=((dt/Ttau)*Xfm[0]+lambda_in[0])/((dt/Ttau)+1.);
      lambda_out[1]=0;
      lambda_out[2]=1.-lambda_out[0]; // (autre phase)
      
    }
    else if( zone==23 ){
      //printf("DlambdaB= %g \n", (Ttau*(lambda[0])+Xfm[0])/(Ttau+1));
      lambda_out[0]=0;
      lambda_out[1]=((dt/Ttau)*Xfm[1]+lambda_in[1])/((dt/Ttau)+1.);
      lambda_out[2]=1.-lambda_out[1]; // (autre phase)
      
    }
  }
  else{printf("Erreur CIN_PHASE= %d\n",CIN_PHASE ); return 1;}


  //printf(" ZONE= %d | In : xA= %g , xB= %g, xC= %g   Out : xA= %g , xB= %g, xC= %g\n",zone,lambda_in[0],lambda_in[1],lambda_in[2],lambda_out[0],lambda_out[1],lambda_out[2] );  
  
 

  /*

   Si methode==2 + cinétique, on ne recalcul P,T et C que si on est en zone bi-phasique sinon on resort ici
   ( apres la modif des fractions massiques )

   if(methode==2){return 0;}

   Si la cinétique est activé alors on lance aussi le Newton et le calcul de C

  */



  /*
  // HYDRO
  //printf("   debut hydro  :  zone= %d\n",zone);
  if(zone==0){  
    *pP=PT;
    *pT=TT;
    *pc=0;

    return 0;
  }
  else if(zone==12 || zone==13 || zone==23){
    phaseY=zone%10;
    phaseX=(zone-phaseY)/10;
    // Calcul de p,T et C grace aux nouvelles valeurs de fractions massiques
    if(CIN_PHASE!=0){
      err= VE_mixte_cin_phase(Nmax, epsilon, phaseX, phaseY, lambda_out[phaseY-1], V0x, E0x, V0y, E0y, V, E, &Vx, &Ex, &Vy, &Ey);
      if(err){return 1;}
      *pP=fP(Vx,Ex,phaseX);
      *pT=fT(Vx,Ex,phaseX);
    }
    else if(CIN_PHASE==0){
      *pP=fP(V0x,E0x,phaseX);
      *pT=fT(V0x,E0x,phaseX);
    } 
  }
  else if(zone==1 || zone==2 || zone==3){
    *pP=fP(V,E,zone);
    *pT=fT(V,E,zone);
  }

  /*
  if(*pP>1e11 || *pT<-1e3){
    printf(" Vx= %g , Ex= %g , Vy= %g , Ey= %g\n",Vx,Ex,Vy,Ey);
    printf(" V= %g , E= %g , zone= %d, xA= %g , xB= %g , xC= %g\n",V,E,zone,lambda[0],lambda[01],lambda[2]);
    printf(" Px= %g , Tx= %g , Py= %g , Ty= %g \n",fP(Vx,Ex,phaseX),fT(Vx,Ex,phaseX),fP(Vy,Ey,phaseY),fT(Vy,Ey,phaseY));
    return 1;
  }
  */
  
  /*
  // Problème de cérité du son 
  // On ne peut calculer P de la même manière, il faut utiliser VE_mixte_cin_phase si on est en zone de mélange mixte
  // Faire une nouvelle fonction celerite avec un test dedans
  err=c_celerite_son_cin_phase(Nmax, epsilon, zone, *pP, V, E, dV, V0x, E0x, V0y, E0y, lambda_out, pc);  
  if(err){return 1;}

  return 0;
}
*/


/*
// Script pour tester le Benchmark
int use_fPTC_METH_cin_phase(int methode, int Nv_test, int Ne_test, double Vmin_test, double Emin_test, double Vmax_test, double Emax_test){
  int err;
  double E,V;
  double P,T;
  double x;
  double dt;
  
  // Pour le calcul de la celerité du son
  double V0x, E0x, V0y, E0y;

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  double* lambda=malloc(3*sizeof(double));
  double Ttau;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  // pour pouvoir utiliser la discretisations juste en lisant les fichiers
  //Organisation du fichier :
  //  nb_points
  //  Nv_dis , Ne_dis
  //  Vmin_dis , Vmax_dis
  //  Emin_dis , Emax_dis
  
  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  // Calcul de la célerité du son
  double Pp,Pm, Tp,Tm;
  double dE,c;
  double dV=1e-8;


  double hv_dis=(Vmax_dis-Vmin_dis)/(Nv_dis-1);
  double he_dis=(Emax_dis-Emin_dis)/(Ne_dis-1);

  double hv_test=(Vmax_test-Vmin_test)/(Nv_test-1);
  double he_test=(Emax_test-Emin_test)/(Ne_test-1);

  printf(" hv_test/hv_dis = %lf percent \n",100*hv_test/hv_dis );
  printf(" he_test/he_dis = %lf percent \n",100*he_test/he_dis );

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // fichier binaire    couple (V,E) pour les zones mixte
  FILE *fileVmx, *fileEmx, *fileVmy, *fileEmy;
  if((fileVmx = fopen("fichiers/Vbench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmx = fopen("fichiers/Ebench_mixte_x.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  if((fileVmy = fopen("fichiers/Vbench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier V bench mixte\n"); return 1;}
  if((fileEmy = fopen("fichiers/Ebench_mixte_y.bin", "rb")) == NULL){printf("erreur ouverture fichier E bench mixte\n"); return 1;}
  
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // lecture des tableaux
  fread(tabVx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmx);
  fread(tabEx,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmx);
  fread(tabVy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileVmy);
  fread(tabEy,sizeof(double),(Nv_dis-1)*(Ne_dis-1),fileEmy);

  fclose(fileVmx);  fclose(fileEmx);
  fclose(fileVmy);  fclose(fileEmy);
  
  // Ouverture du fichier ZONE
  FILE *fileZONE_0,*fileZONE_1,*fileZONE_2,*fileZONE_3,*fileZONE_4;
  if((fileZONE_0 = fopen("fichiers/ZONEbench_0.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 0\n"); return 1;}
  if((fileZONE_1 = fopen("fichiers/ZONEbench_1.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 1\n"); return 1;}
  if((fileZONE_2 = fopen("fichiers/ZONEbench_2.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 2\n"); return 1;}
  if((fileZONE_3 = fopen("fichiers/ZONEbench_3.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 3\n"); return 1;}
  if((fileZONE_4 = fopen("fichiers/ZONEbench_4.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE bench 4\n"); return 1;}
  
  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  
  // Ecriture des tableaux
  fread(tabZONE[0],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_0);
  fread(tabZONE[1],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_1);
  fread(tabZONE[2],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_2);
  fread(tabZONE[3],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_3);
  fread(tabZONE[4],sizeof(int),(Nv_dis-1)*(Ne_dis-1),fileZONE_4);

  fclose(fileZONE_0);  fclose(fileZONE_1);  fclose(fileZONE_2);
  fclose(fileZONE_3);  fclose(fileZONE_4);


  // fichier binaire NOMBRE DE ZONE
  FILE *file_Nb_ZONE;
  if((file_Nb_ZONE = fopen("fichiers/Nb_ZONEbench.bin", "rb")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 

  // Ecriture des tableaux
  fread(tab_Nb_ZONE,sizeof(int),(Nv_dis-1)*(Ne_dis-1),file_Nb_ZONE);

  fclose(file_Nb_ZONE);


  // fichier binaire POLYGONES Frontières de zones
  FILE *fileA_V,*fileB_V,*fileC_V,*fileAB_V,*fileAC_V,*fileBC_V;
  FILE *fileA_E,*fileB_E,*fileC_E,*fileAB_E,*fileAC_E,*fileBC_E;
  if((fileA_V = fopen("fichiers/A_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}     if((fileA_E = fopen("fichiers/A_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE A\n"); return 1;}
  if((fileB_V = fopen("fichiers/B_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}     if((fileB_E = fopen("fichiers/B_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE B\n"); return 1;}
  if((fileC_V = fopen("fichiers/C_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}     if((fileC_E = fopen("fichiers/C_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE C\n"); return 1;}
  if((fileAB_V = fopen("fichiers/AB_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}  if((fileAB_E = fopen("fichiers/AB_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AB\n"); return 1;}
  if((fileAC_V = fopen("fichiers/AC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}  if((fileAC_E = fopen("fichiers/AC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE AC\n"); return 1;}
  if((fileBC_V = fopen("fichiers/BC_V.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}  if((fileBC_E = fopen("fichiers/BC_E.bin", "rb")) == NULL){printf("erreur ouverture fichier ZONE BC\n"); return 1;}
  
  double* SA_V=malloc((NA+1)*sizeof(double));    double* SA_E=malloc((NA+1)*sizeof(double));
  double* SB_V=malloc((NB+1)*sizeof(double));    double* SB_E=malloc((NB+1)*sizeof(double));
  double* SC_V=malloc((NC+1)*sizeof(double));    double* SC_E=malloc((NC+1)*sizeof(double));
  double* SAB_V=malloc((NAB+1)*sizeof(double));  double* SAB_E=malloc((NAB+1)*sizeof(double));
  double* SAC_V=malloc((NAC+1)*sizeof(double));  double* SAC_E=malloc((NAC+1)*sizeof(double));
  double* SBC_V=malloc((NBC+1)*sizeof(double));  double* SBC_E=malloc((NBC+1)*sizeof(double));

  // Ecriture des tableaux
  fread(SA_V,sizeof(double),NA+1,fileA_V);     fread(SA_E,sizeof(double),NA+1,fileA_E);
  fread(SB_V,sizeof(double),NB+1,fileB_V);     fread(SB_E,sizeof(double),NB+1,fileB_E);
  fread(SC_V,sizeof(double),NC+1,fileC_V);     fread(SC_E,sizeof(double),NC+1,fileC_E);
  fread(SAB_V,sizeof(double),NAB+1,fileAB_V);  fread(SAB_E,sizeof(double),NAB+1,fileAB_E);
  fread(SAC_V,sizeof(double),NAC+1,fileAC_V);  fread(SAC_E,sizeof(double),NAC+1,fileAC_E);
  fread(SBC_V,sizeof(double),NBC+1,fileBC_V);  fread(SBC_E,sizeof(double),NBC+1,fileBC_E);

  for(int i=0; i<NA+1; i++){  SA[i][0]=SA_V[i];  SA[i][1]=SA_E[i];  }
  for(int i=0; i<NB+1; i++){  SB[i][0]=SB_V[i];  SB[i][1]=SB_E[i];  }
  for(int i=0; i<NC+1; i++){  SC[i][0]=SC_V[i];  SC[i][1]=SC_E[i];  }

  for(int i=0; i<NAB+1; i++){  SAB[i][0]=SAB_V[i];  SAB[i][1]=SAB_E[i];  }
  for(int i=0; i<NAC+1; i++){  SAC[i][0]=SAC_V[i];  SAC[i][1]=SAC_E[i];  }
  for(int i=0; i<NBC+1; i++){  SBC[i][0]=SBC_V[i];  SBC[i][1]=SBC_E[i];  }

  fclose(fileA_V);   fclose(fileB_V);   fclose(fileC_V);
  fclose(fileA_E);   fclose(fileB_E);   fclose(fileC_E);
  fclose(fileAB_V);  fclose(fileAC_V);  fclose(fileBC_V);
  fclose(fileAB_E);  fclose(fileAC_E);  fclose(fileBC_E);

  free(SA_V); free(SB_V); free(SC_V); free(SAB_V); free(SAC_V); free(SBC_V);
  free(SA_E); free(SB_E); free(SC_E); free(SAB_E); free(SAC_E); free(SBC_E);





 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){   fprintf(fSA,"%.10lf %.10lf\n", SA[i][0], SA[i][1] );  }
  for(int i=0; i<NB+1; i++){   fprintf(fSB,"%.10lf %.10lf\n", SB[i][0], SB[i][1] );  }
  for(int i=0; i<NC+1; i++){   fprintf(fSC,"%.10lf %.10lf\n", SC[i][0], SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10lf %.10lf\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );  }

  fclose(fSA);   fclose(fSB);   fclose(fSC);  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////





  //  ****************************************
  // Pour tester les valeurs
  FILE *fPtest, *fTtest, *fXAtest, *fXBtest, *fXCtest, *fCtest;
  if((fPtest = fopen("fichiers/Ptest.txt", "w+")) == NULL){printf("erreur ouverture fichier P test\n");return 1;}
  if((fTtest = fopen("fichiers/Ttest.txt", "w+")) == NULL){printf("erreur ouverture fichier T test\n");return 1;}
  if((fXAtest = fopen("fichiers/XAtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXBtest = fopen("fichiers/XBtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fXCtest = fopen("fichiers/XCtest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  if((fCtest = fopen("fichiers/Ctest.txt", "w+")) == NULL){printf("erreur ouverture fichier X test\n");return 1;}
  FILE *fVtest, *fEtest;
  if((fVtest = fopen("fichiers/Vtest.txt", "w+")) == NULL){printf("erreur ouverture fichier V test \n");return 1;}
  if((fEtest = fopen("fichiers/Etest.txt", "w+")) == NULL){printf("erreur ouverture fichier E test \n");return 1;}
  for(int j=0; j<Ne_test; j++){
      E=Emin_test+j*he_test;
      fprintf(fEtest, "%.15lf ",E);
  }
  for(int j=0; j<Nv_test; j++){
      V=Vmin_test+j*hv_test;
      fprintf(fVtest, "%.15lf ",V);
  }
  fclose(fVtest);  fclose(fEtest);

  int zone;
  FILE *fZONEtest;
  if((fZONEtest = fopen("fichiers/ZONEtest.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE test\n");return 1;}
  //  ****************************************
  
  //  $$$$$$$   Ecriture texte pour verif   $$$$$$$$$$
  FILE *file_Nb_ZONE_test;
  if((file_Nb_ZONE_test = fopen("fichiers/Nb_ZONEbench.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1)*(Ne_dis-1); i++){
    fprintf(file_Nb_ZONE_test, "%d \n",tab_Nb_ZONE[i] );
  }
  fclose(file_Nb_ZONE_test);

  FILE *file_V_dis;
  if((file_V_dis = fopen("fichiers/Vbench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Nv_dis-1); i++){
    fprintf(file_V_dis, "%.15lf \n",Vmin_dis+i*hv_dis+hv_dis/2 );
  }
  fclose(file_V_dis);
  
  FILE *file_E_dis;
  if((file_E_dis = fopen("fichiers/Ebench_dis.txt", "w+")) == NULL){printf("erreur ouverture fichier Nb ZONE bench\n"); return 1;}
  for(int i=0; i<(Ne_dis-1); i++){
    fprintf(file_E_dis, "%.15lf \n",Emin_dis+i*he_dis+he_dis/2 );
  }
  fclose(file_E_dis);
  //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  double Ntau=1, dx=1;


  // ***** Boucle *****
  for(int i=0; i<Nv_test; i++){
    printf("  - i= %d\n",i );
    V=Vmin_test+i*hv_test;
    for(int j=0; j<Ne_test; j++){
      //printf(" - i= %d, j= %d \n",i,j );
      //printf("   V= %g  E= %g\n",V,E );
      E=Emin_test+j*he_test;
      //err=lec_benchmark_PT(1,V, E, nb_points, SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE, 
      //                     Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, &P, &T, &x, &zone, &V0x, &E0x, &V0y, &E0y);
  
      err=fPTC_METH_cin_phase(methode, 1, 0, V, E, &P, &c, &T, lambda, lambda, Ntau, dx, dt, Ttau, &zone,
                         nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,
                         SA, SB, SC, SAB, SAC, SBC, tabVx, tabEx, tabVy, tabEy, tabZONE, tab_Nb_ZONE );
      if(err){return err;}
      
      //

      fprintf(fPtest,"%.10lf\n",P);
      fprintf(fTtest,"%.10lf\n",T);
      fprintf(fXAtest,"%.10lf\n",lambda[0]);
      fprintf(fXBtest,"%.10lf\n",lambda[1]);
      fprintf(fXCtest,"%.10lf\n",lambda[2]);
      fprintf(fZONEtest,"%d\n",zone);
      fprintf(fCtest,"%.10lf\n",c);      
    }
  }

  fclose(fPtest);  fclose(fTtest);  fclose(fXAtest); fclose(fXBtest); fclose(fXCtest); fclose(fZONEtest); fclose(fCtest);
  

  // Liberation de la mémoire
  // Liberation
  free(tabVx);   free(tabEx);
  free(tabVy);   free(tabEy);
  
  freetab(tabZONE);  free(tab_Nb_ZONE);

  freetab(SA);  freetab(SB);  freetab(SC);
  freetab(SAB);  freetab(SAC);  freetab(SBC);
  
  return 0;
}
*/
