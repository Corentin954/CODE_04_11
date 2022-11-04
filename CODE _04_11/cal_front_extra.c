

//  ****** Fonction pour imprimer dans des fichiers les points des frontières de domaine P,T  et   V,E
/*
   Utile dans la fonction diag_phase  qui calcul ces points
*/
int print_double(int Np12, double* segP12, double* segT12, double* segVa12, double* segEa12, double* segVb12, double* segEb12,
                 int Np13, double* segP13, double* segT13, double* segVa13, double* segEa13, double* segVb13, double* segEb13,
                 int Np23, double* segP23, double* segT23, double* segVa23, double* segEa23, double* segVb23, double* segEb23,
                 int Nb, int Nh, int Ng, int Nd, double Tb, double Th, double Pg, double Pd,
                 double* segPBAS, double* segPHAUT, double* segTGAUCHE, double* segTDROITE){
  
  FILE* fileP12, *fileT12, *fileVa12, *fileEa12, *fileVb12, *fileEb12;
  FILE* fileP13, *fileT13, *fileVa13, *fileEa13, *fileVb13, *fileEb13;
  FILE* fileP23, *fileT23, *fileVa23, *fileEa23, *fileVb23, *fileEb23;
  FILE* filePb, *filePh, *fileTg, *fileTd;


  if((fileP12 = fopen("fichiers/P12.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT12 = fopen("fichiers/T12.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa12 = fopen("fichiers/Va12.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa12 = fopen("fichiers/Ea12.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb12 = fopen("fichiers/Vb12.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb12 = fopen("fichiers/Eb12.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((fileP13 = fopen("fichiers/P13.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT13 = fopen("fichiers/T13.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa13 = fopen("fichiers/Va13.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa13 = fopen("fichiers/Ea13.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb13 = fopen("fichiers/Vb13.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb13 = fopen("fichiers/Eb13.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((fileP23 = fopen("fichiers/P23.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileT23 = fopen("fichiers/T23.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileVa23 = fopen("fichiers/Va23.txt", "w+")) == NULL){printf("erreur ouverture fichier Va\n");return 1;}
  if((fileEa23 = fopen("fichiers/Ea23.txt", "w+")) == NULL){printf("erreur ouverture fichier Ea\n");return 1;}
  if((fileVb23 = fopen("fichiers/Vb23.txt", "w+")) == NULL){printf("erreur ouverture fichier Vb\n");return 1;}
  if((fileEb23 = fopen("fichiers/Eb23.txt", "w+")) == NULL){printf("erreur ouverture fichier Eb\n");return 1;}

  if((filePb = fopen("fichiers/fbas.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((filePh = fopen("fichiers/fhaut.txt", "w+")) == NULL){printf("erreur ouverture fichier P\n");return 1;}
  if((fileTg = fopen("fichiers/fgauche.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}
  if((fileTd = fopen("fichiers/fdroite.txt", "w+")) == NULL){printf("erreur ouverture fichier T\n");return 1;}


  for(int i=0; i<Np12; i++){
    fprintf(fileP12, "%.15lf ",segP12[i] );
    fprintf(fileT12, "%.15lf ",segT12[i] );
    fprintf(fileVa12, "%.15lf ",segVa12[i] );
    fprintf(fileEa12, "%.15lf ",segEa12[i] );
    fprintf(fileVb12, "%.15lf ",segVb12[i] );
    fprintf(fileEb12, "%.15lf ",segEb12[i] );  
  }

  for(int i=0; i<Np13; i++){
    fprintf(fileP13, "%.15lf ",segP13[i] );
    fprintf(fileT13, "%.15lf ",segT13[i] );
    fprintf(fileVa13, "%.15lf ",segVa13[i] );
    fprintf(fileEa13, "%.15lf ",segEa13[i] );
    fprintf(fileVb13, "%.15lf ",segVb13[i] );
    fprintf(fileEb13, "%.15lf ",segEb13[i] );  
  }

  for(int i=0; i<Np23; i++){
    fprintf(fileP23, "%.15lf ",segP23[i] );
    fprintf(fileT23, "%.15lf ",segT23[i] );
    fprintf(fileVa23, "%.15lf ",segVa23[i] );
    fprintf(fileEa23, "%.15lf ",segEa23[i] );
    fprintf(fileVb23, "%.15lf ",segVb23[i] );
    fprintf(fileEb23, "%.15lf ",segEb23[i] );  
  }

  for(int i=0; i<Nb; i++){
    fprintf(filePb, "%15.lf %.15lf\n",segPBAS[i],Tb );
  }
  for(int i=0; i<Nh; i++){
    fprintf(filePh, "%15.lf %.15lf\n",segPHAUT[i], Th );
  }
  for(int i=0; i<Ng; i++){
    fprintf(fileTg, "%15.lf %.15lf\n",Pg, segTGAUCHE[i] );
  }
  for(int i=0; i<Nd; i++){
    fprintf(fileTd, "%15.lf %.15lf\n",Pd, segTDROITE[i] );
  }



  fclose(fileP12);   fclose(fileT12);   
  fclose(fileVa12);  fclose(fileEa12);   
  fclose(fileVb12);  fclose(fileEb12);   

  fclose(fileP13);   fclose(fileT13);   
  fclose(fileVa13);  fclose(fileEa13);   
  fclose(fileVb13);  fclose(fileEb13);   

  fclose(fileP23);   fclose(fileT23);   
  fclose(fileVa23);  fclose(fileEa23);   
  fclose(fileVb23);  fclose(fileEb23);   

}


// *********    Calcul du Point courbe frontières   ************
/*
   Dans le domaine (P,T)    ET    dans le domaine (V,E)
*/
int diag_phase(void){ 

  int Nmax=1e4;
  double epsilon=5e-5;


  int nb_points=20;
  int Np12=2*nb_points; int Np13=nb_points; int Np23=6*nb_points;
  int Nh=3*nb_points;  
  int Ng1=nb_points;    int Ng2=4*nb_points; 
  int Nd1=4*nb_points;  int Nd2=0.5*nb_points; 
  int Nb1=2*nb_points;  int Nb2=4*nb_points; 
  int Ng=Ng1+Ng2;
  int Nd=Nd1+Nd2;
  int Nb=Nb1+Nb2;


  double Ptriple=4.5104217816091e9;
  double Pdeb=Ptriple;


  int phaseA12=1;
  int phaseB12=2;
  double Pfin12=9.4e9;
  
  int phaseA13=1;
  int phaseB13=3;
  double Pfin13=0.;
  
  int phaseA23=2;
  int phaseB23=3;
  double Pfin23=70e9;
  
  double Va012, Vb012, Va013, Vb013, Va023, Vb023;
  double Ea012, Eb012, Ea013, Eb013, Ea023, Eb023;


  // Valeur d'initialisation
  Va012=129.905671875844291e-6;  Ea012=72.168190127265447e3; 
  Vb012=127.666112185063994e-6;  Eb012=99.480256119236515e3;
  
  Va013=129.905671875844291e-6;  Ea013=72.168190127265447e3; 
  Vb013=132.538478805442338e-6;  Eb013=136.154346228414312e3;

  Va023=127.666112185063994e-6;  Ea023=99.480256119236515e3; 
  Vb023=132.538478805442338e-6;  Eb023=136.154346228414312e3;


   
  // Initialisation des tableaux
  double* segP12=malloc(Np12*sizeof(double));
  double* segT12=malloc(Np12*sizeof(double));
  double* segVa12=malloc(Np12*sizeof(double));
  double* segEa12=malloc(Np12*sizeof(double));
  double* segVb12=malloc(Np12*sizeof(double));
  double* segEb12=malloc(Np12*sizeof(double));

  double* segP13=malloc(Np13*sizeof(double));
  double* segT13=malloc(Np13*sizeof(double));
  double* segVa13=malloc(Np13*sizeof(double));
  double* segEa13=malloc(Np13*sizeof(double));
  double* segVc13=malloc(Np13*sizeof(double));
  double* segEc13=malloc(Np13*sizeof(double));

  double* segP23=malloc(Np23*sizeof(double));
  double* segT23=malloc(Np23*sizeof(double));
  double* segVb23=malloc(Np23*sizeof(double));
  double* segEb23=malloc(Np23*sizeof(double));
  double* segVc23=malloc(Np23*sizeof(double));
  double* segEc23=malloc(Np23*sizeof(double));

  double* segPHAUT=malloc(Nh*sizeof(double));
  double* segPBAS=malloc((Nb1+Nb2)*sizeof(double));
  double* segTGAUCHE=malloc((Ng1+Ng2)*sizeof(double));
  double* segTDROITE=malloc((Nd1+Nd2)*sizeof(double));

  // frontiere 12
  //printf("PHASE 12 :\n");
  ligne_double(Nmax, epsilon, phaseA12, phaseB12, Pdeb, Pfin12, Np12, Va012, Ea012, Vb012, Eb012, segP12, segT12, segVa12, segEa12, segVb12, segEb12);

  // frontiere 13
  //printf("PHASE 13 :\n");
  ligne_double(Nmax, epsilon, phaseA13, phaseB13, Pdeb, Pfin13, Np13, Va013, Ea013, Vb013, Eb013, segP13, segT13, segVa13, segEa13, segVc13, segEc13);
  
  // frontiere 23
  //printf("PHASE 23 :\n");
  ligne_double(Nmax, epsilon, phaseA23, phaseB23, Pdeb, Pfin23, Np23, Va023, Ea023, Vb023, Eb023, segP23, segT23, segVb23, segEb23, segVc23, segEc23);
  
  
  //frontières extérieures :
  double Pg=0.;
  double Pd=70e9;
  double Tb=300;
  double Th=2500;

  double Td=2486.138962087622531;
  double Tg=505;
  double Pb=9.4e9;
  // HAUT
  for(int i=0; i<Nh; i++){
    segPHAUT[i]=Pg+i*(Pd-Pg)/(Nh-1);
  }

  // BAS
  for(int i=0; i<Nb1; i++){
    segPBAS[i]=Pg+i*(Pb-Pg)/(Nb1-1);
  }
  for(int i=0; i<Nb2; i++){
    segPBAS[Nb1+i]=Pb+(i+1)*(Pd-Pb)/(Nb2);
  }

  // GAUCHE
  for(int i=0; i<Ng1; i++){
    segTGAUCHE[i]=Tb+i*(Tg-Tb)/(Ng1-1);
  }
  for(int i=0; i<Ng2; i++){
    segTGAUCHE[Ng1+i]=Tg+(i+1)*(Th-Tg)/(Ng2);
  }

  // DROITE
  for(int i=0; i<Nd1; i++){
    segTDROITE[i]=Tb+i*(Td-Tb)/(Nd1-1);
  }
  for(int i=0; i<Nd2; i++){
    segTDROITE[i+Nd1]=Td+(i+1)*(Th-Td)/(Nd2);
  }


  /*    fonction print   */
  int res= print_double(Np12, segP12, segT12, segVa12, segEa12, segVb12, segEb12,  
                        Np13, segP13, segT13, segVa13, segEa13, segVc13, segEc13,    
                        Np23, segP23, segT23, segVb23, segEb23, segVc23, segEc23,
                        Nb, Nh, Ng, Nd, Tb, Th, Pg, Pd,
                        segPBAS, segPHAUT, segTGAUCHE, segTDROITE);    
  if(res){printf("Erreur dans l'impression (diag_phase)"); return 1;}



  int err;
  int phase;
  double P;
  double T;
  double E,V;
  
  FILE* fres1;
  FILE* fres2;
  FILE* fres3;
  if((fres1 = fopen("fichiers/VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("fichiers/VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("fichiers/VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}

  

  //  ### Calcul et impression des frontières dans le domaine V,E
  int N;
  double Tdeb, Tfin;
  double Pfin;

  // frontière 1/A
  phase=1;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=300;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  T=300;
  Pdeb=0;
  Pfin=9.4e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres1);


  // frontière 2/B
  phase=2;
  N=100;
  T=300;
  Pdeb=9.4e9;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  P=70e9;
  Tdeb=300.;
  Tfin=2486.138962087622531;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres2);


  // frontière 3/C
  phase=3;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=2500;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  T=2500;
  Pdeb=0.;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres3);
 

  // Liberation
  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);

  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);

  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);

  free(segPHAUT);    free(segPBAS);
  free(segTGAUCHE);  free(segTDROITE);

}


// Calcul de VE
/*
   Discretise les frontière du diagramme de phase (P,T)
   1 fichier pour chaque phase : "VEfront1.txt"  ,  "VEfront2.txt"  et  "VEfront3.txt"
*/
int routine_VE(void){
  int err;
  int phase;
  double P;
  double T;
  double E,V;

  
  FILE* fres1;
  FILE* fres2;
  FILE* fres3;
  if((fres1 = fopen("fichiers/VEfront1.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres2 = fopen("fichiers/VEfront2.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}
  if((fres3 = fopen("fichiers/VEfront3.txt", "w+")) == NULL){printf("erreur ouverture fichier (routine_VE)\n");return 1;}


  int N;
  double Tdeb, Tfin;
  double Pdeb, Pfin;

  // fronitère de la phase 1/A
  phase=1;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=300;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  T=300;
  Pdeb=0;
  Pfin=9.4e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres1, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres1);

  //  ******************************************************

  // fronitère de la phase 2/B
  phase=2;
  N=100;
  T=300;
  Pdeb=9.4e9;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  P=70e9;
  Tdeb=300.;
  Tfin=2486.138962087622531;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.5lf %.5lf\n",V*1e6,E );
    fprintf(fres2, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres2);


  //  ******************************************************

  // fronitère de la phase 3/C
  phase=3;
  N=100;
  P=0;
  Tdeb=505;
  Tfin=2500;
  for(int i=0; i<N; i++){
    T=Tdeb + i*(Tfin-Tdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  T=2500;
  Pdeb=0.;
  Pfin=70e9;
  for(int i=0; i<N; i++){
    P=Pdeb + i*(Pfin-Pdeb)/(N-1);
    err= fVE(phase, P, T, &V, &E);
    //err= fVE_NR(phase, P, T, &V, &E);
    //printf("%.15lf %.15lf\n",V*1e6,E );
    fprintf(fres3, "%.15lf %.15lf\n",V,E );
  }
  fclose(fres3);

}


// fonction permettant de connaitre le domaine de validité de la phase liquide 
/*    -> calcul des grandeurs thermo une large domaine rectangulaire
         lors de la visualisation de ces grandeurs avec Octave, 
         on peut voir où est-ce que les grandeurs ne sont pas définies
*/
int plot_PT_zone3_domaine_validite(){
  int res; 

  int NP,NT;
  NP=1e2; NT=5e2;
  double Pdeb=-15e9, Pfin=80e9;
  double Tdeb=0, Tfin=25e3;
  double hP=(Pfin-Pdeb)/(NP-1);
  double hT=(Tfin-Tdeb)/(NT-1);

  double P,T,V,E;
  
  FILE* fP3, *fT3, *fV3, *fE3;
  if((fP3 = fopen("fichiers/P_3.txt", "w+")) == NULL){printf("erreur ouverture fichier \n");return 1;}
  if((fT3 = fopen("fichiers/T_3.txt", "w+")) == NULL){printf("erreur ouverture fichier \n");return 1;}
  if((fV3 = fopen("fichiers/V_3.txt", "w+")) == NULL){printf("erreur ouverture fichier \n");return 1;}
  if((fE3 = fopen("fichiers/E_3.txt", "w+")) == NULL){printf("erreur ouverture fichier \n");return 1;}

  for(int iP=0; iP<NP; iP++){
    P=Pdeb+iP*hP;  
    for(int iT=0; iT<NT; iT++){  
      T=Tdeb+iT*hT;
      res=fVE(3, P, T, &V, &E);
      P=fP(V,E,3);
      T=fT(V,E,3);
      fprintf(fP3, "%g ",P );
      fprintf(fT3, "%g ",T );
      fprintf(fV3, "%g ",V );
      fprintf(fE3, "%g ",E );
    }
    fprintf(fP3, "\n ");
    fprintf(fT3, "\n ");
    fprintf(fV3, "\n ");
    fprintf(fE3, "\n ");
  }
  
  fclose(fP3); fclose(fT3); 
  fclose(fV3); fclose(fE3); 
}




/* Fonction qui calcul les grandeurs thermo sur un quadrillage de l'espace
     problème : le calcul des gr,ageurs thermo ne se fait grace à la fonction EOS presente dans les schémas
       cad ne permet de verifier le bon codage de la fonction EOS qu idansle code

*/
int mapsPT(void){
  
  int Nmax=1e4;
  double epsilon=1e-6;
  double Vx, Ex, Vy, Ey, x;
  double P,T;

  double Vat,Vbt,Vct;
  double Eat,Ebt,Ect;
  double PT, TT;
  point_triple(&PT, &TT, &Vat, &Eat, &Vbt, &Ebt, &Vct, &Ect);
  

  int err,res;
  int Nv=1e2;
  int Ne=1e2;

  double Vmin, Vmax;
  double Emin, Emax;
  
  Vmin=80e-6, Vmax=210e-6;
  Emin=-2e3, Emax=16e5;
  
  /*
  Vmin=110e-6, Vmax=150e-6;
  Emin=-5e3, Emax=5e5; 
  */

  double E,V;


  // nb_points est l'argument d'entree de la fonction
  int nb_points=20;
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
  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, 
                       segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );

  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}



  ///////////
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  
  for(int i=0; i<NA+1; i++){
    fprintf(fSA,"%.10lf %.10lf\n",SA[i][0],SA[i][1] );
  }

  for(int i=0; i<NB+1; i++){
    fprintf(fSB," %.10lf %.10lf\n",SB[i][0],SB[i][1] );
  }

  for(int i=0; i<NC+1; i++){
    fprintf(fSC,"%.10lf %.10lf\n",SC[i][0],SC[i][1] );
  }

  for(int i=0; i<NAB+1; i++){
    fprintf(fSAB,"%.10lf %.10lf\n",SAB[i][0],SAB[i][1] );
  }

  for(int i=0; i<NAC+1; i++){
    fprintf(fSAC," %.10lf %.10lf\n",SAC[i][0],SAC[i][1] );
  }

  for(int i=0; i<NBC+1; i++){
    fprintf(fSBC,"%.10lf %.10lf\n",SBC[i][0],SBC[i][1] );
  }
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);
 
  /////////////////

  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);


  FILE *fmaps, *fmapsGNU, *fE, *fV;
  if((fmaps = fopen("fichiers/MAPS.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fmapsGNU = fopen("fichiers/P_GNU.txt", "w+")) == NULL){printf("erreur ouverture fichier MAPS\n");return 1;}
  if((fV = fopen("fichiers/Vmaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Vmaps \n");return 1;}
  if((fE = fopen("fichiers/Emaps.txt", "w+")) == NULL){printf("erreur ouverture fichier Emaps \n");return 1;}


  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fE, "%.15g ",E);
  }
  
  double V_unuse_x, E_unuse_x;
  double V_unuse_y, E_unuse_y;

  for(int i=0; i<Nv; i++){
    V=Vmin+i*hv;
    fprintf(fV, "%.15g ", V);
    for(int j=0; j<Ne; j++){
      //printf("  -j= %d\n",j );
      E=Emin+j*he;
      //printf("V= %.15g, E=%.15g\n",V*1e6,E );
      if(intPOLY(V,E,NC,SC)){ 
        P=fP(V,E,3);
        T=fT(V,E,3);        
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(intPOLY(V,E,NA,SA)){
        P=fP(V,E,1);
        T=fT(V,E,1);
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(intPOLY(V,E,NB,SB)){
        P=fP(V,E,2);
        T=fT(V,E,2);
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(ftri(V, E)){
        P=PT;
        T=TT;
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(intPOLY(V,E,NAB,SAB)){
        /*
        err=VE_mixte(Nmax, epsilon, 1, 2, V, E, &Vx, &Ex, &Vy, &Ey, &x);
        P=fP(Vx,Ex,1);
        T=fT(Vx,Ex,1);
        */
        err=interp_zone_mixte(V, E, 1, 2, N12, segVa12, segVb12, segEa12, segEb12, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
        //printf("P= %g, T= %g   mixte 1/2\n",P*1e-9,T );
        fprintf(fmaps, "%.15g %.15g\n", P,T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(intPOLY(V,E,NAC,SAC)){
        /*
        err=VE_mixte(Nmax, epsilon, 1, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
        P=fP(Vx,Ex,1);
        T=fT(Vx,Ex,1);
        */
        err=interp_zone_mixte(V, E, 1, 3, N13, segVa13, segVc13, segEa13, segEc13, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
        //printf("P= %g, T= %g   mixte 1/3\n",P*1e-9,T );
        fprintf(fmaps, "%.15g %.15g\n", P,T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else if(intPOLY(V,E,NBC,SBC)){
        /*
        err=VE_mixte(Nmax, epsilon, 2, 3, V, E, &Vx, &Ex, &Vy, &Ey, &x);
        P=fP(Vx,Ex,2);
        T=fT(Vx,Ex,2);
        */
        err=interp_zone_mixte(V, E, 2, 3, N23, segVb23, segVc23, segEb23, segEc23, &P, &T, &V_unuse_x, &E_unuse_x, &V_unuse_y, &E_unuse_y);
        //printf("P= %lf, T= %lf   mixte 2/3\n",P*1e-9,T );
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
      else{
        P=-10e9;
        T=0;
        fprintf(fmaps, "%.15g %.15g\n", P, T);
        fprintf(fmapsGNU, "%.15g %.15g %.15g\n", V, E, P);
      }
    }
  }
  

  fclose(fmaps);
  fclose(fV);
  fclose(fE);
  fclose(fmapsGNU);

  // Liberation
  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);


  return 0;

}



// Renvoie la valeur de P et T à partir du benchmark et du newton sur les zones mixtes
/*  Arguments d'entrée :
     - int affichage : affichage==1 -> sort P=-10e9 et T=0 lorsque on est hors domaine ~~ PERMET L'AFFICHAGE DE P ET T ~~
                       afichage==0 --> renvoie un erreur lorsque on est hors domaine ~~ PERMET L'UTILISATION DE LA FONCTION DANS LES SCHEMAS ~~
     - V,E
     - nb_points : nombre de points de discretisation des frontières
     - SA,SB,SC etc.. contiennent les discretisations des frontières
     - tabVx, tabEx, ect.. valeur d'initialisation pour le newton sur les zones mixtes
     - tab_Nb_ZONE : nombre de zones sur chaque case
     - tabZONE : case(s) présente(nt) sur chaque case
     - Nv,Ne, Vmin,Vmax,etc ... paramétrent du cadrillage
    Arguments de sorties :
     - pP,pT,px,pzone : pointeur sur la P,T,fraction massique et la zone de (V,E)
     - pV0x, pE0x, etc...  pointeur sur la valeur initiale choisi pour lancer une zone mixte (On peut y avoir acces à la fonction avant grace à tabVx)
*/
int lec_benchmark_PT(int affichage, double V, double E, int nb_points, double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                     double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE, 
                     int Nv, int Ne, double Vmin, double Vmax, double Emin, double Emax, double* pP, double* pT, double* px, int* pzone, 
                     double* pV0x, double* pE0x, double* pV0y, double* pE0y){
  int Nmax=1e2;
  double epsilon=1e-8;

  double critere;
  int err, n_iter; 
  double P,T;
  double Vmx, Emx, Vmy, Emy;
  double Vx, Ex, Vy, Ey;

  double x;

  // nb_points est un argument de la fonction
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  // Point triple
  double VaT,VbT,VcT;
  double EaT,EbT,EcT;
  double PT,TT;
  point_triple(&PT, &TT, &VaT, &EaT, &VbT, &EbT, &VcT, &EcT);
 // Valeur d'initialisation
  /*  ancien jeu de coeff
  double VaT=129.905671875844291e-6;  double EaT=72.168190127265447e3;
  double VbT=127.666112185063994e-6;  double EbT=99.480256119236515e3;
  double VcT=132.538478805442338e-6;  double EcT=136.154346228414312e3;
  */

  if(V<Vmin || Vmax<V){printf("V hors domaine : Vmin= %.15lf  Vmax= %.15lf  et V= %.15lf (lec_benchmark_PT)\n",Vmin,Vmax,V);   return 1;}  
  if(E<Emin || Emax<E){printf("E hors domaine : Emin= %.15lf  Emax= %.15lf  et E= %.15lf (lec_benchmark_PT)\n",Emin,Emax,E);   return 1;}
 
  double hv=(Vmax-Vmin)/(Nv-1);
  double he=(Emax-Emin)/(Ne-1);
  int i=(V-Vmin)/hv;
  int j=(E-Emin)/he;


 //  printf("  Vrec(i)= %.15lf  Vrec(i+1)= %.15lf\n",Vmin+i*hv, Vmin+(i+1)*hv );
 //  printf("  Erec(j)= %.15lf  Erec(j+1)= %.15lf\n",Emin+j*he, Emin+(j+1)*he );


  int pos_fic=i*(Ne-1) + j;
  //printf("  --i_dis= %d, j_dis= %d \n",i,j );


  Vmx=tabVx[pos_fic];  Emx=tabEx[pos_fic];
  Vmy=tabVy[pos_fic];  Emy=tabEy[pos_fic];
  //printf("Vmx= %lf , Emx= %lf \n",Vmx,Emx );
  //printf("Vmy= %lf , Emy= %lf \n",Vmy,Emy );


  *pV0x=Vmx; *pE0x=Emx;
  *pV0y=Vmy; *pE0y=Emy;


  
  int nb_zone=tab_Nb_ZONE[pos_fic];
  
  /*
  printf("  nb zone= %d\n",nb_zone );
  for(int i=0; i<nb_zone; i++){
    printf("  tabZONE[%d][pos_fic]= %d",i,tabZONE[i][pos_fic]);
  }
  printf("\n");
  */
  

  int zone;
  int phaseA, phaseB;

  
  // On regarde si lorsque qu'il ya 3 zones ou plus, on a -1 dans les zones
  // Dans ce cas on mets comme graines segV[N] 
  // sinon c'est qu'on est au abord du point triple et dans ce cas on met comme graines les sommets du point triple 
  int isout=0;
  if(nb_zone>=3){  
    for(int i=0; i<nb_zone; i++){
      if(tabZONE[i][pos_fic]==-1){isout=1;}
    }
  }
  // on se sert du isout dans la boucle ou nb_zone>=3



  //printf(" nb_zone= %d\n",nb_zone );
  
  //  ##  1 seule zone ##
  if(nb_zone==1){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );

    if(zone==-1){ 
     if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
     else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
     else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else if(zone==3){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
    else if(zone==2){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
    else if(zone==1){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
    else if(zone==23){
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);  if(err){return 1;}
      *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;
      return 0;
    }
    else if(zone==13){
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);  if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);    *pzone=13;
      return 0;
    }
    else if(zone==12){
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);  if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;
      return 0;
    }
    else if(zone==0){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0; }
  }
  //  ##  2 zones   ##   la première zone est forcément la zone pur (ou triple)
  else if(nb_zone==2){
    zone=tabZONE[0][pos_fic];
    //printf("   - zone= %d\n",zone );
    if(zone==3){
      if(intPOLY(V,E,NC,SC)){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone==2){
      if(intPOLY(V,E,NB,SB)){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;   *pzone=2;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;

        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA); *pzone=zone;
        
        return 0;
      } 
    }
    else if(zone==1){
      if(intPOLY(V,E,NA,SA)){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
      else if(tabZONE[1][pos_fic]==-1){  
        if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone==0){ 
      if(ftri(V, E)){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }
      else{
        zone=tabZONE[1][pos_fic];
        phaseB=zone%10;  phaseA=(zone-phaseB)/10;
 
        err=VE_mixte(Nmax, epsilon, phaseA, phaseB, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
        if(err){return 1;}
        *pP=fP(Vx,Ex,phaseA);  *pT=fT(Vx,Ex,phaseA);  *pzone=zone;

        return 0;
      } 
    }
    else if(zone=-1){  
      if(affichage==1){*pP=-10e9; *pT=0.; *px=-1; *pzone=-1; return 0;}
      else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else{printf(" Erreur valeur de zone= %d\n",zone );}
  }
  //  ##  3 zones ou plus  ##
  else{
    for(int ind=0; ind<nb_zone-1; ind++){
      zone=tabZONE[ind][pos_fic];
      //printf("   - zone= %d\n",zone );
      if(zone==3){
        if(intPOLY(V,E,NC,SC)){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
      }
      else if(zone==2){
        if(intPOLY(V,E,NB,SB)){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
      }
      else if(zone==1){
        if(intPOLY(V,E,NA,SA)){  *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
      }
      else if(zone==23){
        if(intPOLY(V,E,NBC,SBC)){
          if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
          else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
          if(err){return 1;}
          *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;

          return 0;
        }
      }
      else if(zone==13){
        if(intPOLY(V,E,NAC,SAC)){
          if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
          err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
          if(err){return 1;}
          *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1); *pzone=13;

          return 0;
        }
      }  
      else if(zone==12){
        if(intPOLY(V,E,NAB,SAB)){
          if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
          else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
          err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
          if(err){return 1;}
          *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;

          return 0;
        }
      }
      else if(zone==0){   if(ftri(V, E)){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }   }  
      else if(zone=-1){  
        if(affichage==1){*pP=-10e9;  *pT=0.;  *px=-1;  *pzone=-1;  return 0;}
        else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E);  return 1;}
        else{printf("Erreur dans la valeur de affichage= %d",affichage);  return 1;}
      }
    }
    
    zone=tabZONE[nb_zone-1][pos_fic];
    if(zone==3){  *pP=fP(V,E,3);  *pT=fT(V,E,3);  *px=-1;  *pzone=3;  return 0;  }
    else if(zone==2){  *pP=fP(V,E,2);  *pT=fT(V,E,2);  *px=-1;  *pzone=2;  return 0;  }
    else if(zone==1){   *pP=fP(V,E,1);  *pT=fT(V,E,1);  *px=-1;  *pzone=1;  return 0;  }
    else if(zone==23){
      if(isout){Vmx=SBC[N23][0]; Emx=SBC[N23][1]; Vmy=SBC[N23-1][0]; Emy=SBC[23-1][1];}
      else{Vmx=VbT; Emx=EbT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 2, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
      if(err){return 1;}
      *pP=fP(Vx,Ex,2);  *pT=fT(Vx,Ex,2);  *pzone=23;

      return 0;
    }
    else if(zone==13){
      if(isout){Vmx=SAC[N13][0]; Emx=SAC[N13][1]; Vmy=SAC[N13-1][0]; Emy=SAC[N13-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VcT; Emy=EcT;}
      err=VE_mixte(Nmax, epsilon, 1, 3, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
      if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=13;

      return 0;
    }  
    else if(zone==12){
      if(isout){Vmx=SAB[N12][0]; Emx=SAB[N12][1]; Vmy=SAB[N12-1][0]; Emy=SAB[N12-1][1];}
      else{Vmx=VaT; Emx=EaT; Vmy=VbT; Emy=EbT;}
      err=VE_mixte(Nmax, epsilon, 1, 2, Vmx, Emx, Vmy, Emy, V, E, &Vx, &Ex, &Vy, &Ey, px, &n_iter, &critere);
      if(err){return 1;}
      *pP=fP(Vx,Ex,1);  *pT=fT(Vx,Ex,1);  *pzone=12;

      return 0;
    }
    else if(zone==0){  *pP=PT;  *pT=TT;  *px=-1;  *pzone=0;  return 0;  }
    else if(zone=-1){  
      if(affichage==1){*pP=-10e9; *pT=0.;  *px=-1;   *pzone=-1;  return 0;}
      else if(affichage==0){printf("Error out of domain : V= %.15lf , E= %.15lf\n",V,E); return 1;}
      else{printf("Erreur dans la valeur de affichage== %d",affichage); return 1;} 
    }
    else{printf("Erreur dans tabZONE[%d][%d]= %d\n",nb_zone-1,pos_fic,tabZONE[nb_zone-1][pos_fic]); return 1;}
    
    
  }
  // Si on arrive la c'est que la bonne zone n'a pas été trouvé
  printf("Erreur choix zone (lec_benchmark_PT) ");
  printf("  nb_zone= %d",nb_zone );
  printf("  tabZONE[0][%d]= %d",pos_fic,tabZONE[0][pos_fic]);
  for(int ind=1; ind<nb_zone; ind++){
    printf("  tabZONE[%d][%d]= %d",ind,pos_fic,tabZONE[ind][pos_fic]); 
  }
  printf("\n V= %.10lf , E=%.10lf\n",V,E);
  return 1; 
  
}

// ZONE 
// repere les zones (triangle pour le moment)
int diag_zone(){
  int err,res; 
  int Nmax=1e2;
  double epsilon=1e-8;

  double E,V;

  // nb_points est un argument de la fonction
  int nb_points=100;
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3; 
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;
  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );


  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pou l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);

  // Initialisation des tableaux
  double* segP12=malloc(N12*sizeof(double));   double* segT12=malloc(N12*sizeof(double));
  double* segVa12=malloc(N12*sizeof(double));  double* segEa12=malloc(N12*sizeof(double));
  double* segVb12=malloc(N12*sizeof(double));  double* segEb12=malloc(N12*sizeof(double));
  
  double* segP13=malloc(N13*sizeof(double));   double* segT13=malloc(N13*sizeof(double));
  double* segVa13=malloc(N13*sizeof(double));  double* segEa13=malloc(N13*sizeof(double));
  double* segVc13=malloc(N13*sizeof(double));  double* segEc13=malloc(N13*sizeof(double));
  
  double* segP23=malloc(N23*sizeof(double));   double* segT23=malloc(N23*sizeof(double));
  double* segVb23=malloc(N23*sizeof(double));  double* segEb23=malloc(N23*sizeof(double));
  double* segVc23=malloc(N23*sizeof(double));  double* segEc23=malloc(N23*sizeof(double));
  

  err=sommets_polygone(NA, NB, NC, NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3, N12, N13, N23, SA, SB, SC, SAB, SAC, SBC, segP12, segT12, segP13, segT13, segP23, segT23, segVa12, segEa12, segVb12, segEb12, segVa13, segEa13, segVc13, segEc13, segVb23, segEb23, segVc23, segEc23 );
  if(err){printf("err= %d (sommets_polygone)\n",err ); return err;}
  
  /*
  for(int i=0; i<NA+1; i++){
    printf("   SA[%d][0]= %lf  SA[%d][1]= %lf \n",i,SA[i][0],i,SA[i][1] );
  }
  printf("\n");
  for(int i=0; i<NB+1; i++){
    printf("   SB[%d][0]= %lf  SB[%d][1]= %lf \n",i,SB[i][0],i,SB[i][1] );
  }
  printf("\n");
  for(int i=0; i<NC+1; i++){
    printf("   SC[%d][0]= %lf  SC[%d][1]= %lf \n",i,SC[i][0],i,SC[i][1] );
  }
  */
   


  // DIagramme (V,E) 
  // Ecriture des sommets des poly pour vérif
  FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  if((fSAB = fopen("fichiers/pSAB.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAB\n");return 1;}
  if((fSAC = fopen("fichiers/pSAC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSAC\n");return 1;}
  if((fSBC = fopen("fichiers/pSBC.txt", "w+")) == NULL){printf("erreur ouverture fichier pSBC\n");return 1;}
  

  for(int i=0; i<NA+1; i++){  fprintf(fSA,"%.10g %.10g\n",SA[i][0],SA[i][1] );  }
  for(int i=0; i<NB+1; i++){  fprintf(fSB,"%.10g %.10g\n",SB[i][0],SB[i][1] );  }
  for(int i=0; i<NC+1; i++){  fprintf(fSC,"%.10g %.10g\n",SC[i][0],SC[i][1] );  }
  for(int i=0; i<NAB+1; i++){  fprintf(fSAB,"%.10g %.10g\n",SAB[i][0],SAB[i][1] );  }
  for(int i=0; i<NAC+1; i++){  fprintf(fSAC,"%.10g %.10g\n",SAC[i][0],SAC[i][1] );  }
  for(int i=0; i<NBC+1; i++){  fprintf(fSBC,"%.10g %.10g\n",SBC[i][0],SBC[i][1] );  }
  
  
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  fclose(fSAB);  fclose(fSAC);  fclose(fSBC);


  // Diagramme (P,T)
  // Ecriture des sommets des poly pour vérif
  //FILE *fSA, *fSB, *fSC, *fSAB, *fSAC, *fSBC;
  if((fSA = fopen("fichiers/pSA_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSA\n");return 1;}
  if((fSB = fopen("fichiers/pSB_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSB\n");return 1;}
  if((fSC = fopen("fichiers/pSC_PT.txt", "w+")) == NULL){printf("erreur ouverture fichier pSC\n");return 1;}
  

  for(int i=0; i<NA+1; i++){  fprintf(fSA,"%.10lf %.10lf\n",fP(SA[i][0],SA[i][1],1 ),fT(SA[i][0],SA[i][1],1 )) ;  }
  for(int i=0; i<NB+1; i++){  fprintf(fSB,"%.10lf %.10lf\n",fP(SB[i][0],SB[i][1],2 ),fT(SB[i][0],SB[i][1],2 )) ;  }
  for(int i=0; i<NC+1; i++){  fprintf(fSC,"%.10lf %.10lf\n",fP(SC[i][0],SC[i][1],3 ),fT(SC[i][0],SC[i][1],3 )) ;  }
  
  
  fclose(fSA);   fclose(fSB);   fclose(fSC);
  

  //  ****  Permet une répresentation des zones (calcul la zone de chaque points)  **********
  /*
  double Vmin, Vmax;
  double Emin, Emax;

  int Nv=1e3;
  int Ne=1e3;

  Vmin=115e-6, Vmax=1450e-6;
  Emin=-5e3,   Emax=4e5;
  
  double he=(Emax-Emin)/(Ne-1);
  double hv=(Vmax-Vmin)/(Nv-1);
  
  
  FILE* fZONE; FILE* fE; FILE* fV;
  if((fZONE = fopen("ZONE.txt", "w+")) == NULL){printf("erreur ouverture fichier ZONE\n"); return 1;}
  if((fV = fopen("Vzone.txt", "w+")) == NULL){printf("erreur ouverture fichier Vzone \n"); return 1;}
  if((fE = fopen("Ezone.txt", "w+")) == NULL){printf("erreur ouverture fichier Ezone \n"); return 1;}
  
  for(int j=0; j<Ne; j++){
      E=Emin+j*he;
      fprintf(fE, "%.15lf ",E);
  }
  
  
  for(int i=0; i<Nv; i++){
    printf(" -- i= %d\n",i );
    V=Vmin+i*hv;
    fprintf(fV, "%.15lf ", V);
    for(int j=0; j<Ne; j++){
      //printf("  -j= %d\n",j );
      E=Emin+j*he;
      //printf("V= %.15lf, E=%.15lf\n",V*1e6,E );
      if(ftri(V, E)){fprintf(fZONE, "0 ");}
      else if(intPOLY(V,E,NA,SA)){fprintf(fZONE, "1 ");}
      else if(intPOLY(V,E,NB,SB)){fprintf(fZONE, "2 ");}
      else if(intPOLY(V,E,NC,SC)){fprintf(fZONE, "3 ");}
      else if(intPOLY(V,E,NAB,SAB)){fprintf(fZONE, "12 ");}
      else if(intPOLY(V,E,NAC,SAC)){fprintf(fZONE, "13 ");}
      else if(intPOLY(V,E,NBC,SBC)){fprintf(fZONE, "23 ");}
      else{fprintf(fZONE, "-1 ");}
    }
    fprintf(fZONE, "\n");
  }  

  fclose(fZONE);
  fclose(fV);
  fclose(fE);
  */


  // Liberation
  free(segP12);   free(segT12);
  free(segVa12);  free(segEa12);
  free(segVb12);  free(segEb12);
  
  free(segP13);   free(segT13);
  free(segVa13);  free(segEa13);
  free(segVc13);  free(segEc13);
  
  free(segP23);   free(segT23);
  free(segVb23);  free(segEb23);
  free(segVc23);  free(segEc23);
}
