#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include "head_multi.h"

//#include "head.h"


/////////////////////  MAIN  ////////////////////
int main(int argc, char* argv[]){
  double T;
  int nx=200;
  double CFL=0.3; 
  double a, b;
  int Nmax=2e6;
  int tst=0; 
  int sch=0;	// choix du calcul de Ustar et Pstar (0:Solveur acoustique  1:Despres)
  int aff=1;
  int kin_fix=0;
  int q=0;
  int so=2; // spatial order
  int dpi=0; //choix du dPI pour l'intégration de rho.u  
  double Cq=1.5; //--> artificial vicosity q
  double Cl=0.15;
  int iu=0;
  int err=0;
  int ind_cond;
  int Ntau=15;
  double dt_first=1.0;
  double Ttau=1e-9;
  int cin_phase=0;
  int sch_cin_phase=0;
  int tst_multi=0;
  int tst_air=0;
  int EOS=0;
  double epsilon=1e-13;
  double zeta=1.0;
  int calcul_err;
  


	for(int i=0; i<argc; i++){
	  if(argv[i][0] == '-'){
	  	switch(argv[i][1]){
        case 'h' : printf("Code d'étude d'hydrodynamique lagrangienne 1D plane\n"); 
	  	             printf("Usage : \n"); 
    	  	  			 printf(" -t test problem \n");
                   printf("                EOS Gaz parfait ou Bizarrium :\n"); 
                   printf("                    000 : Sod\n");
                   printf("                    001 : LeBlanc\n");
                   printf("                    002 : Bizarrium\n"); 
                   printf("                    003 : Onde acoustique\n");
                   printf("                    004 : Sod symétrisé\n");
                   printf("                    005 : Woodward (3 états)\n");
                   printf("                    006 : Shu-Oscher (2 états avec sinusoide)\n");
                   printf("                EOS Etain (multiphase):\n"); 
                   printf("                  Projection d'une plaque d'un milimetre contre un mur a gauche à la vitesse u \n");
                   printf("                    100 : u=320.756,   x=0    beta/gamma     POLE 1\n");
                   printf("                  Projection de 2 plaques d'un milimetre l'une contre l'autre à la vitesse u \n");
                   printf("                   Etat initial aux conditions atmosphèriques et P=0\n");
                   printf("                    101 : u=320.756,   x=0    beta/gamma      POLE 1\n");
                   printf("                    102 : u=331.046,   x=0.2  beta/gamma\n");
                   printf("                    103 : u=346.306,   x=0.5  beta/gamma\n");
                   printf("                    104 : u=361.359,   x=0.8  beta/gamma\n");
                   printf("                    105 : u=371.266,   x=1    beta/gamma      POLE 2\n");
                   printf("                    106 : u=474.570,   phase gamma            POLE 3\n");
                   printf("                    107 : u=1333.77,   x=0    gamma/liquide   POLE 4\n");
                   printf("                    108 : u=1383.13,   x=0.2  gamma/liquide \n");
                   printf("                    109 : u=1453.67,   x=0.5  gamma/liquide\n");
                   printf("                    110 : u=1520.7,    x=0.8  gamma/liquide\n");
                   printf("                    111 : u=1563.72,   x=1    gamma/liquide   POLE 5\n");
                   printf("                    112 : u=372.266,   phase gamma   \n");
                   printf("                    113 : u=1700,   phase liquide   \n");                  
                   printf("                    114 : u=200, choc faible (permet de rester dans la hase beta sauf le wall-heating)   \n");
                   printf("                   Etat initial hors équilibre thermo et vitesse nulle \n");
                   printf("                    200 : Etain - Etat initial aux conditions du POLE 4  \n");
                   printf("                          Air - Etat initial aux conditions atmosphèriques \n");
                   printf("                    201 : Etain - Etat initial aux conditions du POLE 3  \n");
                   printf("                          Air - Etat initial aux conditions atmosphèriques \n");
                   printf("                    202 : Etain - Etat initial aux conditions d'un point d'Hugoniot entre P3 et P4  \n");
                   printf("                          Air - Etat initial aux conditions atmosphèriques \n");
                   
                   printf(" -r initialisation des fluides sur les mailles \n");
                   printf("                    0 : Gaz parfait uniquement (gamma selon le cas-test)  \n");
                   printf("                    1 : Bizarrium uniquement   \n");
                   printf("                    2 : Etain uniquement  \n");
                   printf("                    3 : 1/6 du domaine de chaque coté contient de l'air, le reste est de l'étain   \n");
                   printf("                    4 : 1/6 du domaine à droite contient de l'air, le reste est de l'étain   \n");
                   
                   printf(" -y cinetic scheme \n");
                   printf("                    0 : ordre 1 full-backward implicite   \n");
                   printf("                    1 : Fickett-Rivard explicite  \n");

    	  	  			 printf(" -s hydro scheme\n"); 
    	  	  			 printf("    	 Godunov energie totale\n");
                   printf("                  Solveur ordre 1\n");
    	  	  			 printf("                    000 : Després\n");
    	  	  			 printf("                    001 : Jaouen\n");
                   printf("                    002 : Solveur acoustique ordre 1 à deux états\n");
                   
                   printf("                    003 : RK1 Matsuno + Solveur acoustique ordre 1 à deux états\n");
                   printf("                    004 : Extension faiblement non-linéaire à l'ordre 2 du solveur de Riemann acoustique lagrangien \n");
                   printf("                    005 : Extension faiblement non-linéaire à l'ordre 3 du solveur de Riemann acoustique lagrangien \n");
                   
                   printf("                    006 : Solveur de Dukovicz en choc faible\n");
                   printf("                    007 : Solveur de Dukovicz en choc fort\n");
                   
                   printf("                    008 : Solveur de Riemann en choc \n");
                   printf("                    009 : Solveur de Riemann en detente \n");
                   printf("                    010 : Solveur de Riemann exact \n");
                   
                   printf("                    011 : Solveur Gallice 1 \n");
                   printf("                    012 : Solveur Gallice 2 \n");
                   printf("                    013 : Solveur Gallice 3 \n");
    	  	  			 
                   printf("                  Solveur acoustique GAD ordre 2 \n");
    	  	  			 printf("                    030 : sans limiteur\n");
    	  	  			 printf("                    031 : limiteur MinMod\n");
                   printf("                    032 : limiteur Superbee\n");
            
                   printf("                  Solveur exacte GAD ordre 2\n");
                   printf("                    033 : sans limiteur\n");
                   printf("                    034 : limiteur MinMod\n");
                   printf("                    035 : limiteur Superbee\n");
                  
                   printf("                  Solveur GoHy ordre 2 en temps\n");
                   printf("                    050 : sans limiteur\n");
                   printf("                    051 : limiteur MinMod\n");
                   printf("                    052 : limiteur van Leer\n");

                   printf("                  Solveur GoHy ordre 3 en temps\n");
                   printf("                    053 : sans limiteur\n");
                   printf("                    054 : limiteur MinMod\n");
    	  	  			 printf("                    055 : limiteur van Leer\n");


    	  	  			 printf("      Runge-Kutta energie interne\n");
    	  	  			 printf("                    100: RK1 (Matsuno forward-backward)\n");
                   printf("                    101: RK2 (Heun's method)\n");
                   printf("                    102: RK3 SSP\n");
                   printf("                    103: RK4 Classique\n");
                   printf("                    104: RK5 Cash-Karp\n");
                   printf("                    105: RK5 Dormand-Prince\n");
                   printf("                    106: RK3 HEUN\n");
                   printf("                    107: RK3 RALSTON\n");
                   printf("                    108: RK3 Bogaki-Shampine\n");
                   printf("                    109: RK2 Bogaki-Shampine\n");
                   printf("                    110: KUTTA ordre 3\n");
                   printf("                    111: RK2 Explicit Midpoint Method\n");
                   printf("                    112: RK2 Ralston\n");
                   printf("                    113: EULER forward\n");
                   printf("                    114: SSPRK4(5,4) Spiteri-Ruuth\n");
                   printf("                    115: SSPRK3(4,3) Spiteri-Ruuth\n");
                   printf("                    116: SSPRK3(5,3) Spiteri-Ruuth\n");
                   printf("                    117: SSPRK1(3,1) Spiteri-Ruuth\n");
                   printf("                    118: SSPRK2(3,2) Spiteri-Ruuth\n");
                   printf("                    119: SSPRK2(4,2) Spiteri-Ruuth\n");
                   printf("                    120: SSPRK1(2,1) Spiteri-Ruuth\n");
                   printf("                    121: SSPRK3(4,3) Spiteri-Ruuth\n");   
                   printf("    	 BBC  \n");
                   printf("                    200: BBC JCP 2009\n");
                   printf("                    201: BBC Predictor-Corrector\n");
                   printf("                    202: BBC RK2 average\n");
                   printf("    	 von Neumann-Richtmyer  \n");
                   printf("                    300: vNR\n");
                   printf("    	 Cauchy-Kovalevskaya staggered grid 2nd order (par defaut zeta=1.0)\n");
                   printf("                    400: ordre 1 \n");
                   printf("                    401: ordre 2 \n");

    	  	  			 printf(" -q pseudo viscosity \n");
                   printf("      Ordre 1 usuel : \n");
                   printf("        000 : von Neumann-Richtmyer\n");
                   printf("        001 : Rosenbluth\n");
                   printf("        002 : Landschoff\n");
                   printf("        003 : Magical q combination\n");
                   printf("        004 : Quadratique + lineaire\n");
                   printf("        005 : off\n");
                   printf("      Ordre 2 Christensen : \n");
                   printf("       limiteur MinMod \n");
                   printf("        100 : von Neumann-Richtmyer\n");
                   printf("        101 : Rosenbluth\n");
                   printf("        102 : Landschoff\n");
                   printf("        103 : Magical q combination\n");
                   printf("        104 : Quadratique + lineaire\n");
                   printf("       limiteur Superbee \n");
                   printf("        110 : von Neumann-Richtmyer\n");
                   printf("        111 : Rosenbluth\n");
                   printf("        112 : Landschoff\n");
                   printf("        113 : Magical q combination\n");
                   printf("        114 : Quadratique + lineaire\n");
                   printf("       limiteur Christensen \n");
                   printf("        120 : von Neumann-Richtmyer\n");
                   printf("        121 : Rosenbluth\n");
                   printf("        122 : Landschoff\n");
                   printf("        123 : Magical q combination\n");
                   printf("        124 : Quadratique + lineaire\n");

      	  		     printf(" -o spatial order \n");
                   printf("        2 : 2nd order\n");
                   printf("        3 : 3rd order\n");
                   printf("        5 : 4th and 5th order\n");
                   printf("        7 : 6th and 7th order\n");
                   printf("        9 : 8th and 9th order\n");

                   printf(" -n number of cells \n");
                   printf("        200 : valeur par défault\n");

                   printf(" -c final cycle \n");
                   printf("        2e6 : valeur par défault\n");

                   printf(" -m CFL \n");
                   printf("        0.3 : valeur par défault\n");

                   printf(" -a affichage des cycles \n"); 
                   printf("        0 : off \n");
                   printf("        1 : on (valeur par default)\n");

                   printf(" -z kinetic energy fix\n");
                   printf("        0 : off\n");
                   printf("        1 : CRAS 2016\n");
                   printf("        2 : JCP 2019\n");
                   printf("      A chaque sous-pas (Runge-Kutta) :\n");
                   printf("        11 : CRAS 2016\n");
                   printf("        21 : JCP 2019\n");

                   printf(" -p variante d'ordre élevé énergie interne\n");
                   printf("        0 : dPI (option par default)\n");
                   printf("        1 : AV (variante delta PI bar)\n");

                   printf(" -d options\n");
                   printf("        valeur du coefficient pseudo-quadratique\n");

                   printf(" -l options\n");
                   printf("        valeur du coefficent pseudo-lineaire\n");

                   printf(" -w cinétique de changement de phase\n");
                   printf("        0 : équilibre thermodynamique  (valeur par défault)\n");
                   printf("        1 : cinétique de changement de phases avec Ttau fixé\n");
                   printf("        2 : cinétique de changement de phases avec Ttau=Ntau.dX/c\n");
                   printf("        3 : cinétique de changement de phases avec equation d'évolution de P : dP/dt + rho.c2.div(u) = (1/Ttau)(P - Peq) (schéma Euler implcite)\n");

                   printf(" -i Ntau entier caractéristique de la cinétique de changement de phase\n");
                   printf("        15 : valeur par défault\n");

                   printf(" -u Ttau\n");
                   printf("        1e-9 : valeur par défault\n");

                   printf(" -x coefficent devant le premier pas de temps\n");
                   printf("        1.0 : valeur par défault\n");

                   printf(" -e choix de la méthode d'implémentation de l'EOS multiphase de l'étain\n");
                   printf("        0 : Hybride + Newton pour déterminer les lignes de changement de phase (valeur par défault)\n");
                   printf("        1 : Hybride + détection des zones\n");
                   printf("        2 : Tabulation pure   ! implémenté mais pas optimisé !\n"); 
                   printf("        3 : Newton pure   ! Pas implémenté !\n");

                   printf(" -j valeur du critère sur le Newton à l'équilibre thermodynamique\n");
                   printf("        1e-13 : valeur par défault\n");

                   printf(" -b valeur de zeta pour le schéma Cauchy-Kovalevskaya staggered grid\n");
                   printf("        1.0 : valeur par défault\n");
    	  	  			 
                   printf(" -g calcul de l'erreur pour certains cas-tests\n");
                   printf("        0 : pas de calcul d'erreur (valeur par défaut)\n");
                   printf("        1 : calcul de l'erreur \n");
                   

                   return 0;
	  	  case 't' : tst=atoi(argv[++i]);
	  	  			     break;
	  	  case 's' : sch=atoi(argv[++i]);
	  	  			     break;
        case 'o' : so=atoi(argv[++i]);
                   break;
        case 'm' : CFL=strtod(argv[++i],NULL);
                   break;
	  	  case 'n' : nx=atoi(argv[++i]);
	  	  			     break;
        case 'a' : aff=atoi(argv[++i]);
                   break;
	  	  case 'c' : Nmax=atoi(argv[++i]);
	  	  			     break;
        case 'z' : kin_fix=atoi(argv[++i]);
                   break;
        case 'p' : dpi=atoi(argv[++i]);
                   break;
        case 'q' : q=atoi(argv[++i]);
                   break;
        case 'd' : Cq=strtod(argv[++i],NULL);
                   break;
        case 'l' : Cl=strtod(argv[++i],NULL);
                   break;
        case 'w' : cin_phase=atoi(argv[++i]);
                   break;
        case 'i' : Ntau=atoi(argv[++i]);
                   break;
        case 'x' : dt_first=strtod(argv[++i],NULL);
                   break;
        case 'u' : Ttau=strtod(argv[++i],NULL);
                   break;    
        case 'y' : sch_cin_phase=atoi(argv[++i]);
                   break;
        case 'r' : tst_multi=atoi(argv[++i]);
                   break;
        case 'e' : EOS=atoi(argv[++i]);
                   break;
        case 'j' : epsilon=strtod(argv[++i],NULL);
                   break;
        case 'b' : zeta=strtod(argv[++i],NULL);
                   break;
        case 'g' : calcul_err=atoi(argv[++i]);
                   break;
	  	  default :  printf("Not recognized option : %c\n",argv[i][1] );
	  	  			     return 1;
	  	}
	  }
	}
  
  // ###########
  
  printf("\n  calcul_err=%d (main_multi) \n\n",calcul_err );
  
  // ################
  
  double t_start,t_end;



  int etain=tst/100;
  //printf("etain= %d\n",etain );

// Réécriture de so 
  int ind_so;
  if(so==2 || so==3){ ind_so=so;}
  else if(so==5){ ind_so=4;}
  else if(so==7){ ind_so=5;}
  else if(so==6){ ind_so=6;}

// Test de la valeur de e-->EOS 
  if(EOS==0 || EOS==1){} //ok
  else if(EOS==2 || EOS==3){printf("Pas implémenté (main)\n");}
  else {printf("Erreur valeur EOS=%d (main)\n",EOS ); return 1;}
  
  if(cin_phase==0 || cin_phase==1 || cin_phase==2 || cin_phase==3){}
  else {printf("erreur cin_phase= %d",cin_phase); return 1;}


// ########   Lecture des fichiers binaires contenants les discretisations des frontières  ############
    
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
  

  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

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
  
  if(etain){
    //printf("Entree lec_tabulation\n");
    err=lec_tabulation(nb_points, Nv_dis, Ne_dis, tabZONE, tab_Nb_ZONE, tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    //if(err){printf("Erreur d'ouvertures des fichiers binaires (lec_tabulation)\n"); return 1;}
  }
    

// ##############################################################################
  
  int deux_plaques=-1; // si choc_double==0 --> choc simple sinon choc double 
  if(tst==100){deux_plaques=0;}
  else if(200>tst && tst>100){deux_plaques=1;}
  
  if(1){ // tst 
    if (tst==0){
      a=0.;  b=1.0;
      T=0.2;
      tst_multi=0;
      tst_air=tst;
    }  
    else if(tst==1){
      a=0.0;  b=9.0;  
      T=6;
      tst_multi=0;
      tst_air=tst;
    }
    else if (tst==2){
      //if(cin_phase==0 || cin_phase==1 || cin_phase==2){ printf("Erreur tst= %d et cin_phase= %d\n",tst,cin_phase ); return 1;}
      a=0.;  b=1.0;  
      T=80e-6; 
      tst_multi=1;
      tst_air=tst;
    }
    else if(tst==3){
      a=0.;  b=1.0;  
      T=1.0;
      tst_multi=0;
      tst_air=tst;
    }
    else if(tst==4){
      a=-1.0;  b=1.0;  
      T=0.2;
      tst_multi=0;
      tst_air=tst;
    }
    else if(tst==5){
      a=0.0;  b=1.0;  
      T=0.038;
      tst_multi=0;
      tst_air=tst;
    }
    else if(tst==6){
      a=-5.0;  b=5.0;  
      T=1.8;
      tst_multi=0;
      tst_air=tst;
    }
    else if (deux_plaques==0){
      a=0.;  b=1.0e-3;
      T=2e-7;
      tst_multi=4;
    }
    else if (deux_plaques==1){
      a=-1e-3;  b=1e-3;
      T=0.2e-6;
      tst_multi=2;
    }
    else if (tst==200){
      a=0;  b=1e-3;
      T=0.1e-6;
      tst_multi=4;
    }
    else if(tst==201){
      a=0; b=1e-3;
      T=0.1e-6;
      tst_multi=4;
    }
    else if(tst==202){
      a=0; b=1e-3;
      T=0.1e-6;
      tst_multi=4;
    }
    else {
      printf("Erreur choix tst= %d (main_etain)\n",tst);
      return 1;
    }
  }


  if(tst_multi==1){
    tst_air=2;
  }
  else if(tst_multi==2){
    tst_air=0; // pas d'importance car n'intervient pas
  }
  else if(tst_multi==3){
    tst_air=0;
  }
  else if(tst_multi==4){
    tst_air=0;
  }

  err=print_simulation(tst, sch, so, CFL, nx, aff, Nmax, kin_fix, dpi, q, Cq, Cl,cin_phase, 
                       Ntau, dt_first, Ttau, sch_cin_phase, tst_multi, EOS, epsilon);
  if(err){ printf("Erreur (print_simulation)\n"); return 1; }
  


  // params pour les cond aux limites
  // changement de signe pour u si iu==1
  // Condition de type mur à gauche et libre à droite

  if(1){ // Condition de type mur
    if (tst==0|| tst==1 || tst==4 || tst==5){ //|| tst==4  
      ind_cond=0;
      iu=1;   // indice pour la vitesse
    }
    // Condition de type ? (libre)
    else if (tst==2|| tst==6 ){
      ind_cond=0;
      iu=0; 
    }
    // Condition de type périodique
    else if (tst==3){
      ind_cond=1;
      iu=0; 
    }
    else if (deux_plaques==0){  
      ind_cond=2;
      iu=1;   // indice pour la vitesse
    }
    else if (deux_plaques){  
      ind_cond=0;
      iu=0;   // indice pour la vitesse
    }
    else if (tst==200){  
      ind_cond=2;
      iu=1;   // indice pour la vitesse
    }
    else {
      printf("Erreur dans le choix de tst= %d (main cond lim)\n",tst);
      return 1;
    }
  }
  

  // fonction clock_t
  clock_t temps_initial, temps_final;
  temps_initial = clock ();
  
  t_start = omp_get_wtime();

  // sch = 100*Rsch + 10*Ksch + ??;
  int Rsch=floor(sch/100);
  sch=sch-100*Rsch;
  int Ksch=floor(sch/10);


  if(Rsch==0){
    if(sch<50){
      err=funcGtot_multi(sch, tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, EOS,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                   tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
    }
    else if(sch>=50){
      if(so!=2 && so!=3 && so!=5 && so!=7 && so!=9){printf("Erreur so=%d (main_multi)\n",so ); return 1;}
      err=funcGodGoHy_multi(sch, tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, ind_so, EOS,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                   tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
      /*
      err=funcGodGoHy3_multi(sch, tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, EOS,
                   nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                   tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
      */

    }
    else{printf("Erreur sch=%d\n",sch ); return 1;}
  }
  else if(Rsch==1){
    if(so!=2 && so!=3 && so!=5 && so!=7 && so!=9){printf("Erreur so=%d\n",so ); return 1;}
    err=funcRKint_multi(sch, tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, q, Cq, Cl, kin_fix, dpi, ind_so, EOS,
                  nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                  tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){return 1;}
  }
  else if(Rsch==2){
    if(sch==0){
      err=funcBBC_JCP2009_multi(tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, q, Cq, Cl, EOS,
                          nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,  tabZONE, tab_Nb_ZONE, 
                          tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
    }

    else if(sch==1){
      err=funcBBC_PRED_CORR_multi(tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, q, Cq, Cl, EOS,
                            nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                            tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
    }

    else if(sch==2){
      err=funcBBC_RK2av_multi(tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, q, Cq, Cl, EOS,
                        nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                        tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){return 1;}
    }
    else{printf("Erreur sch (main)\n");}
  }

  else if(Rsch==3){
    err=funcvNR_multi(tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, q, Cq, Cl, EOS,
                nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){return 1;}
  }

  else if(Rsch==4){
    /*
    err=funcKOVA_multi(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, z,
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                 tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){return 1;}
    */
    err=funcKOVAo2_multi(sch,tst, T, a, b, nx, epsilon, Nmax, CFL, aff, ind_cond, iu, Ntau, dt_first, Ttau, cin_phase, sch_cin_phase, tst_multi, tst_air, EOS, zeta,
                          nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,  tabZONE, tab_Nb_ZONE, 
                          tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){return 1;}
    
  }
  else{printf("erreur choix sch (main)\n"); return 1;}
   
  t_end = omp_get_wtime();
  temps_final = clock ();
  double temps_cpu = (temps_final - temps_initial) * 1e-6;
  int min=(temps_cpu)/60;
  double sec=temps_cpu-60*min;
  //printf("     OLD    Temps de calcul = %g secondes   = %d min %g s \n",temps_cpu, min, sec);
  printf("Temps de calcul= %g secondes (ompenMP)\n",t_end-t_start);
  
  printf(" calcul_err=%d\n", calcul_err);

  if(calcul_err==1){
    printf("Calcul de l'erreur :\n");
    err=sol_err_exacte(tst, nx, Rsch); if(err){return err;}  // from func.c
  }

}
