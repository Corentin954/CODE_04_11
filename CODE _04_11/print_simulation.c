#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int print_simulation(int tst, int sch, int so, double CFL, int nx, int aff, int Nmax, int kin_fix, int dpi, 
                      int q, double Cq, double Cl, int cin_phase, int Ntau, double dt_first, double Ttau, 
                      int sch_cin_phase, int tst_multi, int EOS, double epsilon){

int ind_tst=tst/100;
printf("Cas-test : ");
if(ind_tst==0){
  if(tst==0){ printf(" Sod  (gaz parfait)\n"); }
  if(tst==1){ printf(" LeBlanc   (gaz parfait)\n"); }
  if(tst==2){ printf(" Bizarrium\n"); }
  if(tst==3){ printf(" Onde acoustique   (gaz parfait)\n"); }
  if(tst==4){ printf(" Sod symétrisé   (gaz parfait)\n"); }
  if(tst==5){ printf(" Woodward (3 états)   (gaz parfait)\n"); }
  if(tst==6){ printf(" Shu-Oscher (2 états avec sinusoide) (gaz parfait)\n"); }
}
else if (ind_tst==1){
  printf("EOS Etain (multiphase):\n"); 
  if(tst==100){ printf("             Projection d'une plaque d'un milimetre contre un mur a gauche à la vitesse u \n");
                printf("             u=320.756,   x=0    beta/gamma     POLE 1\n"); }
  else {
    printf("            Projection de 2 plaques d'un milimetre l'une contre l'autre à la vitesse u \n");
    printf("            Etat initial aux conditions atmosphèriques et P=0\n"); }
    if(tst==101){ printf("             u=320.756, x=0 (beta/gamma) POLE 1\n"); }
    else if(tst==102){ printf("             u=331.046, x=0.2 (beta/gamma)\n");}
    else if(tst==103){ printf("             u=346.306, x=0.5 (beta/gamma)\n");}
    else if(tst==104){ printf("             u=361.359, x=0.8 (beta/gamma)\n");}
    else if(tst==105){ printf("             u=371.266, x=1 (beta/gamma) POLE 2\n"); }
    else if(tst==106){ printf("             u=474.570, phase gamma POLE 3\n"); }
    else if(tst==107){ printf("             u=1333.77, x=0 (gamma/liquide) POLE 4\n"); }
    else if(tst==108){ printf("             u=1383.13, x=0.2 (gamma/liquide) \n"); }
    else if(tst==109){ printf("             u=1453.67, x=0.5 (gamma/liquide)\n"); }
    else if(tst==110){ printf("             u=1520.7, x=0.8 (gamma/liquide)\n"); }
    else if(tst==111){ printf("             u=1563.72, x=1 (gamma/liquide) POLE 5\n"); }
    else if(tst==112){ printf("             u=372.266, (phase gamma)   \n"); }
    else if(tst==112){ printf("             u=1700, (phase liquide)   \n"); }
}
else if(ind_tst==2){
  printf("Etat initial hors équilibre thermo et vitesse nulle \n");
  if(tst==200){ printf(" Etain - Etat initial aux conditions du POLE 4  \n");
                printf(" Air - Etat initial aux conditions atmosphèriques \n"); }
  if(tst==201){ printf(" Etain - Etat initial aux conditions du POLE 3  \n");
                printf(" Air - Etat initial aux conditions atmosphèriques \n"); }
  if(tst==201){ printf(" Etain - Etat initial aux conditions d'un point d'Hugoniot entre P3 et P4  \n");
                printf(" Air - Etat initial aux conditions atmosphèriques \n"); }
}


if(tst_multi==3){
  printf("Initialisation des fluides sur les mailles :");
  printf(" 1/6 du domaine de chaque coté contient de l'air, le reste est de l'étain   \n");
}
else if(tst_multi==4){
  printf("Initialisation des fluides sur les mailles :");
  printf(" 1/6 du domaine à droite contient de l'air, le reste est de l'étain   \n");
}

int ind_sch1=sch/100;
int ind_sch2=(sch-ind_sch1*100)/10;
printf("Schéma hydro : ");
if(ind_sch1==0){
  printf("Godunov energie totale | ");
  if(sch>=0 && sch<=13){ printf("Solveur ordre 1 | "); }
  if(sch==0){ printf("Després\n"); }
  else if(sch==1){ printf(" Solveur ordre 1 Jaouen\n"); }
  else if(sch==2){ printf(" Solveur acoustique ordre 1 à deux états\n"); }
  else if(sch==3){  printf(" RK1 Matsuno + Solveur acoustique ordre 1 à deux états\n"); }
  else if(sch==4){  printf(" Solveur ordre 1 | Extension faiblement non-linéaire à l'ordre 2 du solveur de Riemann acoustique lagrangien \n"); }
  else if(sch==5){  printf(" Solveur ordre 1 | Extension faiblement non-linéaire à l'ordre 3 du solveur de Riemann acoustique lagrangien \n"); }
  
  else if(sch==6){ printf(" Solveur ordre 1 Dukovicz en choc faible\n"); }
  else if(sch==7){ printf(" Solveur ordre 1 Dukovicz en choc fort\n"); }

  else if(sch==8){ printf(" Solveur ordre 1 | Solveur exacte en double choc\n"); }
  else if(sch==9){ printf(" Solveur ordre 1 | Solveur exacte en double detente\n"); }
  else if(sch==10){ printf(" Solveur ordre 1 | Solveur exacte\n"); }

  else if(sch==11){ printf(" Solveur ordre 1 Gallice 1\n"); }
  else if(sch==12){ printf(" Solveur ordre 1 Gallice 2\n"); }
  else if(sch==13){ printf(" Solveur ordre 1 Gallice 3\n"); }

  else if(sch==30 || sch==31 || sch==32){
    printf("Solveur acoustique GAD ordre 2 |");
    if(sch==30){ printf("sans limiteur\n"); }
    else if(sch==31){ printf("limiteur MinMod\n"); }
    else if(sch==32){printf("limiteur Superbee\n"); }
  }

  else if(sch==33 || sch==34 || sch==35){
    printf("Solveur exact GAD ordre 2 | ");
    if(sch==33){ printf("sans limiteur\n"); }
    else if(sch==34){ printf("limiteur MinMod\n"); }
    else if(sch==35){printf("limiteur Superbee\n"); }
  }

  else if(sch==50 || sch==51 || sch==52){
    printf("Solveur GoHy ordre 2 en temps | ");
    if(sch==50){printf("sans limiteur\n"); }
    else if(sch==51){printf("limiteur MinMod\n"); }
    else if(sch==52){printf("limiteur van Leer\n"); }
  }

  else if(sch==53 || sch==54 || sch==55){
    printf("Solveur GoHy ordre 3 en temps | ");
    if(sch==53){ printf("sans limiteur\n"); }
    else if(sch==54){printf("limiteur MinMod\n"); }
    else if(sch==55){printf("limiteur van Leer\n"); }
  }
}

else if(ind_sch1==1){
  //printf("Runge-Kutta energie interne | ");
  if(sch==100){ printf("RK1 (Matsuno forward-backward)\n"); }
  else if(sch==101){  printf("RK2 (Heun's method)\n"); }
  else if(sch==102){  printf("RK3 SSP\n"); }
  else if(sch==103){  printf("RK4 Classique\n"); }
  else if(sch==104){  printf("RK5 Cash-Karp\n"); }
  else if(sch==105){  printf("RK5 Dormand-Prince\n"); }
  else if(sch==106){  printf("RK3 HEUN\n"); }
  else if(sch==107){  printf("RK3 RALSTON\n"); }
  else if(sch==108){  printf("RK3 Bogaki-Shampine\n"); }
  else if(sch==109){  printf("RK2 Bogaki-Shampine\n"); }
  else if(sch==110){  printf("KUTTA ordre 3\n"); }
  else if(sch==111){  printf("RK2 Explicit Midpoint Method\n"); }
  else if(sch==112){  printf("RK2 Ralston\n"); }
  else if(sch==113){  printf("EULER forward\n"); }
  else if(sch==114){  printf("SSPRK4(5,4) Spiteri-Ruuth\n"); }
  else if(sch==115){  printf("SSPRK3(4,3) Spiteri-Ruuth\n"); }
  else if(sch==116){  printf("SSPRK3(5,3) Spiteri-Ruuth\n"); }
  else if(sch==117){  printf("SSPRK1(3,1) Spiteri-Ruuth\n"); }
  else if(sch==118){  printf("SSPRK2(3,2) Spiteri-Ruuth\n"); }
  else if(sch==119){  printf("SSPRK2(4,2) Spiteri-Ruuth\n"); }
  else if(sch==120){  printf("SSPRK1(2,1) Spiteri-Ruuth\n"); }
  else if(sch==121){  printf("SSPRK3(4,3) Spiteri-Ruuth\n"); }
}
else if(ind_sch1==2){
    printf(" BBC | ");
    if(sch==200){ printf("JCP 2009\n"); }
    else if(sch==200){ printf("Predictor-Corrector\n"); }
    else if(sch==200){printf("RK2 average\n"); }
}
else if(ind_sch1==3){  printf(" von Neumann-Richtmyer  \n"); }
else if(ind_sch1==4){  
  printf(" One-step staggered-grid scheme  "); 
  if(sch==400){ printf(" first order\n"); }
  if(sch==401){ printf(" second order\n"); }
}

if(ind_sch1==1 || ind_sch1==2 || ind_sch1==3){
    int ind_q1=q/100;
    int ind_q2=(q-ind_q1*100)/10;
    int ind_q3=(q-ind_q1*100-ind_q2*10);
    printf("Pseudo-viscosité : ");
    if(ind_q1==0){ printf("Ordre 1 usuel | "); }
    else if(ind_q1==1){printf("Ordre 2 Christensen | ");
      if(ind_q2==0){printf("limiteur MinMod | ");}
      else if(ind_q2==1){ printf("limiteur Superbee |"); }
      else if(ind_q2==2){  printf("limiteur Christensen |"); }
    }
    if(ind_q3==0){ printf("von Neumann-Richtmyer\n"); }
    else if(ind_q3==1){ printf("Rosenbluth\n"); }
    else if(ind_q3==2){printf("Landschoff\n"); }
    else if(ind_q3==3){printf("Magical q combination\n"); }
    else if(ind_q3==4){printf("Quadratique + lineaire\n"); }
    else if(q==005){printf("off\n"); }
}

if(ind_sch1==1 || sch==013 || sch==014 || sch==015){
    if(so>=2){
      printf("Ordre en espace : ");
      if(so==2){printf("2nd order\n"); }
      else if(so==3){ printf("3rd order\n"); }
      else if(so==5){ printf("4th and 5th order\n"); }
      else if(so==7){ printf("6th and 7th order\n"); }
      else if(so==9){ printf("8th and 9th order\n"); }
      else{ printf("Erreur so=%d (print_simulation)\n",so ); return 1; }
    }
}

if(ind_sch1==1){
    if(kin_fix>=0){
      printf("Kinetic energy fix : ");
      if(kin_fix==0){ printf("off\n"); }
      else if(kin_fix==1){ printf("CRAS 2016"); }
      else if(kin_fix==2){  printf("JCP 2019"); }
      else if(kin_fix==11){  printf("CRAS 2016 (a chaque sous-pas)"); }
      else if(kin_fix==12){  printf("JCP 2019 (a chaque sous-pas)"); }
      
      if(kin_fix==1 || kin_fix==2 || kin_fix==11 || kin_fix==12){
        printf(" | Cl= %g Cq= %g\n, ",Cl,Cq);
      }
    }

    if(dpi>=0){
      printf("Variante d'ordre élevé énergie interne | ");
      if(dpi==0){printf("dPI (option par default)\n"); }
      else if(dpi==1){printf("AV (variante delta PI bar)\n"); }
    }
}


if(tst>=100){
    if(cin_phase>=0){
      if(cin_phase==0){ printf("Equilibre thermodynamique"); }
      else if(cin_phase==1){ printf("Cinétique de changement de phases avec Ttau= %g ",Ttau); }
      else if(cin_phase==2){printf("Cinétique de changement de phases avec Ttau=Ntau.dX/c avec Ntau= %d ", Ntau); }
      else if(cin_phase==3){printf("Cinétique de changement de phases avec equation d'évolution de P : dP/dt + rho.c2.div(u) = (1/Ttau)(P - Peq) (schéma Euler implcite) avec Ttau= %g", Ttau); }
    }
    printf("   (Critère sur le Newton à l'éqthermo= %g)\n", epsilon);

    printf("Méthode d'implémentation de l'EOS multiphase de l'étain : ");
    if(EOS==0){ printf("Hybride + Newton pour déterminer les lignes de changement de phase\n");  }
    else if(EOS==1){ printf("Hybride + détection des zones\n"); }
    else if(EOS==2){ printf("Tabulation pure   ! Pas implémenté !\n");  }
    else if(EOS==3){ printf("Newton pure   ! Pas implémenté !\n"); }
    else{ printf("Erreur EOS= %d\n", EOS); return 1;}

    

}

printf("Nombre de mailles : %d  ||  ", nx);
printf("CFL= %g   \n", CFL);
printf("max cycle : %d   #  ",Nmax);
printf("fac*dt0 : %g \n ", dt_first);
printf(" -------------------------\n");


return 0; 

}


