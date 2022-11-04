#!/bin/bash

# rendre un fichier executable sous linux : chmod +x nom_fichier

# Ouvrir un fichier texte avec Sublime text :  subl nom_fichier

gcc main_multi.c EOS_etain.c EOS.c func.c cal_front.c allocfree.c funcGtot_multi.c funcRKint_multi.c funcBBC_JCP2009_multi.c funcBBC_PRED_CORR_multi.c funcBBC_RK2av_multi.c funcvNR_multi.c funcGodGoHy_multi.c print_simulation.c -lgsl -lgslcblas -lm 

tst="1"
nx="2000"
cfl="0.7"
cinetique="0"
dt0="1e-1"
nmax="100000"
crit_newton="1e-12"
affichage_cycle="0"

TITLE="Cas test de Lebalnc, 2000 mailles, CFL=0.3"

scheme="2"
titre1="acoustic"
./a.out  -a ${affichage_cycle} -c ${nmax} -n ${nx} -t ${tst} -s ${scheme} -m ${cfl} -w ${cinetique} -x ${dt0} -j ${crit_newton};
cp output_c.txt data1
#cp output_d.txt acoustic_d.txt
rm output_c.txt 

scheme="11"
titre2="Gallice 1"
./a.out  -a ${affichage_cycle} -c ${nmax} -n ${nx} -t ${tst} -s ${scheme} -m ${cfl} -w ${cinetique} -x ${dt0} -j ${crit_newton};
cp output_c.txt data2
#cp output_d.txt acoustic_d.txt
rm output_c.txt

scheme="12"
titre3="Gallice 2"
./a.out  -a ${affichage_cycle} -c ${nmax} -n ${nx} -t ${tst} -s ${scheme} -m ${cfl} -w ${cinetique} -x ${dt0} -j ${crit_newton}
cp output_c.txt data3
cp output_d.txt data3_u
rm output_c.txt output_d.txt

#            output_c.txt                                              output_d.txt
#
#  C1  : Xc  maillage centrée                                   C1  : Xd  maillage décalée
#  C2  : rho  volume spécifique                                 C2  : u   vitesse (Godunov)  rien (autres types de schémas)
#  C3  : tau  volume spécifique
#  C4  : u   vitesse (Godunov) 0 (autres types de schémas)
#  C5  : e   energie totale
#  C6  : eps   energie interne
#  C7  : P   pression
#  C8  : T   temperature
#  C9  : C   celerite du son
#  C10 : G   derivee fondamentale
#  C11 : S   entropie
rho="2"
tau="3"
u="4"
e="5"
e_interne="6"
P="7"
T="8"
C="9"
G="10"
S="11"

cat <<EOF > density.plot
set term wxt 1 enhanced font "Arial,20"
set title "{/=20 Big} Medium {/=5 Small}"
#set title "${TITLE}"
set xlabel "X"
set ylabel "Density \rho"
set grid
plot "data1" using 1:${rho} w lp title "${titre1}", "data2" using 1:${rho} w lp title "${titre2}" , "data3" using 1:${rho} w lp title "${titre3}"
EOF

cat <<EOF > tau.plot
set terminal x11 2
set title "Etain test problem: 200 lagrangian cells, cfl = 0.95"
set xlabel "X"
set ylabel "Specific volume \tau"
set grid
plot "data1" using 1:${tau} w lp title "${titre1}", "data2" using 1:${tau} w lp title "${titre2}", "data3" using 1:${tau} w lp title "${titre3}"
EOF

cat <<EOF > u.plot
set terminal x11 3
set title "${TITLE}"
set xlabel "X"
set ylabel "Vitesse u"
set grid
plot "data1" using 1:${u_god} w lp title "${titre1}", "data2" using 1:${u_god} w lp "title ${titre2}", "data3" w lp title "${titre3}"
EOF

cat <<EOF > e.plot
set terminal x11
set title "${TITLE}"
set xlabel "X"
set ylabel "Vitesse u"
set grid
plot "data1" using 1:${e} w lp title "${titre1}", "data2" using 1:${e} w lp title "${titre2}", "data3" using 1:${e} w lp title "${titre3}"
EOF

cat <<EOF > eps.plot
set terminal x11
set title "${TITLE}"
set xlabel "X"
set ylabel "Vitesse u"
set grid
plot "data1" using 1:${e_interne} w lp title "${titre1}", "data2" using 1:${e_interne} w lp title "${titre2}", "data3" using 1:${e_interne} w lp title "${titre3}"
EOF

#gnuplot < cmd1.plot
#gnuplot < cmd2.plot


gnuplot "density.plot" -persist  # la fenetre n'est plus interactive !!? 
gnuplot "tau.plot" -persist  # la fenetre n'est plus interactive !!? 

# taper q lorsque qu'on est sur la fenetre pour la supprimer