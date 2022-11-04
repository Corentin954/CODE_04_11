#!/bin/bash

declare -i i j
declare -a liste_schemes liste_name
declare -i nsimu
declare -i scheme test nx cinetique nmax affichage_cycle viscosity dpi spatial_order kinetic_fix prise_erreur
declare cfl dt0 crit_newton Cq Cl

declare zeta

#tst="0" cas_test="Sod"
#tst="2" cas_test="Bizarium"
tst="102" cas_test="Etain (beta=0.8 gamma=0.2)"
#tst="114" cas_test="Etain choc faible (u0=200 m/s)"

#nx="200"
cfl="0.8"
dt0="1e-3"
nmax="245" 
affichage_cycle="0"
prise_erreur="0"

# EOS Etain
crit_newton="1e-12"
cinetique="0"

liste_nx=(500s00)


#liste_name=(
# 'sol riemann acous.'
# 'KOVA ordre 1 zeta=2'
# 'KOVA ordre 1 zeta=3'
#)

#liste_name=(
# 'solveur acoustique'
# 'solveur exten. ordre 3'
# 'BBC JCP 2009'
#)

liste_name=(
  "solveur acoustique nx=${liste_nx[0]}"
)


liste_schemes=(2)
#liste_schemes=(2 400 400)

 
# Cauchy-Kovaleskaya
#liste_zeta=(1.0 2.0 3.0 4.0)
zeta="-99"

viscosity="4"
Cq="1.5"
Cl="0.15"
dpi="0"
spatial_order="2"
kinetic_fix="1"

nsimu=${#liste_schemes[@]}

#TITLE="Cas test de Bizarrium, 1000 mailles, CFL=0.3"
#TITLE="Cas test ${cas_test}, ${nx} mailles, CFL=${cfl}"
TITLE="Cas test ${cas_test}, CFL=${cfl}"
#TITLE="Cas test choc etain, 400 mailles, CFL=0.3"


# Ecriture d'un fichier contenant les données de la simualtion pour creer les PDF
cat <<EOF > "param_sim.txt"
	${TITLE}
	${tst}
	${nsimu}
EOF
for (( i=0 ; i<nsimu ; i++ )) do # boucle sur les simulations
	echo "${liste_name[i]}" >> "param_sim.txt"
done

# -p : affiche le tableau
#declare -p liste_schemes
#declare -p liste_name
#declare -p liste_zeta

# -- executer le fichier donnant la solution exacte ou de référence

# Calcul de la solution de référence
gcc sol_reference.c func_sol_reference.c -lm
if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || 100 -lt ${tst} ]]; then
  ./a.out ${tst}
fi

#echo ${liste_zeta[0]} ${liste_zeta[1]} ${liste_zeta[2]}

#gcc main_multi.c EOS_etain.c EOS.c func.c cal_front.c allocfree.c funcGtot_multi.c funcRKint_multi.c funcBBC_JCP2009_multi.c funcBBC_PRED_CORR_multi.c funcBBC_RK2av_multi.c funcvNR_multi.c funcGodGoHy_multi.c funcKOVAo2.c print_simulation.c -lgsl -lgslcblas -lm
#printf "!!! compil gcc !!! \n"
sh compil_main_multi.sh
 
printf "nsimu=%d" ${nsimu}

# ------------------   Boucle sur les simulations   ---------------------------------------------
for (( i=0 ; i<nsimu ; i++ ))
do
	printf "\ni=%d Schéma= %d, TITLE=%s\n" ${i} ${liste_schemes[i]} "${liste_name[i]}"
  scheme=liste_schemes[i];
  #zeta="${liste_zeta[i]}";

	printf "prise_erreur=%d (test_tab.sh)\n" ${prise_erreur}

  nx=liste_nx[i];

	#printf "\n   	./a.out  -a %d -c %d -n %d -t %d -s %d -m %lf -w %d -x %lf -j %lf -q %d -d %lf -l %lf -p %d -o %d -z %d -b %lf -g %d (test_tab.sh)\n" ${affichage_cycle} ${nmax} ${nx} ${tst} ${scheme} ${cfl} ${cinetique} ${dt0} ${crit_newton} ${viscosity} ${Cq} ${Cl} ${dpi} ${spatial_order} ${kinetic_fix} ${zeta} ${prise_erreur}

#	./a.out  -a ${affichage_cycle} -c ${nmax} -n ${nx} -t ${tst} -s ${scheme} -m ${cfl} -w ${cinetique} -x ${dt0} -j ${crit_newton} -q ${viscosity} -d ${Cq} -l ${Cl} -p ${dpi} -o ${spatial_order} -z ${kinetic_fix} -b ${zeta} -g ${prise_erreur};

  printf "nombre de thread(s) : %d (test_tab.sh)\n" $OMP_NUM_THREADS
  #ccc_mprun -c 32 -x -p rome-bxi -n 1 .
  ./a.out  -a ${affichage_cycle} -c ${nmax} -n ${nx} -t ${tst} -s ${scheme} -m ${cfl} -w ${cinetique} -x ${dt0} -j ${crit_newton} -q ${viscosity} -d ${Cq} -l ${Cl} -p ${dpi} -o ${spatial_order} -z ${kinetic_fix} -b ${zeta} -g ${prise_erreur};

  cp output_c.txt fichiers/GNUPLOT/data$i
  cp output_d.txt fichiers/GNUPLOT/data_u$i
  cp output_lambda.txt fichiers/GNUPLOT/data_lambda$i
  if [[ ${prise_erreur} -eq 1 ]]; then
    cp errEF.txt fichiers/GNUPLOT/errEF$i
  	cp errEP.txt fichiers/GNUPLOT/errEP$i
  	cp err_ligne_abs.txt fichiers/GNUPLOT/err_ligne_abs$i
  	cp err_ligne_rel.txt fichiers/GNUPLOT/err_ligne_rel$i
  fi
  cp nb_iter.txt fichiers/GNUPLOT/nb_iter$i
done
 
rm output_c.txt output_d.txt output_lambda.txt errEF.txt errEP.txt nb_iter.txt
######################################################################################################


#   Organisation des fichies de sorties du code
#            output_c.txt                                              output_d.txt
#
#  C1  : Xc  maillage centrée                                   C1  : Xd  maillage décalée ou Xc maille centrée (Godunov) 
#  C2  : rho  volume spécifique                                 C2  : u   vitesse 
#  C3  : tau  volume spécifique
#  C4  : e   energie totale
#  C5  : eps   energie interne
#  C6  : P   pression
#  C7  : T   temperature
#  C8  : C   celerite du son
#  C9  : G   derivee fondamentale
#  C10 : S   entropie


#************************#************************#************************#*******************
#************************ Creation des scripts PLOT  **********************************************
bash createPLOT.sh



#************************#************************#************************#*******************
#************************ Creation des scripts PDF  **********************************************
#bash createPDF.sh


#************************#************************#************************#*******************
#************************ Creation des diagrammes erreurs  **********************************************
if [[ ${prise_erreur} -eq 1 ]]; then
  bash createERROR.sh
fi


# option -persist  # la fenetre n'est plus interactive !!?
# taper q lorsque qu'on est sur la fenetre pour la supprimer



#  Explication de la creation des plots

# Pour chaque schéma, on enregistre les données de la ième simuation dans le fichier texte : fichiers/data${i} (data1, data2, ... )
# Seule la vitesse est dans un autre fichier data_u${i} car c'est la seule variable qui peut-être frontières aux mailles

# Pour tracer tourtes les varialbles centrées, on utilise 2 boucles imbriquées :

# La première est sur les variables 
# on crée un liste avec les noms et leurs positions dans les fichiers data${i}

# La seconde est sur les simulations : 
# on crée un fichier texte pour l'executer avec gnuplot dans lequel on trace toutes les simulations sur une même figure 
# pour cela on vient chercher les donées dans data${i} 
# on utilise  la focntion "cat" et "echo" pour ecrire dans le fichier "variables.plot"

# Toutes les données sont enregistrées dans le sous-répertoire /fichiers

