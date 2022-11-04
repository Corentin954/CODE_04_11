#!/bin/bash

declare -i i j
declare -a liste_schemes 
declare -i simu nsimu nmaille nsch scheme test nx cinetique nmax affichage_cycle viscosity dpi spatial_order kinetic_fix prise_erreur


# On lit le fichier avec les parametres de la simulation
i="-3"
ind_nx=0;
ind_sch=0;
while read line
do
	if [[ $i -eq -3 ]]; then	TITLE=${line}; 	fi
	if [[ $i -eq -2 ]]; then	nsch=${line}; 	fi
	if [[ $i -eq -1 ]]; then 	nmaille=${line}; 	fi
	if [[ -1 -lt $i && $i -lt ${nsch} ]]; then  liste_schemes_name[ind_sch]=${line}; (( ind_sch=${ind_sch}+1 ));	fi
	if [[ ${nsch}-1 -lt $i ]]; then  liste_nx[ind_nx]=${line}; (( ind_nx=${ind_nx}+1 ));	fi
  (( i=${i}+1 ));
done < param_courbe_conv.txt


((nsimu = nmaille * nsch))


declare -a liste_name
for (( i=0 ; i<nsch ; i++ ))
do
  for (( j=0 ; j<nmaille ; j++ ))
  do
    simu=$((nmaille*i + j))
    liste_name[simu]="${liste_schemes_name[i]}, ${liste_nx[j]} mailles"
  done
done

declare -p liste_name 
declare -p liste_schemes_name
declare -p liste_nx 
  
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
#************************ Creation des diagrammes erreurs  **********************************************

#
# The "histograms" style automates the construction of barcharts with
# various stacking and clustering options.  But it doesn't offer much
# choice for coloring.
#
# This demo show using the "with boxes" plot style to generate clustered
# histograms with individual box coloring generated by a user-provided function.
# This could, for example, depend on data from other columns in the input.
#


rho="2"; tau="3"; e="4"; e_interne="5"; P="6"; T="7"; C="8"; G="9"; S="10"; BETA="11"; GAMMA="12"; LIQ="13";
declare -ra liste_plot=(rho tau e e_interne P T C G S BETA GAMMA LIQ)
declare -ra liste_plot_name=("Densité (kg/m)" "Volume spécifique (m/kg)" "Energie totale (J/kg)" "Energie interne (J/kg)" "Pression (Pa)" "Temperature (K)" "Celerite du son (m/s)" "Entropie" "Dérivée fondamentale" "Fraction massique beta" "Fraction massique gamma" "Fraction massique liquide")
declare -i nplot=${#liste_plot[@]}



## --------------   BARPLOT erreur bar de toutes les variables comparées -----------------------------------------------------------
name_fich="fichiers/GNUPLOT/erreur.plot"
cat <<EOF > $name_fich
  set term x11 font "Arial,20"
  set title "${TITLE}"
  set grid
  red = "#FF0000"; green = "#00FF00"; blue = "#0000FF"; skyblue = "#87CEEB";
  set yrange [0:1]
  set style data histogram
  set style histogram cluster gap 1
  set style fill solid
  set boxwidth 0.9
  set grid ytics
  set xtics ('rho' 0, 'tau' 1, 'etot' 2, 'epsilon' 3, 'P' 4, 'T' 5, 'C' 6, 'G' 7, 'S' 8,'beta' 9, 'gamma' 10, 'liq' 11 )
EOF
#  set xtics ('A' 0, 'B' 1, 'C' 2, 'D' 3)
 
echo -n plot >> $name_fich
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'errEF_courbe_conv${i}' using 2 title '${liste_name[i]}'," >> $name_fich
done
echo " 'errEF_courbe_conv${i}' using 2 title '${liste_name[i]}'" >> $name_fich
###############################################################################################




#################################################################################################
## --------------   boucle BARPLOT sur les erreurs variables ABSOLUES  -----------------------------
for (( j=0 ; j<nplot ; j++ ))  # boucle sur les variables
do
  var=${liste_plot[j]}
  name_fich="fichiers/GNUPLOT/err_${var}_abs.plot"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
  set term x11 font "Arial,20"
  set title "${liste_plot[j]} | Erreur absolue | ${TITLE}"
  set grid
  red = "#FF0000"; green = "#00FF00"; blue = "#0000FF"; skyblue = "#87CEEB";
  set style data histogram
  set style histogram cluster gap 1
  set style fill solid
  set boxwidth 0.9
  set grid ytics
  set xtics ('Etat Final abs' 0, 'Etat Choc rel' 1)
  set ylabel '${liste_phase_name[j]}'
EOF

  echo -n plot >> $name_fich
  for (( i=0 ; i<nsimu-1 ; i++ )) do
    echo -n " 'err_ligne_abs_courbe_conv${i}' using ${!var} title '${liste_name[i]}'," >> $name_fich
  done
  echo " 'err_ligne_abs_courbe_conv${i}' using ${!var} title '${liste_name[i]}'" >> $name_fich

done
###############################################################################################


#################################################################################################
## --------------   boucle BARPLOT sur les erreurs variables RELATIVES  -----------------------------
for (( j=0 ; j<nplot ; j++ ))  # boucle sur les variables
do
  var=${liste_plot[j]}
  name_fich="fichiers/GNUPLOT/err_${var}_rel.plot"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
  set term x11 font "Arial,20"
  set title "${liste_plot[j]} | Erreur relative | ${TITLE}"
  set grid
  red = "#FF0000"; green = "#00FF00"; blue = "#0000FF"; skyblue = "#87CEEB";
  set style data histogram
  set style histogram cluster gap 1
  set style fill solid
  set boxwidth 0.9
  set grid ytics
  set xtics ('Etat Final rel' 0, 'Etat choc rel' 1)
  set ylabel '${liste_phase_name[j]}'
EOF

  echo -n plot >> $name_fich
  for (( i=0 ; i<nsimu-1 ; i++ )) do
    echo -n " 'err_ligne_rel_courbe_conv${i}' using ${!var} title '${liste_name[i]}'," >> $name_fich
  done
  echo " 'err_ligne_rel_courbe_conv${i}' using ${!var} title '${liste_name[i]}'" >> $name_fich

done



#************************#*************************#************************#*******************
#************************ Creation des scripts PLOT **********************************************

point_size="0.4"
liste_type_erreur=(
 "EF_abs"
 "EF_rel"
 "EP_abs"
 "EP_rel"
)
liste_type_erreur_name=(
 "Erreur absolue sur Etat Final"
 "Erreur relative sur Etat Final"
 "Erreur absolue sur Etat Palier"
 "Erreur relative sur Etat Palier"
)
 
## --------------   boucle PLOT sur Etat Final en erreur Absolue   -----------------------------
for (( k=0 ; k<4 ; k++ ))  do # boucle sur les types d'erreurs
  for (( j=0 ; j<nplot ; j++ ))  # boucle sur les variables
  do
    var=${liste_plot[j]}
    #printf "var= %s  et num= %d\n" ${var} ${!var} 
    name_fich="fichiers/GNUPLOT/err_${liste_type_erreur[k]}_${liste_plot[j]}.plot"
    #printf "name_fich= %s\n" ${name_fich}
   
    #echo $j $var $name_fich

cat <<EOF > $name_fich
  set term x11 enhanced font "Arial,20" 
  set title "${TITLE} | ${liste_type_erreur_name[k]}"
  set xlabel "nombre de mailles"
  set ylabel '${liste_plot_name[j]}'
  set grid
EOF
    echo -n plot >> $name_fich
    for (( i=0 ; i<nsch ; i++ )) # boucle sur les simulations
      do
        echo -n " 'err_ligne_${liste_type_erreur[k]}$i' using 1:${!var} w lp pointsize ${point_size}  title '${liste_schemes_name[i]}'," >> $name_fich
      done
    #echo " 'err_ligne_EF_abs$i' using 1:${!var} w lp pointsize ${point_size}  title '${liste_sch[i]}'" >> $name_fich
  done
done
###############################################################################################


#************************#************************#************************#*******************
#************************* Generation des graphes **********************************************

cd fichiers/GNUPLOT
#gnuplot "rho.plot" -persist  # la fenetre n'est plus interactive !!?


gnuplot "erreur.plot" -persist;


# plot variable
#for (( j=0 ; j<nplot ; j++ )) do
#  gnuplot "err_${liste_plot[j]}_abs.plot" -persist
#done

#for (( j=0 ; j<nplot ; j++ )) do
#  gnuplot "err_${liste_plot[j]}_rel.plot" -persist
#done

j=${!rho}
gnuplot "err_${liste_plot[j]}_rel.plot" -persist
j=${!P}
gnuplot "err_${liste_plot[j]}_rel.plot" -persist

# Courbe de convergence
k=1 # Etat Final relatif
for (( j=0 ; j<nplot ; j++ )) do
  gnuplot "err_${liste_type_erreur[k]}_${liste_plot[j]}.plot" -persist
done

#************************#************************#************************#*******************
#************************ fin Creation des diagrammes erreurs  **********************************************


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





#  Supression des fichiers textes
#for (( i=0 ; i<nsch ; i++ ))
#do
#  rm err_ligne_EF_abs$i err_ligne_EF_rel$i err_ligne_EP_abs$i err_ligne_EP_rel$i
#done

# Ils sont supprimer avant de lancer une simyualtion pour pouvoir ré-écrire dessus



