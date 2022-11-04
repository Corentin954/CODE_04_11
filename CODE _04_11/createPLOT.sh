#!/bin/bash

declare -i i nsimu tst
declare -a liste_name

# On lit le fichier avec les parametres de la simulation
i="-3"
while read line
do
  if [[ $i -eq -3 ]]; then	TITLE=${line}; 	fi
  if [[ $i -eq -2 ]]; then	tst=${line}; 	fi
  if [[ $i -eq -1 ]]; then 	nsimu=${line}; 	fi
  if [[ -1 -lt $i ]]; then	liste_name[i]=${line}; 	fi
  (( i=${i}+1 ));
done < param_sim.txt

#printf "TITLE=${TITLE}\n"
#printf "tst=${tst}\n"
#printf "nsimu=${nsimu}\n"
#declare -p liste_name


rho="2"; tau="3"; e="4"; e_interne="5"; P="6"; T="7"; C="8"; G="9"; S="10";
declare -ra liste_plot=(rho tau e e_interne P T C G S)
declare -ra liste_plot_name=("Densité (kg/m)" "Volume spécifique (m/kg)" "Energie totale (J/kg)" "Energie interne (J/kg)" "Pression (Pa)" "Temperature (K)" "Celerite du son (m/s)" "Dérivée fondamentale" "Entropie")
declare -i nplot=${#liste_plot[@]}

beta="2"; gamma="3"; liquide="4"; 
declare -ra liste_phase=(beta gamma liquide)
declare -ra liste_phase_name=("phase BETA" "phase GAMMA" "phase LIQUIDE")


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

point_size="0.4"

## --------------   boucle PLOT sur les variables centrées aux mailles  -----------------------------
for (( j=0 ; j<nplot ; j++ ))  # boucle sur les variables
do
  var=${liste_plot[j]}
  name_fich="fichiers/GNUPLOT/$var.plot"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
  set term x11 enhanced font "Arial,20" 
  set title "${TITLE}"
  set xlabel "X (m)"
  set ylabel '${liste_plot_name[j]}'
  set grid
EOF
  echo -n plot >> $name_fich
  if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || ((100 -lt ${tst} && ${tst} -lt 112)) ]]; then
    echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact'," >> $name_fich
  fi
  for (( i=0 ; i<nsimu-1 ; i++ )) # boucle sur les simulations
  do
    echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size}  title '${liste_name[i]}'," >> $name_fich
  done
  echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size}  title '${liste_name[i]}'" >> $name_fich
done
###############################################################################################



## --------------   PLOT de la vitesse -----------------------------------------------------------
name_fich="fichiers/GNUPLOT/u.plot"
cat <<EOF > $name_fich
	set term x11 font "Arial,20"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel "u (m/s)"
	set grid
EOF
echo -n plot >> $name_fich
if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || ((100 -lt ${tst} && ${tst} -lt 112)) ]]; then
  echo -n " 'sol_reference_u.txt' w l title 'exact'," >> $name_fich
fi
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data_u${i}' w lp title '${liste_name[i]}'," >> $name_fich
done
echo " 'data_u${i}' w lp title '${liste_name[i]}'" >> $name_fich
###############################################################################################""


## --------------   PLOT de les fractions massiques -----------------------------------------------------------
for (( j=0 ; j<3 ; j++ )) do   
	var=${liste_phase[j]}
  name_fich="fichiers/GNUPLOT/$var.plot"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
	set term x11 enhanced font "Arial,20"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel '${liste_phase_name[j]}'
	set grid
EOF

  echo -n plot >> $name_fich
  if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || ((100 -lt ${tst} && ${tst} -lt 112)) ]]; then
    echo -n " 'sol_reference_lambda.txt' using 1:${!var} w l title 'exact'," >> $name_fich
  fi
  for (( i=0 ; i<nsimu-1 ; i++ )) do # boucle sur les phases
    echo -n " 'data_lambda${i}' using 1:${!var} w lp title '${liste_name[i]}'," >> $name_fich
  done
  echo " 'data_lambda${i}' using 1:${!var} w lp title '${liste_name[i]}'" >> $name_fich
done
###############################################################################################""


##--- PLOT : zoom sur les niveaux d'energie interne-----------
var="ezoom"
name_fich="fichiers/GNUPLOT/${var}.plot"
ezoom="5"

declare xrange yrange
if [[ ${tst} -eq 0 ]]; then   xrange="0.67:0.885"; yrange="2.8:2.9"; fi   # Sod
if [[ ${tst} -eq 102 ]]; then xrange="-6e-4:0";   yrange="5e4:5.8e4"; fi  # Etain tst=102

cat <<EOF > $name_fich
	set terminal x11 font "Arial,20"
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m) (zoom)"
	set ylabel "e interne (J/kg)"
	set grid
	set xrange[${xrange}]
	set yrange[${yrange}]
EOF
echo -n plot >> $name_fich
if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || ((100 -lt ${tst} && ${tst} -lt 112)) ]]; then
  echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact'," >> $name_fich
fi
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################""


##--- PDF : zoom sur les niveaux de pressions -----------
# Pression
name_fich="fichiers/GNUPLOT/Pzoom.plot"
var="Pzoom"
Pzoom="6"
declare xrange yrange
if [[ ${tst} -eq 0 ]]; then   xrange="0.45:0.88"; yrange="0.3025:0.3035"; fi  # Sod
if [[ ${tst} -eq 102 ]]; then xrange="-6e-4:0";   yrange="7.5e9:8e9"; fi  # Etain tst=102

cat <<EOF > $name_fich
	set terminal x11 font "Arial,20"
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m) (zoom)"
	set ylabel "P (Pa)"
	set grid
	set xrange[${xrange}]
	set yrange[${yrange}]
EOF
#set xrange[-6e-4:0]
#set yrange[7.5e9:8e9]
echo -n plot >> $name_fich
if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || ((100 -lt ${tst} && ${tst} -lt 112)) ]]; then
  echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact'," >> $name_fich
fi
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################""



##--- PLOT : nombre d'itération du Newton-Raphson sur la dernière itération-----------

var="nb_iteration"
name_fich="fichiers/GNUPLOT/$var.plot"
cat <<EOF > $name_fich
set terminal x11  font "Arial,20"
set output "$var.pdf"
set title "${TITLE}"
set xlabel "maille"
set ylabel "nombre d'itération du Newton-Raphson sur la dernière itération"
set grid
EOF
echo -n plot >> $name_fich
for (( i=0 ; i<nsimu ; i++ )) do
  echo -n " 'nb_iter${i}' using 1 w p pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done


###############################################################################################""

##--- PLOT : valeur du critère (epsilon) pour l'arret du Newton-Raphson-----------

var="critère_Newton_Raphson"
name_fich="fichiers/GNUPLOT/$var.plot"
cat <<EOF > $name_fich
set terminal x11  font "Arial,20"
set output "$var.pdf"
set title "${TITLE}"
set xlabel "maille"
set ylabel "critère du Newton-Raphson sur la dernière itération"
set grid
EOF
echo -n plot >> $name_fich
for (( i=0 ; i<nsimu ; i++ )) do
  echo -n " 'nb_iter${i}' using 2 w p pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done


###############################################################################################""



#************************#************************#************************#*******************
#************************ Generation des graphes **********************************************

cd fichiers/GNUPLOT
#gnuplot "rho.plot" -persist  # la fenetre n'est plus interactive !!?

# plot variable
gnuplot "u.plot" -persist
for (( j=0 ; j<nplot ; j++ )) do
	gnuplot "${liste_plot[j]}.plot" -persist
done

# plot phase
if [[ 100 -lt ${tst} ]]; then
	for (( j=0 ; j<3 ; j++ )) do
		gnuplot "${liste_phase[j]}.plot" -persist;
	done
	gnuplot "GAMMAzoom.plot" -persist;
fi

# PLOT ZOOM
gnuplot "Pzoom.plot" -persist


#if [[ ${tst} -eq 0 ]]; then fi
gnuplot "ezoom.plot" -persist;


gnuplot "nb_iteration.plot" -persist;
gnuplot "critère_Newton_Raphson.plot" -persist;



