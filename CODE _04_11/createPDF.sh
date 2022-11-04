#!/bin/bash

declare -i i
declare -a liste_name
declare -i nsimu

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


#************************#************************#************************#*******************
#************************ Creation des scripts PDF  **********************************************

point_size="0.4"

## --------------   boucle PDF sur les variables centrées aux mailles  -----------------------------
for (( j=0 ; j<nplot ; j++ )) do # boucle sur les variables
	var=${liste_plot[j]}
  name_fich="fichiers/GNUPLOT/$var.tpdf"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
	set terminal pdfcairo color
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel '${liste_plot_name[j]}'
	set grid
EOF

	echo -n plot >> $name_fich
	if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || 100 -lt ${tst} ]]; then
	  echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact'," >> $name_fich
	fi
	for (( i=0 ; i<nsimu-1 ; i++ )) do # boucle sur les simulations
		echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
	done
	echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich

done
###############################################################################################""


## --------------   PDF de la vitesse -----------------------------------------------------------
name_fich="fichiers/GNUPLOT/u.tpdf"
cat <<EOF > $name_fich
	set terminal pdfcairo color
	set output "u.pdf"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel "u (m/s)"
	set grid
EOF
echo -n plot >> $name_fich
if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || 100 -lt ${tst} ]]; then
  echo -n "'sol_reference_u.txt' w l title 'exact' ," >> $name_fich
fi
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data_u${i}' w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data_u${i}' w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################


## --------------   PDF sur les variables centrées aux mailles  -----------------------------
for (( j=0 ; j<nplot ; j++ )) do # boucle sur les variables
	var=${liste_plot[j]}
  name_fich="fichiers/GNUPLOT/$var.tpdf"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
	set terminal pdfcairo color font "Arial,16" linewidth 0.1 size 8,4
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel '${liste_plot_name[j]}'
	set grid
EOF

	echo -n plot >> $name_fich
	if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || 100 -lt ${tst} ]]; then
	  echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact'," >> $name_fich
	fi
	for (( i=0 ; i<nsimu-1 ; i++ )) do # boucle sur les simulations
		echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
	done
	echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich

done
###############################################################################################""

## --------------   PDF sur les fractions massiques  -----------------------------
for (( j=0 ; j<3 ; j++ )) do # boucle sur les variables
	var=${liste_phase[j]}
  name_fich="fichiers/GNUPLOT/$var.tpdf"

  #echo $j $var $name_fich

cat <<EOF > $name_fich
	set terminal pdfcairo color font "Arial,16" linewidth 0.1 size 8,4
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m)"
	set ylabel '${liste_plot_name[j]}'
	set grid
EOF

	echo -n plot >> $name_fich
	if [[ ${tst} -eq 0 || ${tst} -eq 1 || ${tst} -eq 3 || 100 -lt ${tst} ]]; then
	  echo -n " 'sol_reference_lambda.txt' using 1:${!var} w l title 'exact'," >> $name_fich
	fi
	for (( i=0 ; i<nsimu-1 ; i++ )) do # boucle sur les simulations
		echo -n " 'data_lambda${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
	done
	echo " 'data_lambda${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich

done
###############################################################################################""


##--- PDF : zoom sur les niveaux de pressions -----------
# Pression
name_fich="fichiers/GNUPLOT/Pzoom.tpdf"
var="Pzoom"
Pzoom="6"
declare xrange yrange
if [[ ${tst} -eq 0 ]]; then   xrange="0.45:0.88"; yrange="0.3025:0.3035"; fi  # Sod
if [[ ${tst} -eq 102 ]]; then xrange="-6e-4:0";   yrange="7.5e9:8e9"; fi  # Etain tst=102
#if [[ ${tst} -eq 102 ]]; then xrange="-6e-4:-3.5e-4";   yrange="7.6735e9:7.678e9"; fi  # Etain tst=102

cat <<EOF > $name_fich
	set terminal pdfcairo color font "Arial,16" linewidth 0.01 size 8,4
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
echo -n " 'sol_reference.txt' using 1:${!var} w l title 'exact' ," >> $name_fich
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################""


##--- PDF : zoom sur les niveaux des fractions massiques pour le cas-test 102-----------
# Fractions massiques
var="GAMMAzoom"
GAMMAzoom="3"
name_fich="fichiers/GNUPLOT/$var.tpdf"

if [[ ${tst} -eq 102 ]]; then xrange="-6.5e-4:1e-4";   yrange="-1e-2:0.3"; fi  # Etain tst=102

cat <<EOF > $name_fich
	set terminal pdfcairo color font "Arial,16" linewidth 0.01 size 8,4
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m) (zoom)"
	set ylabel "fraction massique de la phase GAMMA"
	set grid
	set xrange[${xrange}]
	set yrange[${yrange}]
EOF
echo -n plot >> $name_fich
echo -n "'sol_reference_lambda.txt' using 1:${!var} w l title 'exact' ," >> $name_fich
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data_lambda${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data_lambda${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################""



##--- PDF : zoom sur les niveaux d'energie interne-----------
var="ezoom"
name_fich="fichiers/GNUPLOT/${var}.tpdf"
ezoom="5"

declare xrange yrange
if [[ ${tst} -eq 0 ]]; then   xrange="0.67:0.885"; yrange="2.8:2.9"; fi   # Sod
#if [[ ${tst} -eq 0 ]]; then   xrange="0.842:0.855"; yrange="1.999:3.2"; fi   # Sod
if [[ ${tst} -eq 102 ]]; then xrange="-6e-4:0";   yrange="5e4:5.8e4"; fi  # Etain tst=102

cat <<EOF > $name_fich
	set terminal pdfcairo color font "Arial,16" linewidth 0.01 size 8,4
	set output "$var.pdf"
	set title "${TITLE}"
	set xlabel "X (m) (zoom)"
	set ylabel "e interne (J/kg)"
	set grid
	set xrange[${xrange}]
	set yrange[${yrange}]
EOF
echo -n plot >> $name_fich
echo -n "'sol_reference.txt' using 1:${!var} w l title 'exact' ," >> $name_fich
for (( i=0 ; i<nsimu-1 ; i++ )) do
  echo -n " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'," >> $name_fich
done
echo " 'data${i}' using 1:${!var} w lp pointsize ${point_size} title '${liste_name[i]}'" >> $name_fich
###############################################################################################""




cd fichiers/GNUPLOT
#gnuplot "rho.plot" -persist  # la fenetre n'est plus interactive !!?


# PDF variable
gnuplot "u.tpdf"
for (( j=0 ; j<nplot ; j++ )) do
	gnuplot "${liste_plot[j]}.tpdf"
done

# PDF ZOOM
#if [[ ${tst} -eq 102 ]]; then	fi
gnuplot "Pzoom.tpdf"; 

gnuplot "ezoom.tpdf";



# PDF phase
if [[ 100 -lt ${tst} ]]; then
	for (( j=0 ; j<3 ; j++ )) do
		gnuplot "${liste_phase[j]}.tpdf";
	done;
	gnuplot "GAMMAzoom.tpdf";
fi