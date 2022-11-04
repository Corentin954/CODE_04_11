*****************************************************************************	
*****               Description des fichiers du projet              ******

*****************************************************************************	




!!! Necessite la bibliothèque     GSL - GNU Scientific Library  !!!
pour l'intégration d'ordre élevé des solutions initiales.
voir : https://www.gnu.org/software/gsl
Dans la fonction func.c pour le calcul des valeurs moyennes initiales : #include <gsl/gsl_integration.h> 


$$ 
Tout les fichiers textes ou binaires se trouvent dans le dossier fichier sauf :
    - les solutions de main_multi.c   sol.txt : P,RHO,U etc.. sur les mailles  
                                      NRJ.txt : diagnocstic de conservation des schémas
                                      frac_mass.txt : fractions massiques sur les mailles
    - la solution de cal_PU.c         PU.txt : valeur de pression et vitesse des tubes à choc.
$$




-------------------------------------------------------------------------------------------
-- Fichier code C --
	
Code langrangien

> main_multi.c : fichier depuis lequel on lance les simulations avec les arguments en ligne de commande.
                 --> compil_main_multi.sh : pour créer l'éxecutable

> func***_multi.c : source des schémas qui sont lancée par le main (c'est d'ici qu'on lance l'appel à l'EOS)

> func.c : contient les routines nécessaires aux schémas : écritures des résultats dans les fichiers textes, conditions aux limites, pseudo-viscosité, limiteurs, reconstruction spatiales, table de Butcher, fonction d'initialisation.

EOS 

> EOS.c : contient l'équation d'état gaz parfait et Bizarium (les source-schémas appellent directement les fonctions contenues dans le fichier)

> EOS_etain : contient toutes les fcts de la fromulation Mie-Gruneisen des 3 phase (P,TS,G et leurs dérivées partielles et les constantes)

> cal_front.c : (6000 lignes de code) contient toutes les fonctions pour la crétaion de l'équation d'état multiphase de l'étain :
                 - les newtons pour calculer le points triple et les lignes de changement de phases
                 - les fonctions qui écrivent les points des frontières de zones dans des fichiers
                 - les fonctions de détections de zones
                 - les fonctions qui calcul la P,T dans les zones mixtes (avec Newton)
                 - fonctions pour tracer P,T en fct de V,E
                 - fonctions pour créer la tabulation et la lire
                 - puis la fonction "fPTC_cin_phase" qui est utilisé dans les script-schémas pour le calcul de P et T.

> hugoniot.c : regroupe les fonctions pour tracer les isentropes, les points de l'onde de compression isentropiques. 
               Et aussi les points caractéristiques de la courbe d'Hugoniot et les fonctions pour écrire les points de la courbe dans un fichier.
               
> main_calfront.c : qui permet de lancer les fonctions de contenues dans cal_front.c et hugoniot.c. 
                    Les lignes pour lancer le calcul des hugoniots, le graphe de la pression ou la création de la tabulation sont mis en commentaire dans le main 
                    --> compil_calfront.sh : pour créer l'éxecutable

Autres fichiers C

> cal_PU.c : permet de calculer la pression et la vitesse de l'état après la discontinuité de contact dans les tubes à choc (Sod et LeBlanc) par un Newton et les écrits dans un fichier texte.
             Utilisé pour tracé les solutions exactes
             --> compil_PU.sh : pour créer l'éxecutable

> allocfree.c : allocation et libération des matrices.



-------------------------------------------------------------------------------------------------------------------------------
-- Header --

> Matrice.h : permet de résoudre un système linéaire (utilisé dans les Newton de calfront.c et hugoniot.c)

> head_multi.h : déclaration des fonctions pour les schémas + EOS

> EOS_etain.h : déclaration des fonctions de la formulation Mie-Gruneisen




-------------------------------------------------------------------------------------------------------------------------------
-- Script Octave --

> trac_main_multi.m :  permet d'obtenir les figures de simulations depuis Octave.
                       Le script compile avec le fichier .sh puis lance dans une boucle l'éxecutable avec les différents arguments à chaque fois. 
                       Les arguments qui peuvent être modifier à chaque itérations de la boucles sont : 
                         le schéma, ordre spatiale, le DeBar fix, la cinétique et la valeur de tau.
                       Les autres arguments sont les mêmes pour chaque itérations de la boucle (CFL, nombre de mailles, pseudo, ect..)

> trac_main_multi_nx.m :  permet d'obtenir les figures de simulations depuis Octave. Cette fois ci en faisant varier le nombre de mailles.

> prise_ordre.m : permet d'obtenir les courbes (log nx, log Err L1) sur le cas-test de l'onde acoustique.

> plotPOLYGONES.m : donne le diagramme de phase en variables V,E et P,T.

> plot_hugoniot.m : graphes des hugoniots 

> trac_PT_cinphase.m : graphe de P,T,C et des fractions massiques à l'équilibre thermo en fonction de VE.

> lec_save_nx_diff.m :  permet de tracer des résultats d'étude de convergence utlisé dans le rapport grace à des sauvegarde de l'environnment octave après les simulations.
                        Les sauvegardes (assez lourdes) sont dans le fichier SAVE.
                        
> lec_save.m :  Idem que le script précedent sauf que simulations sauvgarder ont le même nombre de mailles.
                        