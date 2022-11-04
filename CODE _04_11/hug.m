clear;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    (rho0,E0,Sr) Heuzé du 1er jeu de coeff
% [Point 0 , Point 1 , Point 2 , Point 3 , Point 4 , Point 5]
i=1;
 Ppt(i,:)=[0 7.67399 8.47806 11.3538 44.8666 56.9035];  
 Tpt(i,:)=[300 404.148 354.881 388.436 1965.67 2231.31];  
 RHOpt(i,:)=[7287 8076.01 8288.37 8518.28 10247.9 10609];  
 Vpt(i,:)=[137.231 123.824 120.651 117.395 97.5813 94.2595];  
 Ept(i,:)=[-3.40245e-10 51443 77064.9 112608 889467 1.2226e+06];  
 Upt(i,:)=[-99 320.758 371.266 474.57 1333.77 1563.72]; 
 Dpt(i,:)=[-99 3283.18 3283.18 3283.18 4616.31 4993.81];  
 Cpt(i,:)=[2740.55 1960.14 1926.14 3581.51 4624.85 4829.42];   
 Spt(i,:)=[2.69267e-05 16.0072 16.0167 25.0687 305.216 383.62];   
 G_DIFFpt(i,:)=[2.24 3.07344 2.58963 2.52588 2.05781 1.95532];  

% lecture fichiers
%  fprintf(file, "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n",V,E,P,T,U,D,S,C,Der_Fond,1-tabX[i],tabX[i],0.0 );
file_hug=load("fichiers/hug0.txt");
[~,n0]=size(file_hug);
Vhug(:,i)=file_hug(:,1);  Ehug(:,i)=file_hug(:,2);   Phug(:,i)=file_hug(:,3);    Thug(:,i)=file_hug(:,4);   Uhug(:,i)=file_hug(:,5);   Dhug(:,i)=file_hug(:,6);
Shug(:,i)=file_hug(:,7);  Chug(:,i)=file_hug(:,8);   Der_Fondhug(:,i)=file_hug(:,9);  XBETAhug(:,i)=file_hug(:,10);  XGAMMAhug(:,i)=file_hug(:,11);  XLIQhug(:,i)=file_hug(:,12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    (rho0,E0,Sr) Heuzé du 2nd jeu de coeff convergés
% [Point 0 , Point 1 , Point 2 , Point 3 , Point 4 , Point 5]
i=2;
 Ppt(i,:)=[0 7.25928 8.71925 11.6371 48.246 56.1677];  
 Tpt(i,:)=[300 400.239 331.954 369.935 2322.83 2503.13];  
 RHOpt(i,:)=[7287 8080.87 8391.79 8649.12 10585.6 10833.8];  
 Vpt(i,:)=[137.231 123.749 119.164 115.619 94.4678 92.3033];  
 Ept(i,:)=[-3.59879e-10 48933.7 85564 125751 1.03157e+06 1.26173e+06];  
 Upt(i,:)=[-99 312.838 394.654 501.499 1436.36 1588.54]; 
 Dpt(i,:)=[-99 3184.39 3184.39 3184.39 4609.44 4852.2];  
 Cpt(i,:)=[2712.04 2217.5 2112.51 3417.03 4602.73 4651.91];   
 Spt(i,:)=[2.69267e-05 13.7072 13.6748 25.7704 348.025 396.542];   
G_DIFFpt(i,:)=[1.863 2.16818 2.193 2.64383 2.10714 1.9106];  

file_hug=load("fichiers/hug1.txt");
[~,n0]=size(file_hug);
Vhug(:,i)=file_hug(:,1);  Ehug(:,i)=file_hug(:,2);   Phug(:,i)=file_hug(:,3);    Thug(:,i)=file_hug(:,4);   Uhug(:,i)=file_hug(:,5);   Dhug(:,i)=file_hug(:,6);
Shug(:,i)=file_hug(:,7);  Chug(:,i)=file_hug(:,8);   Der_Fondhug(:,i)=file_hug(:,9);  XBETAhug(:,i)=file_hug(:,10);  XGAMMAhug(:,i)=file_hug(:,11);  XLIQhug(:,i)=file_hug(:,12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    (rho0,E0,Sr) Heuzé du 2nd jeu de coeff convergés
% [Point 0 , Point 1 , Point 2 , Point 3 , Point 4 , Point 5]
i=3;
 Ppt(i,:)=[0 7.21366 8.74846 11.6019 49.9332 57.6626];  
 Tpt(i,:)=[300 399.446 329.73 366.856 2446.58 2626.69];  
 RHOpt(i,:)=[7287 8076.77 8395.26 8646.86 10650 10887.2];  
 Vpt(i,:)=[137.231 123.812 119.115 115.649 93.8969 91.8507];  
 Ept(i,:)=[-3.59879e-10 48399.3 85887.3 125195 1.0819e+06 1.30836e+06];  
 Upt(i,:)=[-99 311.125 396.031 500.39 1470.98 1617.63]; 
 Dpt(i,:)=[-99 3181.8 3181.8 3181.8 4658.36 4891.77];  
 Cpt(i,:)=[2712.04 2245.97 2140.33 3415.31 4637.83 4677.46];   
 Spt(i,:)=[2.69267e-05 13.5091 13.4783 25.4693 358.573 404.415];   
G_DIFFpt(i,:)=[1.863 2.20354 2.23879 2.64479 2.09254 1.89926];  

file_hug=load("fichiers/hug2.txt");
[~,n0]=size(file_hug);
Vhug(:,i)=file_hug(:,1);  Ehug(:,i)=file_hug(:,2);   Phug(:,i)=file_hug(:,3);    Thug(:,i)=file_hug(:,4);   Uhug(:,i)=file_hug(:,5);   Dhug(:,i)=file_hug(:,6);
Shug(:,i)=file_hug(:,7);  Chug(:,i)=file_hug(:,8);   Der_Fondhug(:,i)=file_hug(:,9);  XBETAhug(:,i)=file_hug(:,10);  XGAMMAhug(:,i)=file_hug(:,11);  XLIQhug(:,i)=file_hug(:,12);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Diagramme (V,E)  et  (P,T)   
%    correspondant aux coeff (rho0,E0,Sr) du second jeu de ceoff de Heuzé
SA=load("fichiers/pSA.txt");
SB=load("fichiers/pSB.txt");
SC=load("fichiers/pSC.txt");
SAB=load("fichiers/pSAB.txt");
SAC=load("fichiers/pSAC.txt");
SBC=load("fichiers/pSBC.txt");

SApt=load("fichiers/pSApt.txt");
SBpt=load("fichiers/pSBpt.txt");
SCpt=load("fichiers/pSCpt.txt");
SABpt=load("fichiers/pSABpt.txt");
SACpt=load("fichiers/pSACpt.txt");
SBCpt=load("fichiers/pSBCpt.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Valeur de reference

% Al’tshuler et al. (1962)
D1962=[4.20  6.36   9.02   12.98]*1e3; 
U1962=[1.08  2.46  4.65  7.96]*1e3; 
%Al’tshuler et al. (1981)
D1981=[7.59 8.23  8.28 12.57]*1e3; 
U1981=[3.34   3.86  3.93  7.71]*1e3; 

% 
rho0_tab1=[ 7.291  7.294  7.290  7.289  7.290  7.280  7.291  7.280  7.290  7.292  7.280  7.290  7.280  7.289  7.289  7.280  7.290  7.289  7.289];
Us_tab1=[ 2.755  3.273 3.299 3.263 3.361 3.482 3.476 3.514 3.524 3.658 3.956 4.102 4.382 4.309 4.292 4.500 4.587 4.755 4.766]*1e3;
Up_tab1=[ 0.000 0.304 0.418 0.451 0.577 0.642 0.645 0.677 0.693 0.767 0.896  1.073 1.161 1.171 1.201 1.260 1.361 1.456 1.465]*1e3;
P_tab1=[ 0.000 7.257 10.053 10.727 14.137 16.274 16.347 17.319 17.803 20.459 25.805 32.087 37.037 36.779 37.573 41.278 45.511 50.464 50.893];
V_tab1=[ 0.1372 0.1244 0.1198 0.1182 0.1136 0.1120 0.1117 0.1109 0.1102 0.1084 0.1063 0.1013 0.1010 0.0999 0.0988 0.0989 0.0965 0.0952 0.0950 ];
rho_tab1=[ 7.291 8.041 8.348 8.458 8.801 8.926 8.952 9.017 9.075 9.227 9.412 9.872 9.904 10.009 10.121 10.111 10.366 10.506 10.524];

rho0_tab2=[7.286 7.290 7.289 7.290 7.280 7.280 7.280 7.289 7.289 7.289 7.289 7.291 7.290 7.280 7.280 7.280 7.289 7.290 7.289 7.289 7.289 7.290 ];
Us_tab2=[ 4.793 5.239 5.083 5.331 5.527 5.504 5.481 5.458 5.681 5.717 5.752 6.193 6.500 6.721 6.626 6.664 6.765 6.849 6.858 6.935 7.068 7.236 ]*1e3;
Up_tab2=[ 1.493 1.665 1.669 1.734 1.926 1.929 1.933 1.940 2.002 2.044 2.082 2.461 2.564 2.751 2.764 2.777 2.780 2.859 2.893 2.942 3.038 3.118 ]*1e3;
P_tab2=[ 52.138 63.590 61.836 67.388 77.496 77.293 77.130 77.180 82.609 85.176 87.291 111.122 121.495 134.603 133.328 134.723 137.082 142.748 144.615 148.716 156.514 164.476];
V_tab2=[ 0.0945 0.0936 0.0921 0.0926 0.0895 0.0892 0.0889 0.0884 0.0887 0.0881 0.0875 0.0827 0.0831 0.0811 0.0801 0.0801 0.0808 0.0799 0.0793 0.0790 0.0782 0.0781];
rho_tab2=[ 10.852 10.686 10.852 10.804 11.174 11.208 11.246 11.309 11.277 11.345 11.424 12.099 12.039 12.325 12.490 12.481 12.374 12.514 12.607 12.659 12.784 12.810 ];


TAB1=[transpose(rho0_tab1) , transpose(Up_tab1) , transpose(Us_tab1) , transpose(P_tab1) , transpose(rho_tab1)];
TAB2=[transpose(rho0_tab2) , transpose(Up_tab2) , transpose(Us_tab2) , transpose(P_tab2) , transpose(rho_tab2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Lnombre= ['0','1','2','3','4','5'];
#Lcouleur= ['m','r','blue'];

Lcouleur=[ 0.5,0.5,0 ; 
           1,0.5,0.5 ; 
           0.5,0,0.5 ;
           0,0.5,0.5
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             plot  (V,E)
figure
hold on

% hugoniot
for i = 1:3
  plot(Vhug(:,i)*1e6,Ehug(:,i),"linewidth",2,'color',Lcouleur(i,:));
endfor

% Points caractéristiques 
for i = 1:3
  plot(Vpt(i,:),Ept(i,:),'linestyle','none','marker','+', "linewidth",7,'color',Lcouleur(i,:));
  for j=1:6
    text(Vpt(i,j)-0.5,Ept(i,j)-7e4,Lnombre(j),"fontsize", 40,"color",Lcouleur(i,:));
  endfor
endfor



##% ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');

%xlim([80 160]);
%ylim([-5e3 2e6]);
xlabel('\tau (cm^3/kg)');
ylabel('\epsilon (J/kg)');
grid()
set(gca, "fontsize", 40)
legend('1','2','3');
title('Hugoniot (\tau,\epsilon)');
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             plot (P,T)
figure
hold on

% Hugoniot
for i = 1:3
  plot(Phug(:,i)*1e-9,Thug(:,i),"linewidth",1.5,"color",Lcouleur(i,:));
endfor

% Points caractéristiques
for i = 1:3
  plot(Ppt(i,:),Tpt(i,:),'linestyle','none','marker','+', "linewidth",5,"color",Lcouleur(i,:));
  for j=1:6
    text(Ppt(i,j)-0.5,Tpt(i,j)-90,Lnombre(j),"fontsize", 40,"color",Lcouleur(i,:));
  endfor
endfor

% ploynomes
##S=SApt;
##plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');
##S=SBpt;
##plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');
##S=SCpt;
##plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');

%xlim([-2 100]);
%ylim([0 3000]);
xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40)
legend('1','2','3');
title("Hugoniot (P,T)");
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             plot (V,P)
figure
hold on

% Hugoniot
for i = 1:3
  plot(Vhug(:,i)*1e6,Phug(:,i)*1e-9,"linewidth",1.5,"color",Lcouleur(i,:));
endfor


% Points caractéristiques
for i = 1:3
  plot(Vpt(i,:),Ppt(i,:),'linestyle','none','marker','+', "linewidth",5,"color",Lcouleur(i,:));
  for j=1:6
    text(Vpt(i,j)-0.5,Ppt(i,j)-2.5,Lnombre(j),"fontsize", 40,"color",Lcouleur(i,:));
  endfor
endfor
% Polygones
##P=SApt(:,2);
##V=SA(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');
##P=SBpt(:,2);
##V=SB(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');
##P=SCpt(:,2);
##V=SC(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');

ylabel("P (GPa)");
xlabel('\tau (cm^3/kg)');
grid()
set(gca, "fontsize", 40)
legend('1','2','3');
title('Hugoniot (\tau,P)');
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             plot (u,D)
figure
hold on

% Hugoniot
for i = 1:3
  plot(Uhug(:,i),Dhug(:,i),"linewidth",1.5,"color",Lcouleur(i,:));
endfor

% Points caractéristiques
for i = 1:3
  plot(Upt(i,:),Dpt(i,:),'linestyle','none','marker','+', "linewidth",5,"color",Lcouleur(i,:));
  for j=1:6
    text(Upt(i,j),Dpt(i,j),Lnombre(j),"fontsize", 40,"color",Lcouleur(i,:));
  endfor
endfor
ylabel("D (m/s)");
xlabel('u (m/s)');
grid()
set(gca, "fontsize", 40)
legend('1','2','3','1962','1981','1980 tab1','1980 tab2');
title('Hugoniot (u,D)');
hold off




###############################################################################
#      Comparaison aux valeurs expérimentales
###############################################################################



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             plot (u,D)
figure
hold on

% Hugoniot
for i = 1:3
  plot(Uhug(:,i),Dhug(:,i),"linewidth",1.5,"color",Lcouleur(i,:));
endfor

plot(U1962,D1962,'linestyle','none',"marker",'+','color','r',"linewidth",2);
plot(U1981,D1981,'linestyle','none',"marker",'o','color','blue',"linewidth",2);

plot(Up_tab1,Us_tab1,'linestyle','none',"marker",'d','color','m',"linewidth",2);
plot(Up_tab2,Us_tab2,'linestyle','none',"marker",'s','color','black',"linewidth",2);

##% Points caractéristiques
##for i = 1:3
##  plot(Upt(i,:),Dpt(i,:),'linestyle','none','marker','+', "linewidth",5,"color",Lcouleur(i,:));
##endfor
##
##for j=1:6
##   text(Upt(i,j),Dpt(i,j)-2e2,Lnombre(j),"fontsize", 40,"color",Lcouleur(i,:));
##endfor

% Points caractéristiques
for i = 1:3
 # plot(Upt(i,:),Dpt(i,:),'linestyle','none','marker','+', "linewidth",5,"color",Lcouleur(i,:));
endfor
ylabel("D (m/s)");
xlabel('u (m/s)');
grid()
set(gca, "fontsize", 40)
legend('1','2','3','1962','1981','1980 tab1','1980 tab2');
legend('location','northwest');
title('Hugoniot (u,D)');
hold off
