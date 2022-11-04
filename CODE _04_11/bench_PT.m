clear


P=load("fichiers/Ptest.txt");
T=load("fichiers/Ttest.txt");
X=load("fichiers/Xtest.txt");
V=load("fichiers/Vtest.txt");
E=load("fichiers/Etest.txt");
C=load("fichiers/Ctest.txt");
G=load("fichiers/Gtest.txt");

ZONE_test=load("fichiers/ZONEtest.txt");

[~,ne]=size(E);
[~,nv]=size(V);

P=reshape(P,ne,nv);
T=reshape(T,ne,nv);
Xbeta=reshape(X(:,1),ne,nv);
Xgamma=reshape(X(:,2),ne,nv);
Xliq=reshape(X(:,3),ne,nv);
ZONE_test=reshape(ZONE_test,nv,ne);
C=reshape(C,nv,ne);
G=reshape(G,nv,ne);

##Nb_ZONE=load("fichiers/Nb_ZONEbench.txt");
##V_dis=load("fichiers/Vbench_dis.txt");
##E_dis=load("fichiers/Ebench_dis.txt");
##[nv_dis,~]=size(V_dis);
##[ne_dis,~]=size(E_dis);
##Nb_ZONE=reshape(Nb_ZONE,nv_dis,ne_dis);


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


[na,~]=size(SA); % de taille NA+1 car SA[NA]=S[0] (code C)
[nb,~]=size(SB);
[nc,~]=size(SC);
[nab,~]=size(SAB); % de taille NBA+1 car SA[NA]=S[0] (code C)
[nac,~]=size(SAC);
[nbc,~]=size(SBC);


% Point triple
TRIPLE=load("fichiers/triple.txt");

Va=TRIPLE(1,1);  Ea=TRIPLE(1,2);
Vb=TRIPLE(2,1);  Eb=TRIPLE(2,2);
Vc=TRIPLE(3,1);  Ec=TRIPLE(3,2);
PT=TRIPLE(4,1);  TT=TRIPLE(4,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%          Diagramme de phase (P,T)
figure
hold on
%% ploynomes

S=SABpt;
plot(S(:,1),S(:,2), "linewidth",0.5,"color",'black');
S=SACpt;
plot(S(:,1),S(:,2), "linewidth",0.5,"color",'black');
S=SBCpt;
plot(S(:,1),S(:,2), "linewidth",0.5,"color",'black');

S=SApt;
plot(S(:,1),S(:,2), "linewidth",2,"color",'blue');
S=SBpt;
plot(S(:,1),S(:,2), "linewidth",2,"color",'magenta');
S=SCpt;
plot(S(:,1),S(:,2), "linewidth",2,"color",'red');


text(1.4e9,170,'\beta','color','blue',"fontsize", 30);
text(20e9,850,'\gamma','color','magenta',"fontsize", 30);
text(12.6e9,2000,'liquide','color','red',"fontsize", 30);

xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40)
title("Diagramme de phase (P,T)");
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%          Diagramme de phase (V,E)
figure
hold on
%% point
#Vp=0.000112814; Ep=77653.3;

% rectangle 
##Vmin=9.66856e-05; Vmax=9.71338e-05;  Emin=405485; Emax=428127;
##plot([Vmin Vmax]*1e6 ,[Emin Emin],"linewidth",1,"color",'red');
##plot([Vmin Vmax]*1e6 ,[Emax Emax],"linewidth",1,"color",'red');
##plot([Vmin Vmin]*1e6 ,[Emin Emax],"linewidth",1,"color",'red');
##plot([Vmax Vmax]*1e6 ,[Emin Emax],"linewidth",1,"color",'red');


% Coloriage des zones biphasiques
fill([Va Vb Vc]*1e6,[Ea Eb Ec],[0.8 0.8 0.8]);
S=SAB;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);
S=SAC;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);
S=SBC;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);



S=SAB;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'black');
S=SAC;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'black');
S=SBC;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'black');

S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'blue');
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'magenta');
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'red');


text(131,0,'\beta','color','blue',"fontsize", 45);
text(102,3.8e5,'\gamma','color','magenta',"fontsize", 45);
text(125,5e5,'liquide','color','red',"fontsize", 45);


xlabel(' \tau (cm^3/kg)');
ylabel(' \epsilon (J/kg)');
grid()
set(gca, "fontsize", 40)
legend('zone triphasique','zone biphasique');
title('Diagramme (\tau, \epsilon) ');
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%          P(V,E)  
figure
hold on
%surf(V*1e6,E,P*1e-9);
%%contour(V*1e6,E,P,50);
%%pcolor(V*1e6,E,P);
%contourf(V*1e6,E,P*1e-9,[0:2:70]);
contourf(V*1e6,E,P*1e-9,100);
% Un point 
%% ploynomes
S=SAB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SAC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SBC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("P(V,E)");
hold off


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##               T(V,E)
figure
hold on
%surf(V*1e6,E,T);
%%contour(V*1e6,E,T,50);
%%pcolor(V*1e6,E,T);
%contourf(V*1e6,E,T,[300:50:2500]);
contourf(V*1e6,E,T,100);
%% ploynomes
S=SAB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SAC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SBC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("T(V,E)");
hold off

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%m�               C(V,E)
##figure
##hold on
##surf(V*1e6,E,C);
##%%contour(V*1e6,E,C,50);
##%%pcolor(V*1e6,E,C);
##%contourf(V*1e6,E,C,100);
##% ploynomes
####n=nab-1;
####S=SAB;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
####n=nac-1;
####S=SAC;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
####n=nbc-1;
####S=SBC;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##zlabel("vitesse du son (m/s)");
##grid()
##set(gca, "fontsize", 40)
##title("C(V,E)");
##hold off
##
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%               G(V,E)
##figure
##hold on
##surf(V*1e6,E,G);
##%%contour(V*1e6,E,G,50);
##%%pcolor(V*1e6,E,G);
##%contourf(V*1e6,E,G,100);
##% ploynomes
####n=nab-1;
####S=SAB;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
####n=nac-1;
####S=SAC;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
####n=nbc-1;
####S=SBC;
####plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##zlabel(" Derivee fondamentale");
##grid()
##set(gca, "fontsize", 40)
##title(" Derivee fondamentale G(V,E)");
##hold off


####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######%               X fraction massique beta
##figure
##hold on
##%surf(V*1e6,E,X);
##%contour(V*1e6,E,X,50);
##%pcolor(V*1e6,E,X);
##contourf(V*1e6,E,Xbeta,[0:0.1:1]);
##% ploynomes
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40)
##title("fraction massique de la phase beta");
##hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Nombre de zones par mailles   
##figure
##hold on
##
##%%%% ploynomes
##n=na-1;
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nb-1;
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nc-1;
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');

%%% grille
##he=E_dis(2)-E_dis(1);
##Edeb=E_dis(1)-he/2;
##Efin=E_dis(ne_dis)+he/2;
##
##hv=V_dis(2)-V_dis(1);
##Vdeb=V_dis(1)-hv/2;
##Vfin=V_dis(nv_dis)+hv/2;
##
##for i=1:nv_dis
##  Vplot=Vdeb+i*hv;
##  plot([Vplot Vplot]*1e6,[Edeb Efin]);
##end
##for i=1:ne_dis+1
##  Eplot=Edeb+i*he;
##  plot([V_dis(1) V_dis(nv_dis)]*1e6,[Eplot Eplot]);
##end
##
##% plot
##surf(V_dis*1e6,E_dis,Nb_ZONE) %,'marker','o');
##%contour(V*1e6,E,P,50);
##%pcolor(V*1e6,E,P);
##%contourf(V*1e6,E,P*1e-9,[0:2:70]);
##%contourf(V*1e6,E,Nb_ZONE,100);
##
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40)
##title("Nb zone par case");
##hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             Zones par mailles   
##figure
##hold on
##surf(V*1e6,E,ZONE_test);
##contourf(V*1e6,E,Nb_ZONE,100);
##%%%%% ploynomes
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40)
##title("Nb zone par case");
##hold off


