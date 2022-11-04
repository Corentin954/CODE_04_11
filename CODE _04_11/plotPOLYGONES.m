clear


% Diagramme V,E
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


Nv=10; Ne=9;
Vmin=112; Vmax=142;
Emin=-1e3; Emax=3.2e5;
hv=(Vmax-Vmin)/Nv;
he=(Emax-Emin)/Ne;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Diagramme (V,E)  
figure
hold on
% Un point 
%V=0.000125781; E=94565.1;
%plot(V*1e6,E,'marker','+','markersize',15);

% Element de la tabulation
##V1=0.0001148; V2=0.0001151; E1=38500; E2=52100;
##Vdis=[ V1 V2]; Edis=[ E1 E1];   plot(Vdis*1e6 ,Edis, "linewidth",2.5);
##Vdis=[ V2 V2]; Edis=[ E1 E2];   plot(Vdis*1e6 ,Edis, "linewidth",2.5);
##Vdis=[ V1 V2]; Edis=[ E2 E2];   plot(Vdis*1e6 ,Edis, "linewidth",2.5);
##Vdis=[ V1 V1]; Edis=[ E1 E2];   plot(Vdis*1e6 ,Edis, "linewidth",2.5);


% Triangle
%plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5);

% Ploynomes
S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'b');
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'r');
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'m');
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');

% Coloriage des zones biphasiques
fill([Va Vb Vc]*1e6,[Ea Eb Ec],[0.7 0.7 0.7]);
S=SAB;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);
S=SAC;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);
S=SBC;
fill(S(:,1)*1e6,S(:,2),[0.8 0.8 0.8]);


% Ploynomes
S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'b');
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'r');
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",2,"color",'m');


%  CADRIALLAGE AU CENTRE
for i=0:Ne
  E=Emin+i*he;
  plot([Vmin, Vmax],[E E],"color",[0.3,0.4,0.5], "linewidth",0.5);
endfor

for i=0:Nv
  V=Vmin+i*hv;
  plot([V V],[Emin, Emax],"color",[0.7,0.7,0.6], "linewidth",0.5);
endfor

### DESSIN DE 3 RECTANGLE EN NOIR (PRESENTATION)
i=2; j=4;
plot(Vmin+ hv*[i, i, i+1, i+1, i], Emin+he*[j, j+1, j+1, j, j],"color",[0,0,0], "linewidth",2);
i=4; j=5;
plot(Vmin+ hv*[i, i, i+1, i+1, i], Emin+he*[j, j+1, j+1, j, j],"color",[0,0,0], "linewidth",2);
i=7; j=3;
plot(Vmin+ hv*[i, i, i+1, i+1, i], Emin+he*[j, j+1, j+1, j, j],"color",[0,0,0], "linewidth",2);


Vp=0.000112312; Ep=100342;
plot([Vp]*1e6 ,[ Ep ],'marker','+', 'linewidth',5) 
Vp=0.000112312; Ep=117111;
plot([Vp]*1e6 ,[ Ep ],'marker','+', 'linewidth',5) 

Vl=0.000111975;  Vr=0.000112648;
El=100101;  Er=134121;
plot([Vl Vr]*1e6,[El El]);
plot([Vl Vr]*1e6,[Er Er]);
plot([Vl Vl]*1e6,[El Er]);
plot([Vr Vr]*1e6,[El Er]);

##
##ylim([Emin Emax]);
##xlim([Vmin Vmax]);
xlabel('\tau (cm^3/kg)');
ylabel('\epsilon (J/kg)');
grid()
set(gca,"fontsize", 40)
title('Diagramme (\tau,\epsilon)');
legend('\beta','\gamma','liquide','zone triphasique','zone biphasique')
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Diagramme (P,T)  
figure
hold on

% Courbes mesures
% BETA/GAMMA
##Pf=14;
##Pplot=PT:(Pf - PT)/100:Pf;
##Tplot=621.36 - 4.88*Pplot - 3.14*Pplot.*Pplot;
##plot(Pplot, Tplot)

% BETA/liquide
##Pf=0;
##Pplot=PT:(Pf - PT)/100:Pf;
##Tplot=505 + 32.1*Pplot - 1.77*Pplot.*Pplot;
##plot(Pplot, Tplot)

% GAMMA/LIQUIDE
##Pf=20;
##Pplot=PT:(Pf - PT)/100:Pf;
##Tplot=581.5 + 59.1*(Pplot - 2.87);
##plot(Pplot, Tplot)


% ploynomes
S=SApt;
plot(S(:,1)*1e-9,S(:,2), "linewidth",2,"color",'blue'); %,'marker','.');
S=SBpt;
plot(S(:,1)*1e-9,S(:,2), "linewidth",2,"color",'magenta'); %,'marker','.');
S=SCpt;
plot(S(:,1)*1e-9,S(:,2), "linewidth",2,"color",'red'); %,'marker','.');

text(1,200,'\beta',"fontsize", 60,"color",'blue');
text(25,702,'\gamma',"fontsize", 60,"color",'magenta');
text(8.5,2600,'liquide',"fontsize", 60,"color",'red');

% Point triple
plot([PT]*1e-9 ,[TT],'marker','+', "linewidth",8,"color",'black')
text(6.5 ,TT,'\leftarrow point triple',"fontsize", 40,"color",'black');

xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,T)");
%legend('\beta','\gamma','liquide')
hold off



