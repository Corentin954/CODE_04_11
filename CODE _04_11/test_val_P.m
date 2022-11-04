clear


P=load("Pbench.txt");
T=load("Tbench.txt");
V=load("Vbench.txt");
E=load("Ebench.txt");
[~,ne]=size(Emaps);
[~,nv]=size(Vmaps);

P=MAPS(:,1);
T=MAPS(:,2);

P=reshape(P,ne,nv);
T=reshape(T,ne,nv);


SA=load("pSA.txt");
SB=load("pSB.txt");
SC=load("pSC.txt");
SAB=load("pSAB.txt");
SAC=load("pSAC.txt");
SBC=load("pSBC.txt");

[na,~]=size(SA); % de taille NA+1 car SA[NA]=S[0] (code C)
[nb,~]=size(SB);
[nc,~]=size(SC);
[nab,~]=size(SAB); % de taille NA+1 car SA[NA]=S[0] (code C)
[nac,~]=size(SAC);
[nbc,~]=size(SBC);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             P(V,E)  
figure
hold on
surf(Vmaps*1e6,Emaps,P*1e-9);
%contour(Vmaps*1e6,Emaps,P,50);
%pcolor(Vmaps*1e6,Emaps,P);
%contourf(Vmaps*1e6,Emaps,P*1e-9,[0:2:70]);
% ploynomes
n=nab-1;
S=SAB;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
n=nac-1;
S=SAC;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
n=nbc-1;
S=SBC;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("P(V,E)");
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               T(V,E)
figure
hold on
surf(Vmaps*1e6,Emaps,T);
%contour(Vmaps*1e6,Emaps,T,50);
%pcolor(Vmaps*1e6,Emaps,T);
%contourf(Vmaps*1e6,Emaps,T,[0:2:70]);
% ploynomes
n=nab-1;
S=SAB;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
n=nac-1;
S=SAC;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
n=nbc-1;
S=SBC;
for i=1:n
  plot([S(i,1) S(i+1,1)]*1e6, [S(i,2) S(i+1,2)], "linewidth",1.5,"color",'black');
end
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("P(V,E)");
hold off


