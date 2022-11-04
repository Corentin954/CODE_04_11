clear;

format long;

% On fixe le nombre de points de la discretisation de l'onde de compression isentropiques
N=1e4;

T=2e-7;

u_init=392;

LINE_WIDTH=1;

system("sh compil_onde_composite.sh");  % CHANGER LE .C EN FCT DE TST
system(["a.exe ",num2str(u_init)," ",num2str(N)]);
     
% Onde de compression isentropique
Voci=load("fichiers/Voci.txt");
Eoci=load("fichiers/Eoci.txt");
Poci=load("fichiers/Poci.txt");
Toci=load("fichiers/Toci.txt");
XAoci=load("fichiers/XAoci.txt");
XBoci=load("fichiers/XBoci.txt");
XCoci=load("fichiers/XCoci.txt");
Coci=load("fichiers/Coci.txt");
Soci=load("fichiers/Soci.txt");
Dfoci=load("fichiers/Dfoci.txt");
Uoci=load("fichiers/Uoci.txt");

LAMBDAoci(1,:)=XAoci;
LAMBDAoci(2,:)=XBoci;
LAMBDAoci(3,:)=XCoci;

DS=abs(Soci(N)-Soci(1));
disp(["Variation d'entropie= ",num2str(DS)]);


% Courbe hugoniot
rho0= 7287; V0= 137.230684*1e-6; E0= -0.009371; P0= 0.000000; T0= 300.000000; S0=-4.30968e-06; G0=-0.00937092; Df0=1.15494; C0=2740.55;
rho1= 8076.006891; V1= 123.823569*1e-6; E1= 51443.000837; P1= 7.673987*1e9; T1= 404.147829; D1= 3283.175368; U1= 320.758508; S1= 16.007190; G1= 995194.141539; Df1= 28.698346; C1= 3442.464859;
rho2= 8288.370443; V2= 120.650978*1e-6; E2= 77064.926756; P2= 8.478062*1e9; T2= 354.880625; D2= 3283.175368; U2= 371.265969; S2= 16.016670; G2= 1094267.424242; Df2= -3.031310; C2= 3430.292560;
rho3= 8518.284342; V3= 117.394532*1e-6; E3= 112608.341615; P3= 11.353850*1e9; T3= 388.435466; D3= 3283.175368; U3= 474.570018; S3= 25.068683; G3= 1435750.667237; Df3= 1.047360; C3= 3581.513765;
rho4= 10247.863117; V4= 97.581319*1e-6; E4= 889466.874467; P4= 44.866639*1e9; T4= 1965.674770; D4= 4616.309074; U4= 1333.766759; S4= 305.216384; G4= 4667656.553911; Df4= -5.569062; C4= 4624.846766;
rho4= 10609.008403; V5= 94.259516*1e-6; E5= 1222604.764353; P5= 56.903493*1e9; T5= 2231.310798; D5= 4993.811106; U5= 1563.716582; S5= 383.619791; G5= 5730325.465475; Df5= 25.900309; C5= 4829.418168;
V0*sqrt((P1-P0)/(V0-V1))

lambdaBETA=[1,0,0];

% Soit eps l'epaisseur du choc
eps_choc=1e-8;

x_choc1=V0*sqrt((P1-P0)/(V0-V1))*T - u_init*T;

% Variables qui seront tracées
Xoci=(-Uoci+Coci)*T;

X(1)=0;
X(2:N+1)=flipud(Xoci);
X(N+2)=x_choc1-eps_choc/2;
X(N+3)=x_choc1+eps_choc/2;
X(N+4)=1e-3;

P(1)=Poci(N);
P(2:N+1)=flipud(Poci);
P(N+2)=P1;
P(N+3)=P0;
P(N+4)=P0;

V(1)=Voci(N);
V(2:N+1)=flipud(Voci);
V(N+2)=V1;
V(N+3)=V0;
V(N+4)=V0;

E(1)=Eoci(N);
E(2:N+1)=flipud(Eoci);
E(N+2)=E1;
E(N+3)=E0;
E(N+4)=E0;

U(1)=-Uoci(N);
U(2:N+1)=-flipud(Uoci);
U(N+2)=U1-u_init;
U(N+3)=-u_init;
U(N+4)=-u_init;

C(1)=Coci(N);
C(2:N+1)=flipud(Coci);
C(N+2)=C1;
C(N+3)=C0;
C(N+4)=C0;

S(1)=Soci(N);
S(2:N+1)=flipud(Soci);
S(N+2)=S1;
S(N+3)=S0;
S(N+4)=S0;

Df(1)=Dfoci(N);
Df(2:N+1)=flipud(Dfoci);
Df(N+2)=Df1;
Df(N+3)=Df0;
Df(N+4)=Df0;

for i=1:3
  LAMBDA(i,1)=LAMBDAoci(i,N);
  LAMBDA(i,2:N+1)=LAMBDAoci(i,:);
  LAMBDA(i,N+2)=lambdaBETA(i);
  LAMBDA(i,N+3)=lambdaBETA(i);
  LAMBDA(i,N+4)=lambdaBETA(i);
endfor




%% //////////////////////////////////////////////////
%%               FIGURES


% Diagramme V,E
SA=load("fichiers/pSA.txt");
SB=load("fichiers/pSB.txt");
SC=load("fichiers/pSC.txt");
SAB=load("fichiers/pSAB.txt");
SAC=load("fichiers/pSAC.txt");
SBC=load("fichiers/pSBC.txt");


% hugoniot
Vhug=load("fichiers/Vhug.txt");
Ehug=load("fichiers/Ehug.txt");
[na,~]=size(SA); % de taille NA+1 car SA[NA]=S[0] (code C)
[nb,~]=size(SB);
[nc,~]=size(SC);
[nab,~]=size(SAB); % de taille NBA+1 car SA[NA]=S[0] (code C)
[nac,~]=size(SAC);
[nbc,~]=size(SBC);

% Traingle triple
Va=129.905671875844291e-6;  Ea=72.168190127265447e3;
Vb=127.666112185063994e-6;  Eb=99.480256119236515e3;
Vc=132.538478805442338e-6;  Ec=136.154346228414312e3;



figure
hold on
plot(V*1e6,E,'linewidth',LINE_WIDTH);
%plot(Vhug*1e6,Ehug,"linewidth",1.5);
% Triangle
plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
% Ploynomes
St=SA;
plot(St(:,1)*1e6,St(:,2), "linewidth",0.5,"color",'black','marker','.');
St=SB;
plot(St(:,1)*1e6,St(:,2), "linewidth",0.5,"color",'black','marker','.');
St=SC;
plot(St(:,1)*1e6,St(:,2), "linewidth",0.5,"color",'black','marker','.');
grid();
set(gca, "fontsize", 40);
title(['Energie interne spécinfique, t=',num2str(T),' s, u_unit=',num2str(u_init),' N=',num2str(N)]);
xlabel('V (cm^3/kg')
ylabel("E  (J/Kg)");
xlim([100 180]);
ylim([-2e3 3e5]);
hold off

figure
hold on
plot(X,E,'linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['Energie interne spécifique, t=',num2str(T),' s, u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("E  (J/kg)");
hold off

figure
hold on
plot(X,V,'linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['\tau u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("tau  (m^3/kg)");
hold off

figure
hold on
plot(X,1./V,'linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['\rho u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("rho  (kg/m^3)");
hold off

figure
hold on
plot(X,C,'-o','linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['c u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("c  m/s)");
hold off

figure
hold on
plot(X,S,'-o','linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['S u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("S  J/K/kg)");
hold off

figure
hold on
plot(X,Df,'-o','linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['Dérivee fondamentale u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("G  J/K/kg)");
hold off

figure
hold on
plot(X,U,'linewidth',LINE_WIDTH);
grid();
set(gca, "fontsize", 40);
title(['u u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("u  (m/s)");
hold off

figure
hold on
plot(X,P,'linewidth',LINE_WIDTH);
plot((Coci(N)+Uoci(N))*T,Poci(N),'marker','+','linewidth',4);
grid();
set(gca, "fontsize", 40);
title(['Pression u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("P (Pa)");
hold off


  
figure
hold on
plot(X,LAMBDA(1,:));
plot(X,LAMBDA(2,:));
plot(X,LAMBDA(3,:));
grid();
set(gca, "fontsize", 40);
title(['fraction u_{init}=',num2str(u_init),' N=',num2str(N)]);
xlabel('x (m)')
ylabel("fraction massique");
legend('beta','gamma','liq','location','southwest');
hold off
  

