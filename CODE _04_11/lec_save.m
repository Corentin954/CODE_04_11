clear

% hugoniot
Vhug=load("fichiers/Vhug.txt");
Ehug=load("fichiers/Ehug.txt");
Phug=load("fichiers/Phug.txt");
Thug=load("fichiers/Thug.txt");
Dhug=load("fichiers/Dhug.txt");
Uhug=load("fichiers/Uhug.txt");

XAhug=load("fichiers/XAhug.txt");
XBhug=load("fichiers/XBhug.txt");
XChug=load("fichiers/XChug.txt");

% Première hugoniot
Vhug=load("fichiers/Vhug.txt");
Vhug1=load("fichiers/Vhug1.txt");
Ehug1=load("fichiers/Ehug1.txt");
Phug1=load("fichiers/Phug1.txt");
Thug1=load("fichiers/Thug1.txt");
Dhug1=load("fichiers/Dhug1.txt");
Uhug1=load("fichiers/Uhug1.txt");

XAhug1=load("fichiers/XAhug1.txt");
XBhug1=load("fichiers/XBhug1.txt");
XChug1=load("fichiers/XChug1.txt");


V0= 137.230419*1e-6; E0= -0.286078; P0= 0.000101*1e9; T0= 300.000000;
V1= 123.823506*1e-6; E1= 51442.887611; P1= 7.674026*1e9; T1= 404.145405; D1= 3283.180446; u1= 320.754782;
V2= 120.650923*1e-6; E2= 77064.866965; P2= 8.478099*1e9; T2= 354.878407;
V3= 117.394432*1e-6; E3= 112608.972466; P3= 11.353935*1e9; T3= 388.433657; D3= 3283.180446; u3= 474.567695;
V4= 97.581086*1e-6; E4= 889478.598190; P4= 44.867180*1e9; T4= 1965.687527; D4= 4616.324687; u4= 1333.772744;
V5= 94.259292*1e-6; E5= 1222618.419852; P5= 56.904095*1e9; T5= 2231.323266; D5= 4993.825841; u5= 1563.722707;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAUVEGARDE DES SIMULATIONS DU RAPPORT
%load('BIZ_19_08.mat');
%load('SAVE\ETAIN_PT1.mat')
%load('SAVE\ETAIN_X_GAMMA=08.mat');
%load('SAVE\ETAIN_X_GAMMA=08_n=800.mat');
%load('SAVE\ETAIN_X_GAMMA=08_n=300_TAU=5e9.mat');
%load('SAVE\ETAIN_X_GAMMA=08_n=800_TAU=5e9.mat');
%load('SAVE\ETAIN_RK5_n=800_TAU=5e9.mat');
%load('SAVE\ETAIN_God_n=300_TAU_diff.mat');
%load('SAVE\ETAIN_tst104_RK5_n=300_TAU_diff.mat');
%load('SAVE\ETAIN_tst112_RK5_n=300_TAU_diff.mat');
%load('SAVE\ETAIN_tst104_BBC_nx_diff_TAU_5e9.mat');


dhe=(4.3e3/3)/(6.77e4/5)
dhv=(20*5)/(67*3)

size_line=1.3;

PLOT_NRJ=0;

%%  
if tst!=0 %0
  
  figure
  hold on
  for j=1:ls
    plot(tau(j,:)*1e6,epsilon(j,:),Symleg(j,:),size_line);
  endfor
  plot(Vhug*1e6,Ehug,"linewidth",1.5);
  % Triangle
  plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
  % Ploynomes
  S=SA;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  S=SB;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  S=SC;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique, t=',num2str(T),' s , nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg)')
  ylabel("E  (J/Kg)");
  xlim([100 180]);
  ylim([-2e3 3e5]);
  %legend(Lname,'location','northeast');
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(tau(j,:)*1e6,p(j,:)*1e-9,Symleg(j,:),"linewidth",size_line);
  endfor
  % points
  plot(V0*1e6,P0*1e-9,'marker','+', "linewidth",5);
  plot(V1*1e6,P1*1e-9,'marker','+', "linewidth",5);
  plot(V2*1e6,P2*1e-9,'marker','+', "linewidth",5);
  plot(V3*1e6,P3*1e-9,'marker','+', "linewidth",5);
  plot(V4*1e6,P4*1e-9,'marker','+', "linewidth",5);
  plot(V5*1e6,P5*1e-9,'marker','+', "linewidth",5);
  
  VM=111.82e-6; PM=21.2811e9;
  VF=106.2e-6; PF=26.117*1e9;
  plot((VM)*1e6,(PM)*1e-9,'marker','+', "linewidth",5);
  text((VM)*1e6,(PM)*1e-9,'  M','color','r',"fontsize", 45)
  plot((VF)*1e6,(PF)*1e-9,'marker','+', "linewidth",5);
  text((VF)*1e6,(PF)*1e-9,'  F','color','r',"fontsize", 45)
  %isentrope
  plot(Vhug*1e6,Phug*1e-9,"linewidth",1);
  plot(Vhug1*1e6,Phug1*1e-9,"-.","linewidth",1);
  % droite P0,P1,P3
  vplot=V0:(V3-V0)/100:V3-1e-5;
  plot(vplot*1e6,(P0+(P1-P0)*(vplot-V0)/(V1-V0))*1e-9,"--","linewidth",0.8);
  % droite P0,PA,PF
  vplot=V0:(V4-V0)/100:V4-1e-5; % V4 limite de V
  plot(vplot*1e6,(P0+(PM-P0)*(vplot-V0)/(VM-V0))*1e-9,"--","linewidth",0.8);
  
  %isentrope
  text(V0*1e6,P0*1e-9,' 0',"fontsize", 45);
  text(V1*1e6,P1*1e-9,' 1',"fontsize", 45);
  text(V2*1e6,P2*1e-9,' 2',"fontsize", 45);
  text(V3*1e6,P3*1e-9,' 3',"fontsize", 45);
  text(V4*1e6,P4*1e-9,' 4',"fontsize", 45);
  text(V5*1e6,P5*1e-9,' 5',"fontsize", 45);
  
  ylabel("P (GPa)");
  xlabel("V (cm^3/kg)");
  grid()
  set(gca, "fontsize", 40)
  title("Hugoniot (P,V)");
  hold off

  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:),"linewidth",size_line);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique ,  t=',num2str(T),' s , nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  ylabel("E (J/Kg)");
  %legend(Lname,'location','northeast'); 
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),tau(j,:),Symleg(j,:),"linewidth",size_line);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Volume spécifique nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  ylabel("tau (m^3/Kg)");
  %legend(Lname,'location','southwest');  
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),1./tau(j,:),Symleg(j,:),"linewidth",size_line);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Densité nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname,'location','southwest'); 
  hold off

  
  figure
  hold on
  for j=1:ls
    plot(Xd(j,:),u(j,:),Symleg(j,:),"linewidth",size_line);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Vitesse nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  ylabel("u (m/s)");
  %legend(Lname,'location','southwest');
  hold off

  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:),"linewidth",size_line);
  endfor
  xM=0.0006195; PM=21.2811e9;
  plot([xM],[PM],"color",'r',"linewidth",5,"marker",'+');
  text([xM],[PM],'  M','color','r',"fontsize", 45)
  grid();
  set(gca, "fontsize", 40);
  title(['Pression nx=',num2str(nx),', CFL= ',num2str(numCFL)]);
  xlabel('m')
  ylabel("P  (Pa)");
  %legend(Lname,'location','southwest');
  hold off

endif
 
if tst>=100
  % fraction massique
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \beta nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmB(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \gamma nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('m')
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase liquide nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:)+fmB(j,:)+fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques nx=',num2str(nx)]);
  legend(Lname,'location','southwest');
  hold off

endif

  
  

%% Cas du test de Sod
if tst==0 
  
  Lname=['solution analytique';Lname];
  
  Nana=1e4;
  Xdana=a:(b-a)/Nana:b;
  for i=1:Nana
    Xcana(i)=(Xdana(i)+Xdana(i+1))/2;
  endfor 
  % sol analytique
  [W]=sol_ana_PU(Xcana,xdis,T,wL,wR,p12,u12);
  RHOana=W(1,:);
  Pana=W(3,:);
  [W]=sol_ana_PU(Xdana,xdis,T,wL,wR,p12,u12);
  Uana=W(2,:);

  figure
  hold on
  plot(Xcana,Pana./((gamma-1)*RHOana),'-');
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique ,  t=0.2s , nx=',num2str(nx)]);
  xlabel("m");
  ylabel(["E  (J/kg)"]);
  %legend(Lname,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,RHOana,'-','linewidth',1.5);
  for j=1:ls
    plot(Xc(j,:),1./tau(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx)]);
  %legend(Lname,'location','northwest');
  hold off

  figure
  hold on
  plot(Xdana,Uana,'-');
  for j=1:ls
    plot(Xd(j,:),u(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx)]);
  %legend(Lname,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,Pana,'-');
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['p nx=',num2str(nx)]);
  %legend(Lname,'location','northwest');
  hold off

endif

  
  
  
  
  
if PLOT_NRJ==1
  % Plot de l'energie du système au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Etot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['energie totale Bizarrium']);
  legend(Lname,'location','southwest');
  hold off
  
  
  % Plot de l'impulsion au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,IMPUL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['impulsion']);
  legend(Lname,'location','southwest');
  hold off
  
  % Plot de la masse totale au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Mtot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['masse totale']);
  legend(Lname(1,:),'location','southwest');
  hold off
  
  % Plot du volume total au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,VOL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['volume total']);
  legend(Lname(1,:),'location','southwest');
  hold off
  
endif

