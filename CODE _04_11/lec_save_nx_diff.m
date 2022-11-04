clear all


% SAVE avec des nx différents

%load('SAVE\ETAIN_tst104_nx_diff_equi_termo_RK5.mat');
%load('SAVE\ETAIN_tst104_nx_diff_tau=1e8_RK1_BBC.mat');
%load('SAVE\ETAIN_tst104_nx_diff_tau=8e9_RK1_BBC.mat');
%load('SAVE\ETAIN_tst104_nx_diff_tau=5e9_RK1_BBC.mat');
%load('SAVE\ETAIN_tst104_nx_diff_tau=5e9_RK5_BBC.mat');


%%  
if  (tst>=2 && tst<=100) || (tst>=113) 

  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(tau(j,1:n-1)*1e6,epsilon(j,1:n-1),Symleg(j,:));
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
  title(['Energie interne spécifique, t=',num2str(T),' s, nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([100 180]);
  ylim([-2e3 3e5]);
  %legend(Lname,'location','northeast');
  hold off

  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),epsilon(j,1:n-1),Symleg(j,1:n-1),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique, t=',num2str(T),' s, nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname,'location','northeast'); 
  hold off

  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),tau(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\tau nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("tau  (m^3/Kg)");
  %legend(Lname,'location','southwest');  
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),1./tau(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname,'location','southwest'); 
  hold off

  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("u  (m/s)");
  %legend(Lname,'location','southwest');
  hold off

  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),p(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Pression nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("P  (Pa)");
  %legend(Lname,'location','southwest');
  hold off
endif

if  tst==100 || tst>=113 
  % fraction massique
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmA(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \beta nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmB(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \gamma nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmC(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase liquide nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmA(j,1:n-1)+fmB(j,1:n-1)+fmC(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques nx=',num2str(nx)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  legend(Lname,'location','southwest');
  hold off

endif

  
%%  solution exacte
if  (tst>=101 && tst<=112) 
  
##  Nana=1e4;
##  Xana=a:(b-a)/Nana:b;
##  [Uana,Pana,RHOana,Eana,lambda_ana]=sol_ana_etain(Xana,T,tst);
  
  figure
  hold on
  plot((1./RHOana)*1e6,Eana,'g*');
  for j=1:ls
    n=Ln(j);
    plot(tau(j,1:n-1)*1e6,epsilon(j,1:n-1),Symleg(j,:));
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
  title(['Energie interne spécifique, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([100 180]);
  ylim([-2e3 3e5]);
  %legend(Lname,'location','northeast');
  hold off

  figure
  hold on
  plot(Xana,Eana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),epsilon(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname,'location','northeast'); 
  hold off

  figure
  hold on
  plot(Xana,1./RHOana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),tau(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\tau , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("tau  (m^3/Kg)");
  %legend(Lname,'location','southwest');  
  hold off
  
  figure
  hold on
  plot(Xana,RHOana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),1./tau(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname,'location','southwest'); 
  hold off

  figure
  hold on
  plot(Xana,Uana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("u  (m/s)");
  %legend(Lname,'location','southwest');
  hold off

  figure
  hold on
  plot(Xana,Pana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),p(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Pression , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("P  (Pa)");
  %legend(Lname,'location','southwest');
  hold off

  % fraction massique
  figure
  hold on
  plot(Xana,lambda_ana(:,1),'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmA(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \beta , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  plot(Xana,lambda_ana(:,2),'-');
  for j=1:ls
    n=Ln(j); 
    plot(Xc(j,1:n-1),fmB(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \gamma , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  plot(Xana,lambda_ana(:,3),'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmC(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase liquide , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  plot(Xana,lambda_ana(:,1)+lambda_ana(:,2)+lambda_ana(:,3),'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),fmA(j,1:n-1)+fmB(j,1:n-1)+fmC(j,1:n-1),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques']);
  xlabel('x (m)')
  ylabel("fraction massique");
  legend(Lname,'location','southwest');
  hold off

endif

  

%% Cas du test de Sod
if tst==0 || tst==1
    
##  Nana=1e4;
##  Xdana=a:(b-a)/Nana:b;
##  for i=1:Nana
##    Xcana(i)=(Xdana(i)+Xdana(i+1))/2;
##  endfor 
##  % sol analytique
##  [W]=sol_ana_PU(Xcana,xdis,T,wL,wR,p12,u12,gamma);
##  RHOana=W(1,:);
##  Pana=W(3,:);
##  [W]=sol_ana_PU(Xdana,xdis,T,wL,wR,p12,u12,gamma);
##  Uana=W(2,:);

  figure
  hold on
  plot(Xcana,Pana./((gamma-1)*RHOana),'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),epsilon(j,1:n-1),Symleg(j,1:n-1),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique,  t=',num2str(T),' s , nx=',num2str(nx)]);
  xlabel("m");
  ylabel(["E  (J/kg)"]);
  %legend(Lname,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,RHOana,'-','linewidth',LINE_WIDTH);
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),1./tau(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
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
    n=Ln(j);
    plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx)]);
  legend(Lname,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,Pana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),p(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
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


