clear all


%///////////////             RUNGE-KUTTA               ///////////////////
%///////////////    Plot  schema + sol analytique        ///////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix de la loi d'etat et du cas test
%  tst : 0 (gaz parfait)   2
tst = 102; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deux_plaques=floor(tst/100);

if tst==0  % Sod
  T=0.2; % pas en ligne de commande
elseif tst==1  % LeBlanc
  T=6;
  %a=0; b=9;
elseif tst==2  % Bizarrium
  T=80e-6; % pas en ligne de commande
elseif tst==3  % Onde acoustique
  T=1/2; % pas en ligne de commande
elseif tst==4  % Sod sym�tris�
  %T=0.2;
elseif tst==5  % Woodward
  %T=0.038;
  %a=0; b=1;
elseif tst==6  % Shu-Osher
  %T=1.8;
  %a=-5; b=5;
elseif tst==100   % etain
  T=2e-7;
elseif deux_plaques==1  % etain
  T=2e-7;
  a=-1e-3; b=1e-3; 
elseif tst==200 || tst==201 || tst==202  % etain
  T=1e-7;
else
 disp(["erreur valeur tst= (Script Octave)",num2str(tst)]);
end


if tst==0 || tst==1
  if(tst==0)  % test de Sod
    gamma=7/5; 
    rhoG=1.0; rhoD=0.125;
    uG=0; uD=0;
    pG=1.0; pD=0.1;
    a=0; b=1;
    xdis=(b+a)/2;
  elseif(tst==1)   % test de LeBlanc
    gamma=5/3; 
    rhoG=1; rhoD=1e-3;
    uG=0; uD=0;
    pG=(gamma-1)/10; pD=(gamma-1)*1e-12;
    a=0; b=9;
    xdis=(b+a)/3;
  else 
    disp("Erreur");
  endif
  
  wL=[rhoG; uG; pG];
  wR=[rhoD; uD; pD];
  
  % Calcul de P et U
  system("sh compil_PU.sh");  % CHANGER LE .C EN FCT DE TST
  system(["./a.out ",num2str(tst)]);
  resu=load("PU.txt");
  p12=resu(1)
  u12=resu(2)
  
  Nana=1e4;
  Xdana=a:(b-a)/Nana:b;
  for i=1:Nana
    Xcana(i)=(Xdana(i)+Xdana(i+1))/2;
  endfor 
  % sol analytique
  [W]=sol_ana_PU(Xcana,xdis,T,wL,wR,p12,u12,gamma);
  RHOana=W(1,:);
  Pana=W(3,:);
  [W]=sol_ana_PU(Xdana,xdis,T,wL,wR,p12,u12,gamma);
  Uana=W(2,:);
endif
if tst==3
  gamma=7/5; 
  fac_eps=1e-8;
  kappa=2*pi;
  rho0=1.0; p0=5./7;
  
  a=0; b=1;
  Nana=2e3;
  Xdana=a:(b-a)/Nana:b;
  for i=1:Nana
    Xcana(i)=(Xdana(i)+Xdana(i+1))/2;
  endfor 
  
  rho_bar=@(a,b)  rho0 - (fac_eps/(kappa*(b-a)))*(cos(kappa*b) - cos(kappa*a));
  u_bar=@(a,b)  -(fac_eps/(kappa*(b-a)))*(cos(kappa*b) - cos(kappa*a));
  p_ana=@(x)  p0 + fac_eps*sin(kappa*x);
  
  for i=1:Nana
    xg=Xdana(i); xd=Xdana(i+1);
    RHOana(i)=rho_bar(xg,xd);
    Pana(i)=p_ana((xd+xg)/2);
    Uana(i)=u_bar(xg,xd);
  endfor  
endif



% =======================================================================
% Creation de l'executable
system("sh compil_main_multi.sh");




%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
##Lname=['BBC RK2av, q o1 Q+L, Tau=5e-9, nx=200';
##       'GoHy 2 MinMod, Tau=5e-9, nx=200';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix, Tau=5e-9, nx=200';
##       'BBC RK2av, q o1 Q+L, Tau=5e-9, nx=400';
##       'GoHy 2 MinMod, Tau=5e-9, nx=400';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix, Tau=5e-9, nx=400';
##       'BBC RK2av, q o1 Q+L, Tau=5e-9, nx=800';
##       'GoHy 2 MinMod, Tau=5e-9, nx=800';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix, Tau=5e-9, nx=800';
##       'BBC RK2av, q o1 Q+L, Tau=5e-9, nx=1600';
##       'GoHy 2 MinMod, Tau=5e-9, nx=1600';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix, Tau=5e-9, nx=1600';
##       ];
##
##Lsch=[202,014,105,202,014,105,202,014,105,202,014,105];  % scheme
##Lso=[-1,3,5,-1,5,5,-1,2,5,-1,2,5];           % spatial order
##Lz=[-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1];            % kinetic energy fix
##Lq=(004)*[1,1,1,1,1,1,1,1,1,1,1,1];        % artificial viscosity
##Ldpi=[0,0,0,0,0,0,0,0,0,0,0,0];          % dPI vs delta PI bar
##LTtau=5e-9*[1,1,1,1,1,1,1,1,1,1,1,1];
##Lcin_phase=1*[1,1,1,1,1,1,1,1,1,1,1,1];
##%Lnx=[200, 200, 200, 400, 400, 400, 800, 800, 800, 1600, 1600, 1600];
##Lnx=1e-1*[200, 200, 200, 400, 400, 400, 800, 800, 800, 1600, 1600, 1600];


##Lname=['God acous o1, eq thermo nx=200';
##       'God acous o1, eq thermo nx=400';
##       'God acous o1, eq thermo nx=800';
##       'God acous o1, eq thermo nx=1600';
##       'God acous o1, eq thermo nx=3200';
##];
##
##Lsch=[002,002,002,002,002];  % scheme
##Lso=[-1,-1,-1,-1,-1];           % spatial order
##Lz=[-1,-1,-1,-1,-1];            % kinetic energy fix
##Lq=[-1,-1,-1,-1,-1];        % artificial viscosity
##Ldpi=[0,0,0,0,0];          % dPI vs delta PI bar
##LTtau=0*[1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1];
##Lnx=[200, 400, 800, 1600, 3200];

##Lname=['God acous o1, eq thermo nx=400';
##       'God acous o1, eq thermo nx=800';
##       'God acous o1, eq thermo nx=1600';
##       'God acous o1, eq thermo nx=3200';
##       'God acous o1, eq thermo nx=6400';
##];
##
##Lsch=[002,002,002,002,002];  % scheme
##Lso=[-1,-1,-1,-1,-1];         % spatial order
##Lz=[-1,-1,-1,-1,-1];          % kinetic energy fix
##Lq=[-1,-1,-1,-1,-1];          % artificial viscosity
##Ldpi=[0,0,0,0,0];              % dPI vs delta PI bar
##LTtau=0*[1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1];
##Lnx=[400, 800, 1600, 3200, 6400];


Lname=['God acous o1, eq thermo nx=200';
       'God acous o1, eq thermo nx=400';
       'God acous o1, eq thermo nx=800';
       'God acous o1, eq thermo nx=1600';
       'God acous o1, eq thermo nx=3200';
];

Lsch=[002,002,002,002,002];  % scheme
Lso=[-1,-1,-1,-1,-1];         % spatial order
Lz=[-1,-1,-1,-1,-1];          % kinetic energy fix
Lq=[-1,-1,-1,-1,-1];          % artificial viscosity
Ldpi=[0,0,0,0,0];              % dPI vs delta PI bar
LTtau=0*[1,1,1,1,1];
Lcin_phase=0*[1,1,1,1,1];
Lnx=[200, 400, 800, 1600, 3200];


##Lname=['RK1, q o1 Q+L, so 2, eq thermo, nx=200'; 
##       'RK1, q o1 Q+L, so 2, eq thermo, nx=400'; 
##       ];
##
##Lsch=[100,100];  % scheme
##Lso=[2,2];           % spatial order
##Lz=[1,1];            % kinetic energy fix
##Lq=(004)*[1,1];        % artificial viscosity
##Ldpi=[0,0];          % dPI vs delta PI bar
##LTtau=1e-9*[1,1];
##Lcin_phase=0*[0,0];
##Lnx=[200,400];




ls=length(Lsch);

tabGRILLE=(floor(Lsch/100)>=1);

% CFL
numCFL=0.5;
numCFL=0.9;
%numCFL=0.3;
 
Ntau=15;

fac_dt=1e-2;

epsilon_newton=1e-13;

% Coeff de pseudo
Cq=1.5;
%Cq=0;
Cl=0.15;
%Cl=0;
%Cl=0.29;

% affichage
aff=0;
PLOT_NRJ=0;
LINE_WIDTH=0.5;
PLOT=0;
PLOT_ITER=0;
PLOT_ERR=1;

sch_cin_phase=0;

% Multi fluide
tst_multi=-1;

% Max iter
   Nmax= 2e6;

% symoble pour la legende
Symleg=['.-'; '-.'; '+-';'--', '*-'; 'x-'; 'o-'; 'd-'; 's-';'h-';'v-';':';'^-';'v-' ];


Etot=zeros(ls, 1e5);

maxN=max(Lnx)+10;

Xd = zeros(ls,maxN); % X d�centr�e (fronti�res)
Xc = zeros(ls,maxN);
tau = zeros(ls,maxN);
u = zeros(ls,maxN);
e = zeros(ls,maxN);
p = zeros(ls,maxN);
epsilon = zeros(ls,maxN);


fmA=zeros(ls,maxN);
fmB=zeros(ls,maxN);
fmC=zeros(ls,maxN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ls
  com=["./a.out  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(Lnx(j))," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)," -j ",num2str(epsilon_newton)];
  disp([' ']);
  disp([Lname(j,:)]);
  disp([' ']);
  disp([com]);
  system(com)
  
  resu=load("sol.txt");
  [~,n]=size(resu);
  Ln(j)=n;
  %size(resu)
  Xd(j,1:n) = resu(1,:); % X d�centr�e (fronti�res)
  Xc(j,1:n-1) = resu(2,1:n-1);
  rho(j,1:n-1) = resu(3,1:n-1);
  tau(j,1:n-1) = resu(4,1:n-1);
  u(j,1:n) = resu(5,:);
  e(j,1:n-1) = resu(6,1:n-1);  
  epsilon(j,1:n-1) = resu(7,1:n-1);
  p(j,1:n-1) = resu(8,1:n-1);
  tempe(j,1:n-1) = resu(9,1:n-1);
  c(j,1:n-1) = resu(10,1:n-1);
  G(j,1:n-1) = resu(11,1:n-1);
  Sentropie(j,1:n-1) = resu(12,1:n-1);
  
  % energie du syt�me
  NRJ=load("NRJ.txt");
  Nbiter(j)=length(NRJ(1,:));
  Etot(j,1:Nbiter(j))=NRJ(1,:);
  IMPUL(j,1:Nbiter(j))=NRJ(2,:);
  Mtot(j,1:Nbiter(j))=NRJ(3,:);
  VOL(j,1:Nbiter(j))=NRJ(4,:);
  
  resu=load("frac_mass.txt");
  [~,n]=size(resu);
  fmA(j,1:n-1)=resu(1,1:n-1);
  fmB(j,1:n-1)=resu(2,1:n-1);
  fmC(j,1:n-1)=resu(3,1:n-1);
  
  resu=load("nb_iter.txt");
  tabNiter(j,1:n-1)=resu(1,1:n-1);
  tabCRITERE(j,1:n-1)=resu(2,1:n-1);
endfor


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


if PLOT==1
%%  
if  tst==2 || (tst>=4 && tst<=100) || (tst>=113) 

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
  title(['Energie interne sp�cifique, t=',num2str(T),' s, CFL= ',num2str(numCFL)]);
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
    plot(Xc(j,1:n-1),epsilon(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp�cifique, t=',num2str(T),' s, CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname,'location','northeast'); 
  hold off

  figure
  hold on
  for j=1:ls
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
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),rho(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
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
  for j=1:ls
    n=Ln(j);
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,1:n-1),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille d�cal�e
      plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
    endif
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
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),c(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Celerite du son, CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("c  (m/s)");
  legend(Lname,'location','southwest');
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
  title(['phase \beta, CFL= ',num2str(numCFL)]);
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
  title(['phase \gamma , CFL= ',num2str(numCFL)]);
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
  title(['phase liquide , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname,'location','southwest');
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:)+fmB(j,:)+fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques']);
  xlabel('x (m)')
  ylabel("fraction massique");
  legend(Lname,'location','southwest');
  hold off

endif

  
%%  solution exacte
if  (tst>=101 && tst<=112) 
  
  [Uana,Pana,RHOana,Eana,lambda_ana,Xana, Xchoc1, Xchoc2]=sol_ana_etain(a,b,T,tst);
  
  Lname_ana=['solution exacte';Lname];

  figure
  hold on
  plot((1./RHOana)*1e6,Eana,'r+',10);
  for j=1:ls
    n=Ln(j);
    plot(tau(j,1:n-1)*1e6,epsilon(j,1:n-1),Symleg(j,:),1.5);
  endfor
  plot(Vhug*1e6,Ehug,"linewidth",0.5);
  % Triangle
  plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r',0.5);
  % Ploynomes
  S=SA;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  S=SB;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  S=SC;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp�cinfique, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([110 140]);
  ylim([-1e2 1.4e5]);
  %legend(Lname,'location','northeast');
  hold off
  
%%% plot uniquement les points entre les deux chocs a droite
  
  for j=1:ls
    k=1;
    n=Ln(j);
    
    test=(Xc(j,1:n-1)<=Xchoc1).*(Xc(j,1:n-1)>=Xchoc2);
    tau_P1=test.*tau(j,1:n-1);
    epsilon_P1=test.*epsilon(j,1:n-1);
    for i=1:n-1
      if(tau_P1(i)!=0 && epsilon_P1(i)!=0)
        tau_P1_plot(j,k)=tau_P1(i);
        epsilon_P1_plot(j,k)=epsilon_P1(i);
        k++;
      endif
    endfor
    Lk(j)=k-1;
  endfor
  
  figure("name","diag VE 1 hug")
  hold on
  plot((1./RHOana)*1e6,Eana,'r*',5);
  for j=1:ls
    k=Lk(j);
    plot(tau_P1_plot(j,1:k)*1e6,epsilon_P1_plot(j,1:k),Symleg(j,:),1.5);
  endfor
  plot(Vhug*1e6,Ehug,"linewidth",0.5);
  % Triangle
  plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r',0.5);
  % Ploynomes
  S=SA;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  S=SB;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  S=SC;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([120 130]);
  ylim([-1e2 1e5]);
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
  title(['Energie interne sp�cifique, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname_ana,'location','northeast');
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),c(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Celerite du son, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("c  (m/s)");
  %legend(Lname_ana,'location','northeast');
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),G(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Derivee fondamentale, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("G ");
  %legend(Lname,'location','northeast'); 
  hold off
  
  figure
  hold on
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),Sentropie(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Entropie, t=',num2str(T),' s , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("S  ");
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
  %legend(Lname_ana,'location','southwest');  
  hold off
  
  figure
  hold on
  plot(Xana,RHOana,'-');
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),rho(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname_ana,'location','southwest'); 
  hold off

  figure
  hold on
  plot(Xana,Uana,'-');
  for j=1:ls
    n=Ln(j);
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,1:n-1),u(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille d�cal�e
      plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
    end
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("u  (m/s)");
  %legend(Lname_ana,'location','southwest');
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
  %legend(Lname_ana,'location','southwest');
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
  %legend(Lname_ana,'location','southwest');
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
  %legend(Lname_ana,'location','southwest');
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
  %legend(Lname_ana,'location','southwest');
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
  legend(Lname_ana,'location','southwest');
  hold off

endif

  

%% Cas du test de Sod
if tst==0 || tst==1 || tst==3
  
  Lname_ana=['solution analytique';Lname];
  
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
    plot(Xc(j,1:n-1),epsilon(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp�cifique,  t=',num2str(T),' s ']);
  xlabel("m");
  ylabel(["E  (J/kg)"]);
  %legend(Lnama_ana,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,RHOana,'-','linewidth',LINE_WIDTH);
  for j=1:ls
    n=Ln(j);
    plot(Xc(j,1:n-1),rho(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho ']);
  %legend(Lnama_ana,'location','northwest');
  hold off

  figure
  hold on
  if(tst==0 || tst==1)
    plot(Xdana,Uana,'-');
  elseif(tst==3)
    plot(Xcana,Uana,'-');
  endif
  for j=1:ls
    n=Ln(j);
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,1:n-1),u(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille decentr�e
      plot(Xd(j,1:n),u(j,1:n),Symleg(j,:),'linewidth',LINE_WIDTH);
    end
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u ']);
  legend(Lnama_ana,'location','northwest');
  hold off

  figure
  hold on
  plot(Xcana,Pana,'-');
  for j=1:ls
    n=Ln(j)
    plot(Xc(j,1:n-1),p(j,1:n-1),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['p ']);
  %legend(Lnama_ana,'location','northwest');
  hold off

endif

endif % indice PLOT  
  
  
  
  
if PLOT_NRJ==1
  % Plot de l'energie du syst�me au cours du temps
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
  legend(Lname,'location','southwest');
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
  legend(Lname,'location','southwest');
  hold off
  
endif



if PLOT_ITER 

  for j=1:ls
    n=Ln(j);
    figure("name",["Nb iteration ",Lname(j,:)])
    hold on
    bar(tabNiter(j,1:n-1));
    set(gca, "fontsize", 40);
    ylabel("Nombre d'itération par maille au dernier cycle en temps");
    xlabel('Numéro de la maille');
    title(["Nb iter | ",Lname(j,:)]);
    grid()
    hold off
  endfor

  for j=1:ls
    n=Ln(j);
    figure
    hold on
    plot(tabCRITERE(j,1:n-1),Symleg(j,:));
    set(gca, "fontsize", 40);
    ylabel('Critere max sur les cycles en temps par maille');
    xlabel('Numéro de la maille');
    title(["Critere max sur les cycles par maille | ",Lname(j,:)]);
    grid()
    hold off
  endfor
  
endif


if PLOT_ERR
%%   % Calcul de l'erreur sur P sur l'étain
if (tst>=101 && tst<=112) 
  cas=2;
  [ErrL1,ErrL2,ErrL1_red,ErrL2_red, ErrP1, ErrPf]=err_sol_etain(cas,T,tst,Ln,Xd,Xc,p);
endif

% prise d'ordre
po=polyfit(log10(Ln),log10(ErrL1),1);
ordreL1=-po(1)
po=polyfit(log10(Ln),log10(ErrL2),1);
ordreL2=-po(1)

if  (tst>=101 && tst<=112) % etain

  figure
  hold on
  plot(log10(Ln),log10(ErrL1),'--','marker','+',"markersize", 10);
  plot(log10(Ln),-log10(Ln) + (log10(ErrL1(1))+log10(Ln(1))));
  plot(log10(Ln),log10(ErrL1_red),'--','marker','o',"markersize", 10);
  xlabel('log nombre de mailles');
  ylabel('log erreur L1');
  title(['Erreur L1 log-log sur la pression | ordre exp= ',num2str(ordreL1,"%.2f")]);
  legend('Erreur','droite ordre 1','erreur reduite');
  set(gca, "fontsize", 40);
  grid()
  hold off

  figure
  hold on
  plot(log10(Ln),log10(ErrL2),'--','marker','+',"markersize", 10);
  plot(log10(Ln),-log10(Ln) + (log10(ErrL2(1))+log10(Ln(1))));
  plot(log10(Ln),log10(ErrL2_red),'--','marker','o',"markersize", 10);
  xlabel('log nombre de mailles');
  ylabel('log erreur L2');
  title(['Erreur L2 log-log sur la pression | ordre exp= ',num2str(ordreL2,"%.2f")]);
  legend('Erreur','Droite ordre 1','erreur reduite');
  set(gca, "fontsize", 40);
  grid()
  hold off

endif

if tst==102
  
  ErrP1_abs=ErrP1(:,1);
  ErrP1_rel=ErrP1(:,2);
  
  ErrPf_abs=ErrPf(:,1);
  ErrPf_rel=ErrPf(:,2);
  
  % prise d'ordre
  po=polyfit(log10(Ln),log10(ErrP1_abs),1);
  ordreP1=-po(1)
  otracP1=round(ordreP1);
  po=polyfit(log10(Ln),log10(ErrPf_abs),1);
  ordrePf=-po(1)
  otracPf=round(ordrePf);

  figure
  hold on
  plot(log10(Ln),log10(ErrP1_abs),'--','marker','+',"markersize", 10);
  plot(log10(Ln), otracP1*log10(Ln/Ln(1)) + log10(ErrP1_abs(1)) ));
  plot(log10(Ln),log10(ErrPf_abs),'--','marker','o',"markersize", 10);
  plot(log10(Ln), otracPf*log10(Ln/Ln(1)) + log10(ErrPf_abs(1)) );
  xlabel('log nombre de mailles');
  ylabel('log erreur L1');
  title(['Erreur log-log absolue sur la pression sur P1 et etat final ']);
  legend(['Erreur P1 | ordre= ',num2str(ordreP1)],['droite ordre ',num2str(otracP1)],['Errreur Etat final | ordre= ',num2str(ordrePf)],['droite ordre ',num2str(otracPf)]);
  set(gca, "fontsize", 40);
  grid()
  hold off
  
  % prise d'ordre
  po=polyfit(log10(Ln),log10(ErrP1_rel),1);
  ordreP1=-po(1)
  otracP1=round(ordreP1);
  po=polyfit(log10(Ln),log10(ErrPf_rel),1);
  ordrePf=-po(1)
  otracPf=round(ordrePf);
  
  
  figure
  hold on
  plot(log10(Ln),log10(ErrP1_rel),'--','marker','+',"markersize", 10);
  plot(log10(Ln), otracP1*log10(Ln/Ln(1)) + log10(ErrP1_rel(1)) ));
  plot(log10(Ln),log10(ErrPf_rel),'--','marker','o',"markersize", 10);
  plot(log10(Ln), otracPf*log10(Ln/Ln(1)) + log10(ErrPf_rel(1)) );
  xlabel('log nombre de mailles');
  ylabel('log erreur L1');
  title(['Erreur log-log relative sur la pression sur P1 et etat final ']);
  legend(['Erreur P1 | ordre= ',num2str(ordreP1)],['droite ordre ',num2str(otracP1)],['Errreur Etat final | ordre= ',num2str(ordrePf)],['droite ordre ',num2str(otracPf)]);
  set(gca, "fontsize", 40);
  grid()
  hold off
  
endif
  

endif


%save('SAVE\ETAIN_tst104_nx_diff_tau=5e9_RK5_BBC_1600.mat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('BIZ_19_08.mat');