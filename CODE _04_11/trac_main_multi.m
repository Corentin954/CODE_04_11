clear all

%///////////////    Plot  schema + sol analytique        ///////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix de la loi d'etat et du cas test
%  tst : 0 (gaz parfait)   2
tst = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_init=361.3;
deux_plaques=floor(tst/100);

if tst==0  % Sod
  T=0.2; % pas en ligne de commande
  nx=400;
  elseif tst==1  % LeBlanc
    T=6;
    nx=1000;
    %a=0; b=9;
  elseif tst==2  % Bizarrium
    T=80e-6; % pas en ligne de commande
    nx=1000;
  elseif tst==3  % Onde acoustique
    T=1/2; % pas en ligne de commande
    nx=16;
  elseif tst==4  % Sod symétrisé
    T=0.2;
    nx=200;
  elseif tst==5  % Woodward
    %T=0.038;
    nx=300;
    %a=0; b=1;
  elseif tst==6  % Shu-Osher
    %T=1.8;
    nx=300;
    %a=-5; b=5;
  elseif tst==100   % etain
    T=2e-7;
    nx=200;
  elseif deux_plaques==1  % etain
    T=2e-7;
    a=-1e-3;   b=1e-3; 
    nx=400;
  elseif tst==200 || tst==201 || tst==202  % etain
    T=1e-7;
    nx=300;
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
  
  Nana=2e3;
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

Lcfl=[];

%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
##Lname=['Godunov acous o1, dP/dt Ttau=1e-9';
##       'Godunov acous o2, GAD MinMod, dP/dt Ttau=1e-9';
##       "vNR, q o1 Q+L, dP/dt Ttau=1e-9";
##       'BBC RK2av, q o1 Q+L, dP/dt Ttau=1e-9';
##       'RK1, q o1 Q+L, so 2, kin fix, dP/dt Ttau=1e-9';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix, dP/dt Ttau=1e-9';
##        ];
##
##Lsch=[002,011,300,202,100,105];  % scheme
##Lso=[-1,-1,-1,-1,2,5];           % spatial order
##Lz=[-1,-1,-1,-1,1,1];            % kinetic energy fix
##Lq=(004)*[-1,-1,1,1,1,1];        % artificial viscosity
##Ldpi=[-1,-1,-1,-1,0,0];          % dPI vs delta PI bar
##LTtau=1e-9*[1,1,1,1,1,1];
##Lcin_phase=3*[1,1,1,1,1,1 ];


##Lname=['Solveur acoustique';
##       'Ext. faiblement non-linéaire, equi thermo'; 
##       'Solveur Dukovicz, equi thermo';  
##       'Solveur double choc avec 5 itérations'];
##
##Lsch=[002,004,006,007];  % scheme
##Lso=[-1,-1,-1,-1];           % spatial order
##Lz=[-1,-1,-1,-1];            % kinetic energy fix
##Lq=(004)*[-1,-1,-1,-1];      % artificial viscosity
##Ldpi=[0,0,0,0];          % dPI vs delta PI bar
##LTtau=1e-9*[1,1,1,1];
##Lcin_phase=0*[0,0,0,0];



##Lname=['acoustic';
##       'Gallice 1';
##       'Gallice 2'; 
##       'Gallice 3'; 
##       ];
##
##Lsch=[2, 11, 12, 13];         % scheme
##Lso=[-1,-1,-1,-1];           % spatial order
##Lz=[-1,-1,-1,-1];            % kinetic energy fix
##Lq=(004)*[-1,-1,-1,-1];      % artificial viscosity
##Ldpi=[0,0,0,0];           % dPI vs delta PI bar
##LTtau=1e-9*[1,1,1,1];
##Lcin_phase=0*[0,0,0,0];

##Lname=['ext. o3 riemann';];
##
##Lsch=[5];         % scheme
##Lso=[-1];           % spatial order
##Lz=[-1];            % kinetic energy fix
##Lq=(004)*[-1];      % artificial viscosity
##Ldpi=[0];           % dPI vs delta PI bar
##LTtau=1e-9*[1];
##Lcin_phase=0*[0];


##Lname=['acoustic';
##       'BBC JCP 2009, q 4';
##       'BBC JCP 2009, q off';
##       'KOVA ordre 2';
##      ];
##
##Lsch=[2,200,200,401];         % scheme
##Lso=[-1,-1,-1,-1];           % spatial order
##Lz=[-1,-1,-1,-1];            % kinetic energy fix
##Lq=(004)*[-1,1,0,-1];      % artificial viscosity
##Ldpi=[0,0,0,0];           % dPI vs delta PI bar
##LTtau=1e-9*[1,1,1,1];
##Lcin_phase=0*[0,0,0,0];
##Lzeta=[-1,-1,-1,2.5];

Lname=['acoustic';
       'KOVA ordre 1 zeta=1';
       'KOVA ordre 1 zeta=2';
       'KOVA ordre 1 zeta=2.5';
       'KOVA ordre 1 zeta=3';
      ];

Lsch=[2,400,400,400,400];   
Lso=[-1,-1,-1,-1,-1];           % spatial order
Lz=[-1,-1,-1,-1,-1];            % kinetic energy fix
Lq=(004)*[-1,1,0,-1,-1];      % artificial viscosity
Ldpi=[0,0,0,0,0];           % dPI vs delta PI bar
LTtau=1e-9*[1,1,1,1,1];
Lcin_phase=0*[0,0,0,0,0];
Lzeta=[-1,1,2,2.5,3];


ls=length(Lsch);

% Permet de creer une liste Lzeta si il n'en n'existe pas (car le schema KOVA n'est pas lancée)
res=1;
for sch=Lsch
  res*=(1-(int8(sch/100)==4));
endfor
if res==1
  Lzeta=zeros(1,ls);
endif


tabGRILLE=(floor(Lsch/100)>=1);

% CFL
%numCFL=0.01;
%numCFL=0.6; %0.55
%numCFL=0.45;   % RK1 + sol acous
%numCFL=0.99;
%numCFL=0.9;
numCFL=0.7;

Ntau=15;

fac_dt=1e-1;

epsilon_newton=1e-12;

% Coeff de pseudo
Cq=1.5;
%Cq=0.75;
%Cq=0;
Cl=0.15;


% affichage
aff=0;
PLOT_NRJ=0;
LINE_WIDTH=0.5;
           PLOT=1   ;
PLOT_ITER=0;
PLOT_ERR=0;

sch_cin_phase=0;

% Multi fluide
tst_multi=-1;

% Max iter
   Nmax= 2005;


% symbole pour la legende
Symleg=['.-'; '-.'; '+-';'--', '*-'; 'x-'; 'o-'; 'd-'; 's-';'h-';'v-';':' ];

Etot=zeros(ls, 1e5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ls
  if length(Lcfl)>=1
    com=["./a.out  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(Lcfl(j))," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)," -j ",num2str(epsilon_newton)," -b ",num2str(Lzeta(j))];
  elseif length(Lcfl)==0
    com=["./a.out  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)," -j ",num2str(epsilon_newton)," -b ",num2str(Lzeta(j))];
  endif
  disp([' ']);
  disp([Lname(j,:)]);
  disp([' ']);
  disp([com]);
  system(com)

  resu=load("sol.txt");
  [~,n]=size(resu)
  %size(resu)
  Xd(j,:) = resu(1,:);  % X decentree (frontieres)
  Xc(j,:) = resu(2,1:n-1);
  rho(j,:) = resu(3,1:n-1);
  tau(j,:) = resu(4,1:n-1);
  u(j,:) = resu(5,:);
  e(j,:) = resu(6,1:n-1);
  epsilon(j,:) = resu(7,1:n-1);
  p(j,:) = resu(8,1:n-1);
  tempe(j,:) = resu(9,1:n-1);
  c(j,:) = resu(10,1:n-1);
  G(j,:) = resu(11,1:n-1);
  Sentropie(j,:) = resu(12,1:n-1);
  
  % energie du syteme
  NRJ=load("NRJ.txt");
  Nbiter(j)=length(NRJ(1,:));
  Etot(j,1:Nbiter(j))=NRJ(1,:);
  IMPUL(j,1:Nbiter(j))=NRJ(2,:);
  Mtot(j,1:Nbiter(j))=NRJ(3,:);
  VOL(j,1:Nbiter(j))=NRJ(4,:);
  
  resu=load("frac_mass.txt");
  [~,n]=size(resu);
  fmA(j,:)=resu(1,1:n-1);
  fmB(j,:)=resu(2,1:n-1);
  fmC(j,:)=resu(3,1:n-1);
  
  resu=load("nb_iter.txt");
  tabNiter(j,:)=resu(1,1:n-1);
  tabCRITERE(j,:)=resu(2,1:n-1);
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


 %% 
if PLOT==1 
  
if (tst>=113) 

  figure("name","Diag VE")
  hold on
  for j=1:ls
    plot(tau(j,:)*1e6,epsilon(j,:),Symleg(j,:));
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
  title(['Energie interne sp�cinfique,  s, nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([100 180]);
  ylim([-2e3 3e5]);
  %legend(Lname,'location','northeast');
  hold off
endif

if  tst==2 || (tst>=4 && tst<=100) || (tst>=113) 

  figure("name","E")
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique,  s, nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname,'location','northeast'); 
  hold off

  figure("name","V")
  hold on
  for j=1:ls
    plot(Xc(j,:),tau(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\tau nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("tau  (m^3/Kg)");
  %legend(Lname,'location','southwest');  
  hold off
  
  figure("name","rho")
  hold on
  for j=1:ls
    plot(Xc(j,:),rho(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname,'location','southwest'); 
  hold off
  
  figure("name","c")
  hold on
  for j=1:ls
    plot(Xc(j,:),c(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['c nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("c  (m/s)");
  %legend(Lname,'location','southwest'); 
  hold off
  
  figure("name","T")
  hold on
  for j=1:ls
    plot(Xc(j,:),tempe(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Temperature nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("T (K)");
  %legend(Lname,'location','southwest'); 
  hold off
  

  figure("name","u")
  hold on
  for j=1:ls
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,:),u(j,1:nx),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille d�cal�e
      plot(Xd(j,:),u(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
    end
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("u  (m/s)");
  legend(Lname,'location','southwest');
  hold off

  figure("name","p")
  hold on
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Pression nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("P  (Pa)");
  %legend(Lname,'location','southwest');
  hold off
endif


if  (tst==100 || tst>=113) 
  % fraction massique
  figure("name","beta")
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \beta nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  %legend(Lname,'location','southwest');
  hold off
  
  figure("name","gamma")
  hold on
  for j=1:ls
    plot(Xc(j,:),fmB(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \gamma nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  %legend(Lname,'location','southwest');
  hold off
  
  figure("name","liq")
  hold on
  for j=1:ls
    plot(Xc(j,:),fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase liquide nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  %legend(Lname,'location','southwest');
  hold off
  
  figure("name","phases tot")
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:)+fmB(j,:)+fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques nx=',num2str(nx)]);
  xlabel('x (m)')
  legend(Lname,'location','southwest');
  hold off
endif


 %%  solution exacte
if  (tst>=101 && tst<=112)
  
  [Uana,Pana,RHOana,Eana,Sana,Cana,Tana,Gana,lambda_ana,Xana, Xchoc1, Xchoc2]=sol_ana_etain(a,b,T,tst);
  
  %[Xoci,Uoci,Poci,Voci,Eoci,Coci,Soci,Dfoci,LAMBDAoci]=sol_oci_etain(T, u_init);
  
  %Lname=['solution onde composite';Lname];
  Lname_ana=['solution exacte';Lname];
  
  figure("name","diag VE")
  hold on
  plot((1./RHOana)*1e6,Eana,'r*',15);
  %plot((Voci)*1e6,Eoci);
  for j=1:ls
    plot(tau(j,:)*1e6,epsilon(j,:),Symleg(j,:),1.5);
  endfor
  plot(Vhug*1e6,Ehug,"linewidth",0.5);
  % Triangle
  plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r',0.5);
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
  xlim([110 140]);
  ylim([-1e2 1.4e5]);
  %legend(Lname_ana,'location','northeast');
  hold off
  
  %%% plot uniquement les points entre les deux chocs a droite
  
  tau_P1_plot=zeros(ls,n);
  epsilon_P1_plot=zeros(ls,n);
  for j=1:ls
    k=1;
    
    test=(Xc(j,:)<=Xchoc1).*(Xc(j,:)>=Xchoc2);
    tau_P1=test.*tau(j,:);
    epsilon_P1=test.*epsilon(j,:);
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
  plot((1./RHOana)*1e6,Eana,'r*',15);
  for j=1:ls
    plot(tau_P1_plot(j,1:Lk(j))*1e6,epsilon_P1_plot(j,1:Lk(j)),Symleg(j,:),1.5);
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
  title(['Energie interne spécifique, t=',num2str(T),' s, CFL= ',num2str(numCFL)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  xlim([120 130]);
  ylim([-1e2 1e5]);
  legend(Lname,'location','northeast');
  hold off
  
  figure("name","c")
  hold on
  plot(Xana,Cana,'-');
  for j=1:ls
    plot(Xc(j,:),c(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['c nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("c  (m/s)");
  %legend(Lname_ana,'location','southwest'); 
  hold off
  
  figure("name","G")
  hold on
  plot([-1e-3 1e-3], [0 0]);
  plot(Xana,Gana,'-');
  for j=1:ls
    plot(Xc(j,:),G(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['G nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("G  (m/s)");
  legend(Lname_ana,'location','southwest'); 
  hold off
  
  figure("name","S")
  hold on
  plot(Xana,Sana,'-');
  for j=1:ls
    plot(Xc(j,:),Sentropie(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['S nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("S  (m/s)");
  %legend(Lname_ana,'location','southwest'); 
  hold off
  
  figure("name","T")
  hold on
  plot(Xana,Tana,'-');
  for j=1:ls
    plot(Xc(j,:),tempe(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['temperature nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("T (K)");
  %legend(Lname_ana,'location','southwest'); 
  hold off
  
  figure("name","E")
  hold on
  plot(Xana,Eana,'-');
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp�cifique, t=',num2str(T),' s, nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("E  (J/Kg)");
  %legend(Lname_ana,'location','northeast'); 
  hold off

  figure("name","V")
  hold on
  plot(Xana,1./RHOana,'-');
  for j=1:ls
    plot(Xc(j,:),tau(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\tau nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("tau  (m^3/Kg)");
  %legend(Lname_ana,'location','southwest');  
  hold off
  
  figure("name","rho")
  hold on
  plot(Xana,RHOana,'-');
  for j=1:ls
    plot(Xc(j,:),rho(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("rho  (Kg/m^3)");
  %legend(Lname_ana,'location','southwest'); 
  hold off

  figure("name","u")
  hold on
  plot(Xana,Uana,'-');
  for j=1:ls
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,:),u(j,1:nx),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille d�cal�e
      plot(Xd(j,:),u(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
    end
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("u  (m/s)");
  %legend(Lname_ana,'location','southwest');
  hold off

  figure("name","P")
  hold on
  plot(Xana,Pana,'-');
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Pression nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("P  (Pa)");
  %legend(Lname_ana,'location','southwest');
  hold off

  % fraction massique
  figure("name","beta")
  hold on
  plot(Xana,lambda_ana(:,1),'-');
  for j=1:ls
    plot(Xc(j,:),fmA(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \beta nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  ylabel("fraction massique");
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname_ana,'location','southwest');
  hold off
  
  figure("name","gamma")
  hold on
  plot(Xana,lambda_ana(:,2),'-');
  for j=1:ls
    plot(Xc(j,:),fmB(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase \gamma nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')  
  ylabel("fraction massique");
  %legend(Lname_ana,'location','southwest');
  hold off
  
  figure("name","liq")
  hold on
  plot(Xana,lambda_ana(:,3),'-');
  for j=1:ls
    plot(Xc(j,:),fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase liquide nx=',num2str(nx),' , CFL= ',num2str(numCFL)]);
  xlabel('x (m)')
  ylabel("fraction massique");
  %legend(Lname_ana,'location','southwest');
  hold off
  
  figure("name","phases tot")
  hold on
  plot(Xana,lambda_ana(:,1)+lambda_ana(:,2)+lambda_ana(:,3),'-');
  for j=1:ls
    plot(Xc(j,:),fmA(j,:)+fmB(j,:)+fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques nx=',num2str(nx)]);
  xlabel('x (m)')
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

  figure("name","E")
  hold on
  plot(Xcana,Pana./((gamma-1)*RHOana),'-');
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp�cifique,  t=',num2str(T),' s , nx=',num2str(nx)]);
  xlabel("x (m)");
  ylabel(["E  (J/kg)"]);
  %legend(Lname_ana,'location','northwest');
  hold off

  figure("name","rho")
  hold on
  plot(Xcana,RHOana,'-','linewidth',LINE_WIDTH);
  for j=1:ls
    plot(Xc(j,:),rho(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx)]);
  xlabel("x (m)");
  ylabel(["tau (m^3/kg)"]);
  %legend(Lname_ana,'location','northwest');
  hold off

  figure("name","u")
  hold on
  if(tst==0 || tst==1)
    plot(Xdana,Uana,'-');
  elseif(tst==3)
    plot(Xcana,Uana,'-');
  endif
  for j=1:ls
    if(tabGRILLE(j)==0) % grille centr�e
      plot(Xc(j,:),u(j,1:nx),Symleg(j,:),'linewidth',LINE_WIDTH);
    elseif(tabGRILLE(j)==1) % grille d�cal�e
      plot(Xd(j,:),u(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
    end
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx)]);
  xlabel("x (m)");
  ylabel(["u (m/s)"]);
  legend(Lname_ana,'location','northwest');
  hold off

  figure("name","P")
  hold on
  plot(Xcana,Pana,'-');
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:),'linewidth',LINE_WIDTH);
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['p nx=',num2str(nx)]);
  xlabel("x (m)");
  ylabel(["P (Pa)"]);
  %legend(Lname_ana,'location','northwest');
  hold off

endif

  
  
   
  

endif % (PLOT=1)
  
if PLOT_NRJ==1
  % Plot de l'energie du syst�me au cours du temps
  figure("name","E tot")
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Etot(j,1:Nbiter(j)),Symleg(j,:));
  end
  grid();
  set(gca, "fontsize", 40);
  title(['energie totale']);
  legend(Lname,'location','southwest');
  hold off
  
  
  % Plot de l'impulsion au cours du temps
  figure("name","impulsion")
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,IMPUL(j,1:Nbiter(j)),Symleg(j,:));
  end
  grid();
  set(gca, "fontsize", 40);
  title(['impulsion']);
  legend(Lname,'location','southwest');
  hold off
  
  % Plot de la masse totale au cours du temps
  figure("name","masse tot")
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Mtot(j,1:Nbiter(j)),Symleg(j,:));
  end
  grid();
  set(gca, "fontsize", 40);
  title(['masse totale']);
  legend(Lname,'location','southwest');
  hold off
  
  % Plot du volume total au cours du temps
  figure("name","vol tot")
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,VOL(j,1:Nbiter(j)),Symleg(j,:));
  end
  grid();
  set(gca, "fontsize", 40);
  title(['volume total']);
  legend(Lname,'location','southwest');
  hold off
  
  for j=1:ls
    TV_Etot(j)=0;
    TV_IMPUL(j)=0;
    TV_Mtot(j)=0;
    TV_VOL(j)=0;
    for i=1:Nbiter(j)-1
      TV_Etot(j)=TV_Etot(j)+ abs(Etot(j,i+1) - Etot(j,i));
      TV_IMPUL(j)=TV_IMPUL(j)+ abs(IMPUL(j,i+1) - IMPUL(j,i));
      TV_Mtot(j)=TV_Mtot(j)+ abs(Mtot(j,i+1) - Mtot(j,i));
      TV_VOL(j)=TV_VOL(j)+ abs(VOL(j,i+1) - VOL(j,i));
    end
  end
  
  figure("name","conservation")
  hold on
  
  subplot(2,2,1);
  bar(TV_Etot);
  title("TV de energie totale sur le temps");
  xlabel('schemas');
  set(gca, "fontsize", 40);
  grid();
  
  subplot(2,2,2);
  bar(TV_IMPUL);
  title("TV de l'implusion sur le temps");
  xlabel('schemas');
  set(gca, "fontsize", 40);
  grid();
  
  subplot(2,2,3);
  bar(TV_Mtot);
  title("TV de la masse totale sur le temps");
  xlabel('schemas');
  set(gca, "fontsize", 40);
  grid();
  
  subplot(2,2,4);
  bar(TV_VOL);
  title("TV du volume totale sur le temps");
  xlabel('schemas');
  set(gca, "fontsize", 40);
  grid();
  
  hold off

endif


if PLOT_ITER 

  for j=1:ls
    figure("name",["Nb iteration ",Lname(j,:)])
    hold on
    bar(tabNiter(j,:));
    set(gca, "fontsize", 40);
    ylabel("Nombre d'itération par maille au dernier cycle en temps");
    xlabel('Numéro de la maille');
    title(["Nb iter | ",Lname(j,:)]);
    grid()
    hold off
  endfor

  for j=1:ls
    figure("name",["Critere ",Lname(j,:)])
    hold on
    plot(tabCRITERE(j,:),Symleg(j,:));
    set(gca, "fontsize", 40);
    ylabel('Critere max sur les cycles en temps par maille');
    xlabel('Numéro de la maille');
    title(["Critere max sur les cycles par maille | ",Lname(j,:)]);
    grid()
    hold off
  endfor
  
endif






if tst>100
for j=1:ls  Ln(j)=n;   end

%%   % Calcul de l'erreur sur P sur l'étain
if (tst>=101 && tst<=112) 
  cas=2;
  [ErrL1,ErrL2,ErrL1_red,ErrL2_red, ErrP1, ErrPf]=err_sol_etain(cas,T,tst,Ln,Xd,Xc,p);
endif

minL1=min(ErrL1);
minL2=min(ErrL2);
for j=1:ls
  ErrL1_pc(j)=((ErrL1(j)/minL1)-1)*100;
  ErrL2_pc(j)=((ErrL2(j)/minL2)-1)*100;
endfor

disp('');
disp('Erreur sur la pression ');
for j=1:ls
  disp(Lname(j,:))
  disp(['    L1 : ',num2str(ErrL1(j),' %.6g'),'  (',num2str(ErrL1_pc(j),' %.6g'),' +% du min)']);
  disp(['    L2 : ',num2str(ErrL2(j),' %.6g'),'  (',num2str(ErrL2_pc(j),' %.6g'),' +% du min)'])
endfor 

disp('');

if tst==102
  
  ErrP1_abs=ErrP1(:,1);
  ErrP1_rel=ErrP1(:,2);
  
  ErrPf_abs=ErrPf(:,1);
  ErrPf_rel=ErrPf(:,2);
  
  minP1_abs=min(ErrP1_abs);
  minPf=min(ErrPf_abs);
  for j=1:ls
    ErrL1_pc(j)=((ErrL1(j)/minL1)-1)*100;
    ErrL2_pc(j)=((ErrL2(j)/minL2)-1)*100;
  endfor

  disp('');
  disp('Erreur sur la pression de l etat P et l etat final');
  for j=1:ls
    disp(Lname(j,:))
    disp(['    P1 : ',num2str(ErrP1_abs(j),' %.6g'),'  (',num2str(ErrP1_rel(j),' %g'),' relatif)']);
    disp(['    Ef : ',num2str(ErrPf_abs(j),' %.6g'),'  (',num2str(ErrPf_rel(j),' %g'),' relatif)'])
  endfor 
endif

if PLOT_ERR
  
figure
hold on
bar(ErrL1)
ylabel('Erreur L1');
title(['Erreur L1 sur la pression ']);
legend(Lname);
set(gca, "fontsize", 40);
grid()
hold off

figure
hold on
bar(ErrL2)
ylabel('Erreur L2');
title(['Erreur L2 sur la pression ']);
legend(Lname);
set(gca, "fontsize", 40);
grid()
hold off

  
endif

endif


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%