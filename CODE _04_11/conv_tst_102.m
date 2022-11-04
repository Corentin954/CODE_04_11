clear all


%///////////////             RUNGE-KUTTA               ///////////////////
%///////////////    Plot  schema + sol analytique        ///////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix de la loi d'etat et du cas test
%  tst : 0 (gaz parfait)   2
tst = 102; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=2e-7;

% =======================================================================
% Creation de l'executable
system("sh compil_main_multi.sh");


##Lname=['Godunov solveur acoustique, equi thermo';
##       'Ext. faiblement non-linéaire solveur acous, equi thermo'; 
##      ];
##
##Lsch=[002,004];      % scheme
##Lso=[-1,-1];           % spatial order
##Lz=[-1,-1];            % kinetic energy fix
##Lq=(004)*[-1,-1];      % artificial viscosity
##Ldpi=[-1,-1];          % dPI vs delta PI bar
##LTtau=1e-9*[1,1];
##Lcin_phase=0*[0,0];

##Lname=['Solveur Dokovicz, equi thermo';
##      ];
##
##Lsch=[005];      % scheme
##Lso=[-1];           % spatial order
##Lz=[-1];            % kinetic energy fix
##Lq=(004)*[-1];      % artificial viscosity
##Ldpi=[-1];          % dPI vs delta PI bar
##LTtau=1e-9*[1];
##Lcin_phase=0*[0];

Lname=['Solveur acoustique, equi thermo';
       'Ext. faiblement non-linéaire, equi thermo'; 
       'Solveur Dukovicz, equi thermo';  ];

Lsch=[002,004,005];  % scheme
Lso=[-1,-1,-1];           % spatial order
Lz=[-1,-1,-1];            % kinetic energy fix
Lq=(004)*[-1,-1,-1];      % artificial viscosity
Ldpi=[0,0,0];          % dPI vs delta PI bar
LTtau=1e-9*[1,1,1];
Lcin_phase=0*[0,0,0];

ls=length(Lsch);

tabGRILLE=(floor(Lsch/100)>=1);

###
Ln=[200, 400, 800, 1600, 3200];

% CFL
numCFL=0.5;
numCFL=0.90;
numCFL=0.4;
 
Ntau=15;

fac_dt=1e-3;

epsilon_newton=1e-13;

% Coeff de pseudo
Cq=1.5;
%Cq=0;
Cl=0.15;
%Cl=0;
%Cl=0.29;

% affichage
aff=0;

sch_cin_phase=0;

% Multi fluide
tst_multi=-1;

% Max iter
   Nmax= 2e6;

% symoble pour la legende
Symleg=['o-'; '*-'; '+-';'x-'; 'd-'; 's-';'h-';'v-';':';'^-';'v-' ];



maxN=max(Ln)+10;

Xd = zeros(ls,maxN); % X d�centr�e (fronti�res)
Xc = zeros(ls,maxN);
tau = zeros(ls,maxN);
u = zeros(ls,maxN);
e = zeros(ls,maxN);
p = zeros(ls,maxN);
epsilon = zeros(ls,maxN);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ls
    for kn=1:length(Ln)

      nx=Ln(kn);

      com=["./a.out  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)," -j ",num2str(epsilon_newton)];
      disp([' ']);
      disp([Lname(j,:)]);
      disp([' ']);
      disp([com]);
      system(com)
      
      resu=load("sol.txt");
      [~,n]=size(resu)

      Xd(kn,1:n) = resu(1,:);
      Xc(kn,1:n-1) = resu(2,1:n-1);
      rho(kn,1:n-1) = resu(3,1:n-1);
      tau(kn,1:n-1) = resu(4,1:n-1);
      if tabGRILLE(j)==0
        u(kn,1:n-1) = resu(5,1:n-1);
      elseif tabGRILLE(j)==1
        u(kn,1:n) = resu(5,:);
      end
      e(kn,1:n-1) = resu(6,1:n-1);  
      epsilon(kn,1:n-1) = resu(7,1:n-1);
      p(kn,1:n-1) = resu(8,1:n-1);
      tempe(kn,1:n-1) = resu(9,1:n-1);
      c(kn,1:n-1) = resu(10,1:n-1);
      G(kn,1:n-1) = resu(11,1:n-1);
      Sentropie(kn,1:n-1) = resu(12,1:n-1);
      
      resu=load("frac_mass.txt");
      [~,n]=size(resu);
      fmA(kn,1:n-1)=resu(1,1:n-1);
      fmB(kn,1:n-1)=resu(2,1:n-1);
      fmC(kn,1:n-1)=resu(3,1:n-1);
      
  endfor
  
  % Calcul de l'erreur 
  cas=2;
  [ErrL1,ErrL2,ErrL1_red,ErrL2_red, ErrP1, ErrPf]=err_sol_etain(cas,T,tst,Ln+1,Xd,Xc,p);

  ErrP1_abs(j,:) = ErrP1(1,:) ;
  ErrP1_rel(j,:) = ErrP1(2,:) ;

  ErrPf_abs(j,:) = ErrPf(1,:) ;
  ErrPf_rel(j,:) = ErrPf(2,:) ;


endfor


%%%%%%%%%%%%%%%%   figures  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% prise d'ordre
for j=1:ls
  po=polyfit(log10(Ln),log10(ErrP1_abs(j,:)),1);
  ordreP1(j)=-po(1)
  list=[Lname(j,:),'| P1 |',num2str(ordreP1(j),2)];
  LnameP1(j,1:length(list)) = list;
  Lname_P(2*j-1,1:length(list)) = list
  po=polyfit(log10(Ln),log10(ErrPf_abs(j,:)),1);
  ordrePf(j)=-po(1)
  list=[Lname(j,:),'| EF |',num2str(ordrePf(j),2)];
  LnamePf(j,1:length(list)) = list; 
  Lname_P(2*j,1:length(list)) = list
endfor

figure
hold on
for j=1:ls
  plot(log10(Ln),log10(ErrP1_abs(j,:)),Symleg(2*j-1,:));
  plot(log10(Ln),log10(ErrPf_abs(j,:)),Symleg(2*j,:));
endfor
xlabel('log nombre de mailles');
ylabel('log erreur');
title(['Erreur log-log absolue sur la pression de l etat P1 et etat final']);
legend(Lname_P);
set(gca, "fontsize", 40);
grid()
hold off

figure
hold on
for j=1:ls
  plot(log10(Ln),log10(ErrP1_rel(j,:)),Symleg(2*j-1,:));
  plot(log10(Ln),log10(ErrPf_rel(j,:)),Symleg(2*j,:));
endfor  
set(gca, "fontsize", 40);
xlabel('log nombre de mailles');
ylabel('log erreur');
title(['Erreur log-log relative sur la pression sur P1 et etat final']);
legend(Lname_P);
grid()
hold off

figure
hold on
for j=1:ls
  loglog(Ln,ErrP1_rel(j,:),Symleg(2*j-1,:));
  loglog(Ln,ErrPf_rel(j,:),Symleg(2*j,:));
endfor  
set(gca, "fontsize", 40);
xlabel('nombre de mailles');
ylabel('erreur');
title(['Erreur log-log relative sur la pression sur P1 et etat final']);
legend(Lname_P);
grid()
hold off


figure
hold on
for j=1:ls
  plot(log10(Ln),log10(ErrP1_rel(j,:)),Symleg(j,:));
endfor
xlabel('log nombre de mailles');
ylabel('log erreur L1');
title(['Erreur log-log relative sur la pression sur P1']);
legend(LnameP1);
set(gca, "fontsize", 40);
grid()
hold off


figure
hold on
for j=1:ls
  plot(log10(Ln),log10(ErrPf_abs(j,:)),Symleg(j,:));
endfor
xlabel('log nombre de mailles');
ylabel('log erreur L1');
title(['Erreur log-log absolue sur la pression sur etat final']);
legend(LnamePf);
set(gca, "fontsize", 40);
grid()
hold off

figure
hold on
for j=1:ls
  plot(log10(Ln),log10(ErrPf_rel(j,:)),Symleg(j,:));
endfor
xlabel('log nombre de mailles');
ylabel('log erreur L1');
title(['Erreur log-log relative sur la pression sur etat final']);
legend(LnamePf);
set(gca, "fontsize", 40);
grid()
hold off


%save('SAVE\Interval_tst_104_God.mat');

% Afichage des erreurs pour les copier dans un tableur^
for j=1:ls
  disp(Lname(j,:))    
  disp('nx ; err absolue ; err relative:'); 
  disp('Etat 1');
  for i=1:length(Ln)
    disp([num2str(Ln(i)),' ; ',num2str(ErrP1_abs(j,i),'%.5g'),' ; ',num2str(ErrP1_rel(j,i),'%.5g')]);
  endfor
  disp('Etat final');
  for i=1:length(Ln)
    disp([num2str(Ln(i)),' ; ',num2str(ErrPf_abs(j,i),'%.5g'),' ; ',num2str(ErrPf_rel(j,i),'%.5g')]);
  endfor
  disp('');
endfor

