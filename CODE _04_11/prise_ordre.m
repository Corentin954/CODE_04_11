clear all

format long % permet d'avoir une dizaine de chiffre significatif


%//  Prise d'ordre  sur le cas-test limite acoustique  ///////////////////

tst=3;
a=0; b=1;
fac_eps=1e-8;

%%%%%%%%%%%%%%%%%%%%   
Lname=['sol acous o1';
       'GAD 2'; 
       'GoHy 2 so 2';
       'GoHy 3 so 2';
       'GoHy 3 so 3';
       'GoHy 3 so 5'];

Lsch=[02,10,13,20,20,20];  % scheme
Lso=[-1,-1,2,2,3,5];           % spatial order
Lz=[-1,-1,-1,-1,-1,-1];            % kinetic energy fix
Lq=(004)*[-1,-1,-1,-1,-1,-1];        % artificial viscosity  
Ldpi=[0,0,0,0,0,0];          % dPI vs delta PI bar
LTtau=[1,1,1,1,1,1];
Lcin_phase=0*[1,1,1,1,1,1];

##%%%%%%%%%%%%%%%%%%%%   
####Lname=['GoHy 3 so 3'];
####
####Lsch=[20];  % scheme
####Lso=[3];           % spatial order
####Lz=[-1];            % kinetic energy fix
####Lq=(004)*[-1];        % artificial viscosity  
####Ldpi=[0];          % dPI vs delta PI bar
####LTtau=[1];
####Lcin_phase=0*[1];


%%%%%%%%%%%%%%%%%%%%   
##Lname=['sol acous o1'; 
##       'GAD ordre 2'; 
##       'RK1, so 2';
##       'RK1, so 3';
##       'RK2, so 2';
##       'RK2, so 3'; 
##       'RK3, so 3'; 
##       'RK5, so 2'; 
##       'RK5 Cash-Karp, so 5'; 
##        ];
##
##Lsch=[02,10,100,100,101,101,102,104,104];  % scheme
##Lso=[-1,-1,2,3,2,3,3,2,5];           % spatial order
##Lz=[-1,-1,1,1,1,1,1,1,1];            % kinetic energy fix
##Lq=(004)*[-1,-1,1,1,1,1,1,1,1];        % artificial viscosity 
##Ldpi=[0,0,0,0,0,0,0,0,0];          % dPI vs delta PI bar
##LTtau=[1,1,1,1,1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1,1,1,1,1];


%%%%%%%%%%%%%%%%%%%%   
##Lname=['Sol acous o1';  
##       'Forward Euler, so 2';  
##       'GAD';  
##       'GoHy 2 so 2';  
##       'GoHy 3 so 2';  
##       'GoHy 3 so 3';  
##       'RK3, so 2';
##       'RK3, so 3';
##       'RK5, so 5';
##       ];
##
##Lsch=[002,113,10,13,20,20,102,102,104];  % scheme
##Lso=[-1,2,-1,2,2,3,2,3,5];           % spatial order
##Lz=[-1,1,-1,-1,-1,-1,1,1,1];            % kinetic energy fix
##Lq=(005)*[-1,1,-1,-1,-1,-1,1,1,1];        % artificial viscosity 
##Ldpi=[0,0,0,0,0,0,0,0,0];          % dPI vs delta PI bar
##LTtau=[1,1,1,1,1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1,1,1,1,1];

%%%%%%%%%%%%%%%%%%%%   
##Lname=['RK2 Heun, so 2'; 
##       'SPPRK3, so 3'; 
##       'RK5 Kash-Krap, so 5';  ];
##
##Lsch=[101,102,104];  % scheme
##Lso=[2,3,5];           % spatial order
##Lz=[1,1,1];            % kinetic energy fix
##Lq=(005)*[1,1,1];        % artificial viscosity 
##Ldpi=[0,0,0];          % dPI vs delta PI bar
##LTtau=[1,1,1];
##Lcin_phase=0*[1,1,1];


%%%%%%%%%%%%%%%%%%%%   
##Lname=['BBC JCP2009'; 
##       'BBC Pred-Corr'; 
##       'BBC RK2av';
##       'vNR';  ];
##
##Lsch=[200,201,202,300];  % scheme
##Lso=-1*[1,1,1,1];           % spatial order
##Lz=0*[1,1,1,1];            % kinetic energy fix
##Lq=(005)*[1,1,1,1];        % artificial viscosity 
##Ldpi=[0,0,0,0];          % dPI vs delta PI bar
##LTtau=[1,1,1,1];
##Lcin_phase=0*[1,1,1,1];


%%%%%%%%%%%%%%%%%%%%%%%%%
##Lname=['Despr�s'; 
##       'Jaouen'; 
##       'Sol acous o1'; ];
##
##Lsch=[000,001,002];  % scheme
##Lso=-1*[1,1,1];           % spatial order
##Lz=0*[1,1,1];            % kinetic energy fix
##Lq=(005)*[1,1,1];        % artificial viscosity 
##Ldpi=[0,0,0];          % dPI vs delta PI bar
##LTtau=[1,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%
Lname=['GoHy 3 so 5';  ];

Lsch=[20];  % scheme
Lso=[5];           % spatial order
Lz=[-1];            % kinetic energy fix
Lq=(005)*[1];        % artificial viscosity
Ldpi=[0];          % dPI vs delta PI bar
LTtau=[1];
Lcin_phase=0*[1];

ls=length(Lsch);
tabGRILLE=(floor(Lsch/100)>=1);


Ntau=15;

fac_dt=1;

% Coeff de pseudo
Cq=1.5;
Cq=0;
Cl=0.15;
Cl=0;

% affichage
aff=0;
PLOT_NRJ=0;
LINE_WIDTH=0.5;


sch_cin_phase=0;

% Multi fluide
tst_multi=-1;

% Max iter
   Nmax=2e6;


% nx 
     %Lnx=[ 40,80, 160,320, 640, 1280, 2560]; 
     %Lnx=[ 10, 20 ,40,80, 160]; 
     Lnx=[40, 80, 160]; 

lnx=length(Lnx);

% CFL
Lcfl=0.1*ones(1,lnx);
 

Lcfl(1)=0.4;
%facCFL=(1/7)^(1/3); % pour GoHy 3 so 5 --> ordre 6 ??
facCFL=(1/2)^(2/3); % pour GoHy 3 so 5 --> ordre 6 ??
%facCFL=(1/2)^(4/3); % pour GoHy 3 so 7 -->  
%facCFL=(1/2); % pour RK1 so 2 -->  
##facCFL=(1/3); % pour RK1so 3 -->  
for i=2:lnx
  Lcfl(i)=Lcfl(i-1)*facCFL;
endfor


ErrTAU=zeros(ls, lnx);
ErrU=zeros(ls, lnx);

% symoble pour la legende
Symleg=['x-';'o-'; '+-';'*-'; 'd-'; 's-';'h-';'v-';':' ];


% =======================================================================
% Creation de l'executable
system("sh compil_main_multi.sh");


j=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sch=Lsch
  k=1;
  for nx=Lnx
    nx
    disp([' ']);
    Lname(j,:)
    %["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)]
    %system(["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)])
    system(["./a.out  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(Lcfl(k))," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)])
    
    resu=load("sol.txt");
    [~,n]=size(resu);
    %size(resu);
    Xd = resu(1,:); % X decentree (frontieres)
    Xc = resu(2,1:n-1);
    rho = resu(3,1:n-1);
    tau = resu(4,1:n-1);
    if (tabGRILLE(j)==0)
      u = resu(5,1:n-1);
    elseif (tabGRILLE(j)==1)
      u = resu(5,:);
    end

    % Calcul de la solution analytique sur la grille initiale
    [TAUana, RHOana, Uana, RTAUana, RUana, REPSana]=sol_ordre(a,b,nx,fac_eps,tabGRILLE(j));
    
    dx=(b-a)/nx;
    % Calcul de l'erreur
    ErrRHO(j,k)= norm(RHOana - rho,1)*dx;
    ErrTAU(j,k)= norm(TAUana - tau,1)*dx;
    ErrU(j,k)= norm(Uana-u,1)*dx;
    
    Xd_ana=a:dx:b;
    Xc_ana=(a+dx/2):dx:(b-dx/2);
    ErrXd(j,k)=norm(Xd-Xd_ana,1)*dx;
    ErrXc(j,k)=norm(Xc-Xc_ana,1)*dx;
    
    k=k+1;
  endfor
  j=j+1;
endfor


% Courbes d'erreur sur les varaibles Euleriennes  (log- log)  %%%%%%%%%%%%%%%%

logN=log10(Lnx);


for j=1:ls
  p=polyfit(logN,log10(ErrTAU(j,:)),1);
  ordreTAU(j)=-p(1);
  p=polyfit(logN,log10(ErrRHO(j,:)),1);
  ordreRHO(j)=-p(1);
  p=polyfit(logN,log10(ErrU(j,:)),1);
  ordreU(j)=-p(1);  
  
  p=polyfit(logN,log10(ErrXd(j,:)),1);
  ordreXd(j)=-p(1);
  p=polyfit(logN,log10(ErrXc(j,:)),1);
  ordreXc(j)=-p(1);
end


##j=1;
##list=1:6;
##p=polyfit(logN(list),log10(ErrTAU(j,list)),1);
##ordreTAU(j)=-p(1);
##p=polyfit(logN(list),log10(ErrU(j,list)),1);
##ordreU(j)=-p(1);
##
##j=2;
##list=1:ls;
##p=polyfit(logN(list),log10(ErrTAU(j,list)),1);
##ordreTAU(j)=-p(1);
##p=polyfit(logN(list),log10(ErrU(j,list)),1);
##ordreU(j)=-p(1);
##
##j=3;
##list=1:ls;
##p=polyfit(logN(list),log10(ErrTAU(j,list)),1);
##ordreTAU(j)=-p(1);
##p=polyfit(logN(list),log10(ErrU(j,list)),1);
##ordreU(j)=-p(1);
##
##j=4;
##list=3:5;
##p=polyfit(logN(list),log10(ErrTAU(j,list)),1);
##ordreTAU(j)=-p(1);
##p=polyfit(logN(list),log10(ErrU(j,list)),1);
##ordreU(j)=-p(1);
##
##j=5;
##list=3:5;
##p=polyfit(logN(list),log10(ErrTAU(j,list)),1);
##ordreTAU(j)=-p(1);
##p=polyfit(logN(list),log10(ErrU(j,list)),1);
##ordreU(j)=-p(1);



for j=1:ls
  Lname_TAU(j,:) =[Lname(j,:),' | ',num2str(abs(ordreTAU(j)),"%.2f")];
  Lname_RHO(j,:) =[Lname(j,:),' | ',num2str(abs(ordreRHO(j)),"%.2f")];
  Lname_U(j,:) =[Lname(j,:),' | ',num2str(abs(ordreU(j)),"%.2f")];

  Lname_Xd(j,:) =[Lname(j,:),' | ',num2str(abs(ordreXd(j)),"%.2f")];
  Lname_Xc(j,:) =[Lname(j,:),' | ',num2str(abs(ordreXc(j)),"%.2f")];
end

for j=1:ls 
  disp([Lname(j,:)])
  disp([' -TAU:', num2str(ordreTAU(j),2),' -RHO:', num2str(ordreRHO(j),2), '  -U :', num2str(ordreU(j),2)]);
  disp('')
end

% figure
figure("name","erreur tau");
hold on
for j=1:ls
  plot(logN,log10(ErrTAU(j,:)),Symleg(j,:),'linewidth',1);
endfor
xlabel('log10 nx');
ylabel('log10 Err L1 ');
set(gca, "fontsize", 40);
grid();
legend(Lname_TAU,'location','southwest');
title("Erreur L1 tau (log-log)");
hold off

figure("name","erreur rho");
hold on
for j=1:ls
  plot(logN,log10(ErrRHO(j,:)),Symleg(j,:),'linewidth',1);
endfor
xlabel('log10 nx');
ylabel('log10 Err L1 ');
set(gca, "fontsize", 40);
grid();
legend(Lname_RHO,'location','southwest');
title("Erreur L1 rho (log-log)");
hold off

figure("name","erreur u");
hold on
for j=1:ls
  plot(logN,log10(ErrU(j,:)),Symleg(j,:),'linewidth',1);
endfor
xlabel('log10 nx');
ylabel('log10 Err L1');
set(gca, "fontsize", 40);
grid();
legend(Lname_U,'location','southwest');
title("Erreur L1 u (log-log)");
hold off


figure("name","erreur Xd");
hold on
for j=1:ls
  plot(logN,log10(ErrXd(j,:)),Symleg(j,:),'linewidth',1);
endfor
xlabel('log10 nx');
ylabel('log10 Err L1');
set(gca, "fontsize", 40);
grid();
legend(Lname_Xd,'location','southwest');
title("Erreur L1 grill� d�cal�e (log-log)");
hold off


figure("name","erreur Xc");
hold on
for j=1:ls
  plot(logN,log10(ErrXc(j,:)),Symleg(j,:),'linewidth',1);
endfor
xlabel('log10 nx');
ylabel('log10 Err L1');
set(gca, "fontsize", 40);
grid();
legend(Lname_Xc,'location','southwest');
title("Erreur L1 grill� centr�e (log-log)");
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('comp_err_RKint_SK2_workspace.mat');