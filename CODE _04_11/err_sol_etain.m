
% Erreur sur la solution référence HUGONIOT etain

% CAS==1 
%    calcul de la pression de reference par intégration du creneau sur les mailles
% CAS==2
%    pression ponctuelle 
% !! différence sur le mailles comportant la discontinuité !!

function [ErrL1,ErrL2,ErrL1_red,ErrL2_red, ErrP1, ErrPf]=err_sol_etain(cas,T,tst,Ln,Xd,Xc,p)

xa_P1=-56e-5;
xb_P1=-47.5e-5;

xa_Pf=-20e-5;
xb_Pf=20e-5;

ls=length(Ln);

P0=0; P1=7.673987e9; P2=8.478062e9 ; P3=11.353850e9  ; P4=44.866639e9 ; P5=56.903493e9; 	
rho0=7287.00 ; rho1=8076.00 ; rho2=8288.37 ; rho3=8518.28 ; rho4=10247.86 ; rho5=10609.00;  
E0=-0.009371 ; E1=51443.0 ; E2=77064.92 ; E3=112608.341 ; E4=889466.87 ; E5=1222604.76; 
u0=0; u1=320.75; u2=371.265; u3=474.57; u4=1333.766; u5=1563.716;

lambda0=[1,0,0];
lambda1=[1,0,0];

 
lambda(101,:)=[1,0,0];
p_hug(101)= P1; u_hug(101)= u1; rho_hug(101)= rho1; E_hug(101)= E1;

% ZONE MIXTE beta/gamma (1/2)
% lambda=0.2
lambda(102,:)=[0.8,0.2,0];
p_hug(102)= 7.83711e+09; u_hug(102)= 331.046; rho_hug(102)= 8118.55; E_hug(102)= 56474.8;
% lambda=0.5
lambda(103,:)=[0.5,0.5,0];
p_hug(103)= 8.07994e+09; u_hug(103)= 346.31; rho_hug(103)= 8182.28; E_hug(103)= 64111.2;
% lambda=0.8
lambda(104,:)=[0.2,0.8,0];
p_hug(104)= 8.31994e+09; u_hug(104)= 361.359; rho_hug(104)= 8245.95; E_hug(104)= 71850.5;

lambda(105,:)=[0,1,0];
p_hug(105)= P2; u_hug(105)= u2; rho_hug(105)= rho2; E_hug(105)= E2;

lambda(106,:)=[0,1,0]; 
p_hug(106)= P3; u_hug(106)= u3; rho_hug(106)= rho3; E_hug(106)= E3;

lambda(107,:)=[0,1,0];
p_hug(107)= P4; u_hug(107)= u4; rho_hug(107)= rho4; E_hug(107)= E4;


% ZONE MIXTE gamma/liquide (2/3)
% lambda=0.2
lambda(108,:)=[0,0.8,0.2];
p_hug(108)= 4.73658e+10; u_hug(108)= 1383.13; rho_hug(108)= 10326.1; E_hug(108)= 956520;
% lambda=0.5
lambda(109,:)=[0,0.5,0.5];
p_hug(109)= 5.10185e+10; u_hug(109)= 1453.67; rho_hug(109)= 10437.2; E_hug(109)= 1.05658e+06;
% lambda=0.8
lambda(110,:)=[0,0.2,0.8];
p_hug(110)= 5.45762e+10; u_hug(110)= 1520.7; rho_hug(110)= 10542; E_hug(110)= 1.15626e+06;

lambda(111,:)=[0,0,1];
p_hug(111)= P5; u_hug(111)= u5; rho_hug(111)= rho5; E_hug(111)= E5;

lambda(112,:)=[0,1,0];
u_hug(112)=900; rho_hug(112)=1./0.000106217; E_hug(112)=405000; p_hug(112)=2.61173e+10;


% simple choc
if(tst==101 || tst==106 || tst==107 || tst==108 || tst==109 || tst==110 || tst==111 || tst==112)
  V0=1./rho0;
  V=1./rho_hug(tst);
  choc1=V0*sqrt((p_hug(tst)-P0)/(V0-V))*T - u_hug(tst)*T;
   
  f_pe=@(x) (p_hug(tst)-P0)*(abs(x)<=choc1) + P0;
    
  for j=1:ls
    n=Ln(j);
    ps=p(j,1:n-1);
    xc=Xc(j,1:n-1);
    
    if cas==1
      for i=1:n-1
        xa=Xd(j,i); xb=Xd(j,i+1);
        if(xb<=-choc1 || xa>=choc1)
          pe(i)=P0;
        elseif(xa>=-choc1 && xb<=choc1)
          pe(i)=p_hug(tst);
        elseif(xa<=-choc1 && xb>=-choc1)
          pe(i)=(p_hug(tst)*(xb+choc1) + P0*(-choc1-xa))/(xb-xa);
        elseif(xa<=choc1 && xb>=choc1)
          pe(i)=(P0*(xb-choc1) + p_hug(tst)*(choc1-xa))/(xb-xa);
        else
          disp("Erreur sur le test sur Xd");
        endif
      endfor
    elseif cas==2
      pe=f_pe(xc);  
    endif
    
    ErrL1(j)=norm(ps-pe,1)/(n-1);
    ErrL2(j)=norm(ps-pe,2)/(n-1);
  endfor
else 
% double choc
  V0=1./rho0;
  V1=1./rho1;
  V=1./rho_hug(tst);
  choc1=V0*sqrt((P1-P0)/(V0-V1))*T - u_hug(tst)*T;
  choc2=V1*sqrt((p_hug(tst)-P1)/(V1-V))*T  + u1*T - u_hug(tst)*T;
  
  f_pe=@(x) (p_hug(tst)-P1)*(abs(x)<choc2) + (P1-P0)*(abs(x)<choc1) + P0;
    
  for j=1:ls
    n=Ln(j);
    ps=p(j,1:n-1);
    xc=Xc(j,1:n-1);
    xd=Xd(j,1:n);
        
    if cas==1
      for i=1:n-1
        xa=Xd(j,i); xb=Xd(j,i+1);
        if(xb<=-choc1 || xa>=choc1)
          pe(i)=P0;
        elseif(xa>=-choc1 && xb<=-choc2)
          pe(i)=P1;
        elseif(xa>=choc2 && xb<=choc1)
          pe(i)=P1;
        elseif(xa>=-choc2 && xb<=choc2)
          pe(i)=p_hug(tst);
        elseif(xa<=-choc1 && xb>=-choc1)
          pe(i)=(P0*(-choc1-xa) + P1*(xb+choc1))/(xb-xa);
        elseif(xa<=-choc2 && xb>=-choc2)
          pe(i)=(P1*(-choc2-xa) + p_hug(tst)*(xb+choc2))/(xb-xa);
        elseif(xa<=choc2 && xb>=choc2)
          pe(i)=(p_hug(tst)*(choc2-xa) + P1*(xb-choc2))/(xb-xa);
        elseif(xa<=choc1 && xb>=choc1)
          pe(i)=(P1*(choc1-xa) + P0*(xb-choc1))/(xb-xa);
        else
          disp("Erreur sur le test sur Xd");
        endif
      endfor
    elseif cas==2
      pe=f_pe(xc);  
    endif
    
    err=ps-pe;
    err_plot=err;
 
 % On enleve les erreurs autours des discontinuitées
    l=length(err_plot);
    nb_err_1=7; %*((n-1)/1600) % le nombre d'erreur que l'on souhaite enlever de chaque coté du choc 1
    nb_err_tot_1=1+2*nb_err_1;

    nb_err_2=5*((n-1)/200); % le nombre d'erreur que l'on souhaite enlever de chaque coté du choc 2
    nb_err_tot_2=1+2*nb_err_2;
    for i=1:n-1
      xa=Xd(j,i); xb=Xd(j,i+1);
      test_1=( (xa<=-choc1 && xb>=-choc1) || (xa<=+choc1 && xb>=+choc1) );
      test_2=( (xa<=-choc2 && xb>=-choc2) || (xa<=+choc2 && xb>=+choc2) );
      if test_1
        err_plot(i-nb_err_1:i+nb_err_1)=0;
        l=l-nb_err_tot_1;
      elseif test_2
        err_plot(i-nb_err_2:i+nb_err_2)=0;
        l=l-nb_err_tot_2;
      endif
    endfor


    ErrL1_red(j)=norm(err_plot,1)/l;
    ErrL2_red(j)=norm(err_plot,2)/l;

    ErrL1(j)=norm(ps-pe(1:n-1),1)/(n-1);
    ErrL2(j)=norm(ps-pe(1:n-1),2)/(n-1);
    
    % Erreur sur des plages fixée pour l'état 1
    if tst==102

      % P1
      [fmax,ia]=min(abs(xd-xa_P1));
      [fmax,ib]=min(abs(xd-xb_P1));
      
      if(xd(ia)>xa_P1)  ia+=1;  endif
      if(xd(ib)<xb_P1)  ib-=1;  endif
      
      S_P1=(xd(ia)-xa_P1)*ps(ia-1);
      S_P1+=(xb_P1-xd(ib))*ps(ib+1);      
      for i=ia:ib-1
        S_P1+=ps(i)*(xd(i+1)-xd(i));
      endfor
      S_P1=S_P1/(xb_P1 - xa_P1);

##      S_ex=(xd(ia)-xa_P1)*(ps(ia-1)-P1);
##      S_ex+=(xb_P1-xd(ib))*(ps(ib+1)-P1);
##      for i=ia:ib-1
##        S_ex+=(ps(i)-P1)*(xd(i+1)-xd(i));
##      endfor
##      S_ex=S_ex/(xb_P1-xa_P1)
      
      ErrP1(1,j)=abs(P1 - S_P1);
      ErrP1(2,j)=abs(P1 - S_P1)/abs(P1);

      % Etat final
      [fmax,ia]=min(abs(xd-xa_Pf));
      [fmax,ib]=min(abs(xd-xb_Pf));
      
      if(xd(ia)>xa_Pf)  ia+=1;  endif
      if(xd(ib)<xb_Pf)  ib-=1;  endif

      S_Pf=(xd(ia)-xa_Pf)*ps(ia-1);
      S_Pf+=(xb_Pf-xd(ib))*ps(ib+1);
      for i=ia:ib-1
        S_Pf+=ps(i)*(xd(i+1)-xd(i));
      endfor
      S_Pf = S_Pf/(xb_Pf - xa_Pf);

      ErrPf(1,j)=abs(p_hug(tst) - S_Pf);
      ErrPf(2,j)=abs(p_hug(tst) - S_Pf)/abs(p_hug(tst));
    endif
    

    Pmax=max(ps);  Pmin=min(ps);

    ##    figure("name","P")
##    hold on
##    plot(xc,pe(1:n-1),'-+');
##    plot(xc,ps,'--');
##    plot(xc,ps-pe(1:n-1),'-*');
##    plot(xc,err_plot,'-o');
##    plot([xa_P1 xa_P1], [Pmin Pmax], 'g');
##    plot([xb_P1 xb_P1], [Pmin Pmax], 'g');
##    plot([xa_Pf xa_Pf], [Pmin Pmax], 'm');
##    plot([xb_Pf xb_Pf], [Pmin Pmax], 'm');
##    grid();
##    set(gca, "fontsize", 40);
##    title(['p simulation ',num2str(j)]);
##    xlabel("x (m)");
##    ylabel(["P (Pa)"]);
##    legend('ref','sim','err','err reduce');
##    hold off
##    
##    figure("name","P")
##    hold on
##    plot(xc,pe(1:n-1),'-');
##    plot(xc,ps,'.-');
##    plot([xa_P1 xa_P1], [Pmin Pmax], 'g','linewidth',2);
##    plot([xb_P1 xb_P1], [Pmin Pmax], 'g','linewidth',2);
##    plot([xa_Pf xa_Pf], [Pmin Pmax], 'm','linewidth',2);
##    plot([xb_Pf xb_Pf], [Pmin Pmax], 'm','linewidth',2);
##    grid();
##    set(gca, "fontsize", 40);
##    title(['p simulation ',num2str(j)]);
##    xlabel("x (m)");
##    ylabel(["P (Pa)"]);
##    legend('ref','sim');
##    hold off

  endfor
endif


endfunction