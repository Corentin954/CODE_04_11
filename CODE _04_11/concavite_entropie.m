clear


% Condition de la concavitede l'entropie
cond3=load("fichiers/cond3_concavite.txt");
cond1_beta=load("fichiers/cond1_beta_concavite.txt");
cond1_gamma=load("fichiers/cond1_gamma_concavite.txt");
cond1_liq=load("fichiers/cond1_liq_concavite.txt");

cond3_test=load("fichiers/cond3_test_concavite.txt");
def_entropie=load("fichiers/def_entropie.txt");

V=load("fichiers/cond_V.txt");
E=load("fichiers/cond_E.txt");

[nv,ne]=size(cond1_beta);
cond1_beta=transpose(reshape(cond1_beta,nv,ne));
cond1_gamma=transpose(reshape(cond1_gamma,nv,ne));
cond1_liq=transpose(reshape(cond1_liq,nv,ne));


% Diagramme V,E
SA=load("fichiers/pSA.txt");
SB=load("fichiers/pSB.txt");
SC=load("fichiers/pSC.txt");
SAB=load("fichiers/pSAB.txt");
SAC=load("fichiers/pSAC.txt");
SBC=load("fichiers/pSBC.txt");

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             cond3  
figure
hold on

plot(V*1e6,0*V,'linewidth',3);
plot(V*1e6,cond3(:,1),'linewidth',3);
plot(V*1e6,cond3(:,2),'linewidth',3);
plot(V*1e6,cond3(:,3),'linewidth',3);
legend('x=0','\beta','\gamma','liquide')
set(gca, "fontsize", 40)
title("condition 3");
xlabel("V (cm^3/kg)");
grid();
hold off



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%             cond3 test 
##figure
##hold on
##
##plot(V*1e6,0*V,'linewidth',3);
##plot(V*1e6,cond3(:,1),'linewidth',3);
##plot(V*1e6,cond3_test(:,1),'linewidth',3);
##legend('x=0','\beta','\beta test')
##set(gca, "fontsize", 40)
##title("condition 3");
##xlabel("V (cm^3/kg)");
##grid();
##hold off
##
##
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%             Ensemble de définition de l'entropie
##figure
##hold on
##
##plot(V*1e6,def_entropie(:,1),'linewidth',3);
##plot(V*1e6,def_entropie(:,2),'linewidth',3);
##plot(V*1e6,def_entropie(:,3),'linewidth',3);
##legend('\beta','triangle','\gamma','liquide')
##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5);
##% Ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'r');
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'g');
##set(gca, "fontsize", 40)
##title("Ensemble de def de l'entropie \{(V,E) | E>f(V)\}");
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid();
##hold off
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%     Ensemble de définition de l'entropie BETA
##figure
##hold on
##plot(V*1e6,def_entropie(:,1),'--','linewidth',3,'color','r');
##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5,'color','g');
##% Ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##legend('limite de def','triangle','domaine')
##set(gca, "fontsize", 40)
##title("beta : Ensemble de def de l'entropie \{(V,E) | E>f(V)\}");
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid();
##hold off
##
##
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%     Ensemble de définition de l'entropie GAMMA
##figure
##hold on
##plot(V*1e6,def_entropie(:,2),'--','linewidth',3,'color','r');
##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5,'color','g');
##% Ploynomes
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##legend('limite de def','triangle','domaine')
##set(gca, "fontsize", 40)
##title("\gamma : Ensemble de def de l'entropie \{(V,E) | E>f(V)\}");
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid();
##hold off
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%     Ensemble de définition de l'entropie LIQUIDE
##figure
##hold on
##plot(V*1e6,def_entropie(:,3),'--','linewidth',3,'color','r');
##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5,'color','g');
##% Ploynomes
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##legend('limite de def','triangle','domaine')
##set(gca, "fontsize", 40)
##title("liquide  : Ensemble de def de l'entropie \{(V,E) | E>f(V)\}");
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid();
##hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             condition 1   BETA
figure
hold on

%surf(V*1e6,E,cond1_beta); %*1e-20);
contourf(V*1e6,E,cond1_beta,1)
text
##
##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5);
##% Ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'r');
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'g');

xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("Condition 1 BETA");
%title("Condition 1 (*1e-20) BETA");

%xlim([min(V) max(V)]);
%ylim([min(E) max(E)]);

hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             condition 1   BETA
figure
hold on

%surf(V*1e6,E,cond1_gamma); %*1e-20);
contourf(V*1e6,E,cond1_gamma,1);

##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5);
##% Ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'r');
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'g');

xlabel("V (cm^3/kg)");
ylabel("E ");
grid()
set(gca, "fontsize", 40)
title("Condition 1 GAMMA");
%title("Condition 1 (*1e-20) GAMMA");

##xlim([0.9*min(V) 1.1*max(V)]);
##ylim([0.9*min(E) 1.1*max(E)]);

hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             condition 1   LIQUIDE
figure
hold on

%surf(V*1e6,E,cond1_liq); %*1e-20);
contourf(V*1e6,E,cond1_liq,1);

##% Triangle
##plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'c', "linewidth",1.5);
##% Ploynomes
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'b');
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'r');
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",1.5,"color",'g');

xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
%title("Condition 1 (*1e-20) LIQUIDE");
title("Condition 1 LIQUIDE");

##xlim([min(V) max(V)]);
##ylim([min(E) max(E)]);

hold off