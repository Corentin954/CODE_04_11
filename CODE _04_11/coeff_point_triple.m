

clear


temp=load("fichiers/PARAMS_triple.txt");
dPdT12=load("fichiers/dPdT12.txt");
dPdT13=load("fichiers/dPdT13.txt");

PT=temp(:,1);
TT=temp(:,2);
distance=temp(:,3);

[n12,~]=size(dPdT12);
[n13,~]=size(dPdT13);

PT=reshape(PT,n12,n13);
TT=reshape(TT,n12,n13);
distance=reshape(distance,n12,n13);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          distance au point triple en fct de dPT12 et dPdT13 
figure
hold on
%surf(V*1e6,E,P*1e-9);
%%contour(V*1e6,E,P,50);
%%pcolor(V*1e6,E,P);
%contourf(V*1e6,E,P*1e-9,[0:2:70]);
contourf(dPdT12*1e-7,dPdT13*1e-7,distance,[0:0.001:0.05]);
# point
x=-2.13; y=3.125;
plot(x,y,'marker','+','markersize',15);
text(x+0.005,y,'Heuzé','color','blue',"fontsize", 30);


xlabel('dP/dT \gamma (1e^7 Pa/K)');
ylabel('dP/dT liq (1e^7 Pa/K)');
grid()
set(gca, "fontsize", 40)
title('Distance |P-P_{ref}|/P_{ref} + |T-T_{ref}|/T_{ref}');
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          distance au point triple en fct de dPT12 et dPdT13 
figure
hold on
surf(dPdT12*1e-7,dPdT13*1e-7,distance);
%contourf(dPdT12*1e-7,dPdT13*1e-7,distance);

xlabel("dPdT12 (1e7 Pa/K)");
ylabel("dPdT13 (1e7 Pa/K)");
grid()
set(gca, "fontsize", 40)
title("Distance |P-Pref|/Pref + |T-Tref|/Tref");
hold off








