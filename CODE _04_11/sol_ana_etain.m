
% Solution exacte HUGONIOT etain

function [Uana,Pana,RHOana,Eana,Sana,Cana,Tana,Gana,LAMBDA,x,choc1,choc2]=sol_ana_etain(a,b,time,tst)

n_ana=1e4;

dx=1e-8;

P0=0; P1=7.673987e9; P2=8.478062e9 ; P3=11.353850e9  ; P4=44.866639e9 ; P5=56.903493e9; 	
rho0=7287.00 ; rho1=8076.00 ; rho2=8288.37 ; rho3=8518.28 ; rho4=10247.86 ; rho5=10609.00;  
E0=-0.009371 ; E1=51443.0 ; E2=77064.92 ; E3=112608.341 ; E4=889466.87 ; E5=1222604.76; 
u0=0; u1=320.75; u2=371.265; u3=474.57; u4=1333.766; u5=1563.716;

u_tst=331.046034536;
V0=1./rho0;
V1=1./rho1;
choc1=V0*sqrt((P1-P0)/(V0-V1))*time - u_tst*time

rho0= 7287;  E0= 3.1279315986251e-10; P0= 0; T0= 300; S0=2.6926731550247e-05; G0=3.1279315986251e-10; G0_diff=2.240097804912; C0=2740.5538450273;
rho1= 8076.00681801076; E1= 51443.002832343; P1= 7.6739859076668e9; T1= 404.14787213729; D1= 3283.1753189253; u1= 320.75848494575; S1= 16.007216133979; G1= 995194.05593801; G1_diff= 2.019367150304; C1= 3442.4648018356;
rho2= 8288.3703764792; E2= 77064.927809945; P2= 8.4780616126608e9; T2= 354.88066418806; D2= 3283.1753189253; u2= 371.26594786928; S2= 16.016695936291; G2= 1094267.3454656; G2_diff= 1.7223706495414; C2= 3430.2925202527;
rho3= 8518.2842765723; E3= 112608.3412405; P3= 11.353849215428e9; T3= 388.4355077502; D3= 3283.1753189253; u3= 474.56999745137; S3= 25.068708561843; G3= 1435750.5849554; G3_diff= 1.6743445292304; C3= 3581.513727052;
rho4= 10247.863012022; E4= 889466.82426922; P4= 44.866636800633e9; T4= 1965.6747151802; D4= 4616.3090031665; u4= 1333.7667144364; S4= 305.21638123923; G4= 4667656.3444466; G4_diff= 1.3493141278849; C4= 4624.8467142991;
rho5= 10609.008312322; E5= 1222604.7150191; P5= 56.903490402702e9; T5= 2231.3107527235; D5= 4993.8110485992; u5= 1563.7165440188; S5= 383.61978969721; G5= 5730325.2776706; G5_diff= 1.2127819931885; C5= 4829.4181283246;

V0=1./rho0;
V1=1./rho1;
choc1=V0*sqrt((P1-P0)/(V0-V1))*time - u_tst*time


lambda0=[1, 0, 0];

lambda1=[1,0,0];

lambda(101,:)=[1,0,0];
p_hug(101)= P1; u_hug(101)= u1; rho_hug(101)= rho1; E_hug(101)= E1; T_hug(101)= T1; S_hug(101)= S1; C_hug(101)= C1; G_hug(101)= G1_diff;   


% ZONE MIXTE beta/gamma (1/2)
% lambda=0.2
lambda(102,:)=[0.8, 0.2, 0];
num=102;
%p_hug(102)= 7.83711e+09; u_hug(102)= 331.046; rho_hug(102)= 8118.55; E_hug(102)= 56474.8;
P= 7837109295.2405; u= 331.046034536; rho= 8118.5453550428; V= 0.00012317477531598; E= 56474.765160832; T= 394.02227652445; S=16.007357604096; G=1021625.0939154; G_diff=1.9570238747402; D=1963.3946949723; C=3442.6037758954;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;   

% lambda=0.5 
num=103;
lambda(103,:)=[0.5, 0.5, 0];
%p_hug(103)= 8.07994e+09; u_hug(103)= 346.31; rho_hug(103)= 8182.28; E_hug(103)= 64111.2;
P= 8079944236.04; u= 346.31015149475; rho= 8182.2810751876; V= 0.00012221530778654; E= 64111.232454084; T= 379.07885359786; S=16.009038155771; G=1051485.1118034; G_diff=1.8680657939052; D=1967.2771487405; C=3440.3507108391;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;   

% lambda=0.8
lambda(104,:)=[0.2, 0.8, 0];
num=104;
%p_hug(104)= 8.31994e+09; u_hug(104)= 361.359; rho_hug(104)= 8245.95; E_hug(104)= 71850.5;
P= 8319940699.7071; u= 361.35917843614; rho= 8245.9493636798; V= 0.00012127166392806; E= 71850.5064501; T= 364.44976977374; S=16.013105491071; G=1080763.317928; G_diff=1.7803046192314; D=1970.0261716925; C=3435.2332470752;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;   
  
lambda(105,:)=[0, 1, 0];
p_hug(105)= P2; u_hug(105)= u2; rho_hug(105)= rho2; E_hug(105)= E2; T_hug(105)= T2; S_hug(105)= S2; C_hug(105)= C2; G_hug(105)= G2_diff;   

lambda(106,:)=[0,1,0];
p_hug(106)= P3; u_hug(106)= u3; rho_hug(106)= rho3; E_hug(106)= E3; T_hug(106)= T3; S_hug(106)= S3; C_hug(106)= C3; G_hug(106)= G3_diff;   

lambda(107,:)=[0, 1, 0];
p_hug(107)= P4; u_hug(107)= u4; rho_hug(107)= rho4; E_hug(107)= E4; T_hug(107)= T4; S_hug(107)= S4; C_hug(107)= C4; G_hug(107)= G4_diff;   


% ZONE MIXTE gamma/liquide (2/3)
% lambda=0.2
lambda(108,:)=[0, 0.8, 0.2];
num=108;
%p_hug(108)= 4.73658e+10; u_hug(108)= 1383.13; rho_hug(108)= 10326.1; E_hug(108)= 956520;
P= 47365759158.054; u= 1383.1267402852; rho= 10326.095622402; V= 9.6842023991191e-05; E= 956519.78984598; T= 2023.7388172682; S=321.74209642994; G=4965690.1012973; G_diff=1.3211817565948; D=4699.5227372271; C=4669.9996611123;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;   

% lambda=0.5
lambda(109,:)=[0, 0.5, 0.5];
num=109;
%p_hug(109)= 5.10185e+10; u_hug(109)= 1453.67; rho_hug(109)= 10437.2; E_hug(109)= 1.05658e+06;
P= 51018537544.691; u= 1453.6692554144; rho= 10437.177697807; V= 9.5811341816103e-05; E= 1056577.1520685; T= 2105.7330981251; S=345.65235546068; G=5291859.4850195; G_diff=1.2788079077289; D=4816.3011068108; C=4733.3837556223;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;

% lambda=0.8
lambda(110,:)=[0, 0.2, 0.8];
num=110;
%p_hug(110)= 5.45762e+10; u_hug(110)= 1520.7; rho_hug(110)= 10542; E_hug(110)= 1.15626e+06;
P= 54576171339.318; u= 1520.6957243034; rho= 10542.025726704; V= 9.4858429103141e-05; E= 1156257.7429573; T= 2182.5586666218; S=368.68221712734; G=5605149.5720574; G_diff=1.2373877250962; D=4925.0650514303; C=4792.2995940202;
p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;
 
lambda(111,:)=[0, 0, 1];
num=111;
p_hug(111)= P5; u_hug(111)= u5; rho_hug(111)= rho5; E_hug(111)= E5; T_hug(111)= T5; S_hug(111)= S5; C_hug(111)= C5; G_hug(111)= G5_diff;

lambda(112,:)=[0, 1, 0];
num=112;
%u_hug(112)=900; rho_hug(112)=1./0.000106217; E_hug(112)=405000; p_hug(112)=2.61173e+10;
%p_hug(num)= P; u_hug(num)= u; rho_hug(num)= rho; E_hug(num)= E; T_hug(num)= T; S_hug(num)= S; C_hug(num)= C; G_hug(num)= G_diff;
 
 
% simple choc
if(tst==101 || tst==106 || tst==107 || tst==108 || tst==109 || tst==110 || tst==111 || tst==112)
  V0=1./rho0;
  V=1./rho_hug(tst);
  choc1=V0*sqrt((p_hug(tst)-P0)/(V0-V))*time - u_hug(tst)*time

  x=[a,-choc1-dx,-choc1+dx,choc1-dx,choc1+dx,b];

  Uana=-sign(x)*((u0-u_hug(tst))*(abs(x)<=choc1) + u_hug(tst));
  Pana=(p_hug(tst)-P0)*(abs(x)<=choc1) + P0;
  RHOana=(rho_hug(tst)-rho0)*(abs(x)<=choc1) + rho0;
  Eana=(E_hug(tst)-E0)*(abs(x)<=choc1) + E0;
  Sana=(S_hug(tst)-S0)*(abs(x)<=choc1) + S0;
  Cana=(C_hug(tst)-C0)*(abs(x)<=choc1) + C0;
  Tana=(T_hug(tst)-T0)*(abs(x)<=choc1) + T0;
  Gana=(G_hug(tst)-G0_diff)*(abs(x)<=choc1) + G0_diff;
  for i=1:3
    LAMBDA(:,i)=(lambda(tst,i)-lambda0(i))*(abs(x)<choc1) + lambda0(i);
  endfor
else 
% double choc
  V0=1./rho0
  V1=1./rho1
  V=1./rho_hug(tst)
  choc1=V0*sqrt((P1-P0)/(V0-V1))*time - u_hug(tst)*time
  choc2=V1*sqrt((p_hug(tst)-P1)/(V1-V))*time  + u1*time - u_hug(tst)*time
  
  x=[a,-choc1-dx,-choc1,-choc1+dx,-choc2-dx,-choc2+dx,choc2-dx,choc2+dx,choc1-dx,choc1,choc1+dx,b];

  Uana=(-u_hug(tst) + u1)*(abs(x)<choc2) - (u1)*(abs(x)<choc1) + u_hug(tst);
  Uana=-Uana.*sign(x);
  Pana=(p_hug(tst)-P1)*(abs(x)<choc2) + (P1-P0)*(abs(x)<choc1) + P0;
  RHOana=(rho_hug(tst)-rho1)*(abs(x)<choc2) + (rho1-rho0)*(abs(x)<choc1) + rho0;
  Eana=(E_hug(tst)-E1)*(abs(x)<choc2) + (E1-E0)*(abs(x)<choc1) + E0;
  Tana=(T_hug(tst)-T1)*(abs(x)<choc2) + (T1-T0)*(abs(x)<choc1) + T0;
  Sana=(S_hug(tst)-S1)*(abs(x)<choc2) + (S1-S0)*(abs(x)<choc1) + S0;
  Cana=(C_hug(tst)-C1)*(abs(x)<choc2) + (C1-C0)*(abs(x)<choc1) + C0;
  Gana=(G_hug(tst)-G1_diff)*(abs(x)<choc2) + (G1_diff-G0_diff)*(abs(x)<choc1) + G0_diff;
  for i=1:3
    LAMBDA(:,i)=(lambda(tst,i)-lambda1(i))*(abs(x)<choc2) + (lambda1(i)-lambda0(i))*(abs(x)<choc1) + lambda0(i);
  endfor
endif


endfunction