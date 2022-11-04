

function [Xoci,Uoci,Poci,Voci,Eoci,Coci,Soci,Dfoci,LAMBDAoci]=sol_oci_etain(T, u_init)
     
system("sh compil_onde_compressible.sh");  % CHANGER LE .C EN FCT DE TST
system(["a.exe ",num2str(u_init)]);
     
% Onde de compression isentropique
Voci=load("fichiers/Voci.txt");
Eoci=load("fichiers/Eoci.txt");
Poci=load("fichiers/Poci.txt");
Toci=load("fichiers/Toci.txt");
XAoci=load("fichiers/XAoci.txt");
XBoci=load("fichiers/XBoci.txt");
XCoci=load("fichiers/XCoci.txt");
Coci=load("fichiers/Coci.txt");
Soci=load("fichiers/Soci.txt");
Dfoci=load("fichiers/Dfoci.txt");
Uoci=load("fichiers/Uoci.txt");

lambda(1,:)=XAoci;
lambda(2,:)=XBoci;
lambda(3,:)=XCoci;

Xoci=(U+C)*T;