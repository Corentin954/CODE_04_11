% Solution exacte AVERAGE
% n : nx+1

function [ RHObar,Ubar,Pbar]=sol_acous(Xd,Xu,Xc,n,fac_eps)
     
   format long % permet d'avoir une dizaine de chiffre significatif 
     
    
   % Planar Kidder 
   kappa=2*pi;
   rho0=1.0; p0=5./7;
   
   rho_bar=@(a,b)  rho0 - (fac_eps/(kappa*(b-a)))*(cos(kappa*b) - cos(kappa*a));
   u_bar=@(a,b)  -(fac_eps/(kappa*(b-a)))*(cos(kappa*b) - cos(kappa*a));
   p_bar=@(a,b)  p0 - (fac_eps/(kappa*(b-a)))*(cos(kappa*b) - cos(kappa*a));

   rho_ana=@(x)  rho0 + fac_eps*sin(kappa*x);
   u_ana=@(x)  fac_eps*sin(kappa*x);
   p_ana=@(x)  p0 + fac_eps*sin(kappa*x);
   
   %%  RHO et P
    RHObar=zeros(1,n-1); 
    Pbar=zeros(1,n-1); 
    for i=1:n-1 
      %RHObar(i)=rho_bar(Xd(i),Xd(i+1));
      %Pbar(i)=p_bar(Xd(i),Xd(i+1));
      
      RHObar(i)=rho_ana(Xc(i));
      Pbar(i)=p_ana(Xc(i));
    endfor

   
    %%  U
    Ubar=zeros(1,n); 
    for i=1:n
      %Ubar(i)=u_bar(Xu(i),Xu(i+1));
      
      Ubar(i)=u_ana(Xd(i));
    endfor
    
    
endfunction