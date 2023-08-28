function [R_i,tau_i] = DRT_values(E,A,B,C)

% Eigenvalue decomposition of 
[U,T] = eig(A,E); % eigenvalue decomposition

% Determination of the poles
pol   = diag(T);     

% Determination of the residues
Bt = U\(E\B); 
Ct = C*U;
res = Bt.*Ct.';   

% Determination resistances R_i and related time constants values tau_i 
% constituting the DRT
R_i   = (-res./pol);
tau_i = abs(-1./pol);
 
 
end
