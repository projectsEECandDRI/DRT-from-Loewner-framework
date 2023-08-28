function [Ek, Ak, Bk, Ck,s_LLs1] = state_space_mod(L,Ls,V,W)

% rank of the Loewner pencil [L Ls]
k_LLs1 = rank([L Ls]);
k = k_LLs1;

% Singular value decomposition of the Loewner pencil [L Ls]
[Y_LLs1,svd_LLs1, X_LLs1] = svd([L Ls],'econ');
s_LLs1 = diag(svd_LLs1);

% Singular value decomposition of the Loewner matrix
[Y_L,svd_L,X_L]=svd(L);
s_L=diag(svd_L);

Yk = Y_L(:,1:k);
Xk = X_L(:,1:k);

% Reduced state space model interpolating the data
Ek=-Yk.'*L*Xk;
Ak=-Yk.'*Ls*Xk;
Bk=Yk.'*V;
Ck=W.'*Xk;

end