% Ex_GE.m             
%
% Extraction of the tempered distribution of relaxation times (DRT) from 
% impedance spectra related to the Gerischer element. An analytical DRT is 
% also caluclated in order to compare it with an analytical one.

% The numbers of the algorithm steps and the equations mentioned in the 
% comments refer to the supplementary information document of the paper 
% "Sorrentino A., B. Patel,I. V. Gosea, A. Anthoulas, T. Vidakovic-Koch,
% Determination of the distribution function of relaxation times through
% Loewner framework: a versatile and regularization free approach".

% Main steps of the algorithm:
% -Creating an impedance data set related to the Zarc circuit <=> (Step 1)
% -Construction of the Loewner matrix L, shifted Loewner marix, vectors V
%  and W using Loewner_Framework.m <=> (Steps 2-4)
% -Determination of the reduced version of the state space model Ek, Ak,
%  Bk, Ck interpolating impedance data <=> (Steps 5-6)
% -Determination of the resistances R_i and time constants time_i from the 
%  poles and residues of the impedance transfer function using the fuction
%  find_poles_residues.m <=> (Steps 5-6) 
% -Plots of the diagrams
%  
%
% This example script is part of the DRT using Loewner framework package
% (DRT_LF.zip) 
% Last revised: 23.08.2023. 
%

clear all; 
close all;


%% Creating an impedance spectra data set of Zarc circuit

% Parameters
Z0   = 100; % DC resistance [Ohm]
tau0 = 0.01; % characteristic time constant of the process [s]

f = logspace(-3,5,150); % frequency range of the impedance spectra [Hz]

% Calculation of the impedance dataset 
w   = 2*pi*f;
Z_G = (Z0./(sqrt(1+1i*w*tau0)));


%% Construction of the Loewner framework matrices and vectors (Steps 2-4)
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*f,Z_G,'real');

  
%% Construction of the reduced state space model (Steps 5-6)
% NOTE: In this case, since there is no ohmic immediate response (R0=0),it
% is better to perform the SVD of the Loenwer pencil [L Ls] to determine
% the state space model. In case the state_space_mod.m function (
% based on SVD of Loewner matrix) is used,the algorithm would still give 
% the state space model interpolating the data and the correct DRT.  
% However, some negligible impulse functions(different order of magnitude 
% smaller of the other ones) would be also obtained which do not compromise
% the use of the algorithmfor practical applications in any case. These
% additional quantities are due to numerical approximations in the 
% interpolation of the transfer function at high frequencies in case of 
% ohmic resistance equal to zero. This is also mentioned in the
% Supplementary Information

% Singular value decomposition of the Loewner pencils [L Ls] and [L; Ls]
[Y_LLs1,svd_LLs1, X_LLs1] = svd([L Ls],'econ'); 
s_LLs1=diag(svd_LLs1);

[Y_LLs2,svd_LLs2, X_LLs2] = svd([L;Ls],'econ');
s_LLs2=diag(svd_LLs2);

% rank of the Loewner pencil [L Ls] and model order
k_LLs1 = rank([L Ls]);
k = k_LLs1;

% Reduced state space model interpolating the data
Yk = Y_LLs1(:,1:k);
Xk = X_LLs2(:,1:k);

Ek=-Yk'*L*Xk;
Ak=-Yk'*Ls*Xk;
Bk=Yk'*V;
Ck=W.'*Xk;


% Transfer function interpolating the data
for ii=1:length(f)
    Hk(ii)=Ck*((1i*2*pi*f(ii)*Ek-Ak)\Bk);
end

% Residula between the data and the obtained impedance transfer function
ReZ_rel=abs(((real(Z_G)-real(Hk))./(abs(Z_G))))*100;
ImZ_rel=abs(((imag(Z_G)-imag(Hk))./(abs(Z_G))))*100;
   

%% Determination fo the resistance and time constants (Steps 7-9)

% Determination poles and residues of the impedance transfer function
[R_i,tau_i] = DRT_values(Ek,Ak,Bk,Ck);
 

%% Analytical DRT of Gerischer element (see Section 3.4 in the paper)

% parameters 
tau = tau_i(1:end); % time vector [s]
tau = [tau;0.01];


% Analytical DRT function of Gerischer element (see equation 10 in the
% Supplementary Information)
fun   = @(x) real((sqrt(2).*Z0./pi).*(sqrt((sqrt(1-(tau0.*exp(x)).^2)-1)./...
    (1-(tau0.*exp(x)).^2))));

% Integration fo analytical DRT over log(tau) (Equation 19 in the paper)
yy= length(tau); 
for i=1:length(tau)
      
q = integral(fun,log(1/tau(1)),log(1/tau(yy)));

yy=yy-1;
qqq(i,1)=q;
end

% Calculation of integrals of the DRT qq_c in the intervals [(log(tau_i+1) 
% +log(tau_i))/2, (log(tau_i) +log(tau_i-1))/2] (See equation 20 in the 
% paper)
qq=qqq(1:end-1)-qqq(2:end);
qq=-qq(end:-1:1);
tau_aa=(tau(1:end-1)+tau(2:end))./2;

for i=1:length(qq)-2

    qq_c(i)=qq(i)+(qq(i+1)-qq(i))*((tau_aa(i)-tau(i+1))/...
    (tau_aa(i)-tau_aa(i+1)));
    
end
qq_c=qq_c(:);
%
qq_c_end=qq_c(end-1)-(qq_c(end-1)-qq_c(end))*((tau(end-3)-tau(end-1))/...
         (tau(end-3)-tau(end-2)));
qq_c=[qq_c;qq_c_end];    


%% Plotting diagrams 

h=figure(1);hold on;
set(h,'Units','centimeters');
set(h,'Position',[5 0 20 20]);

% figure 1a: Computed Nyquist plot vs Data
subplot('Position',[0.165 0.6 0.3 0.3]);
plot(real(Z_G),-imag(Z_G),'ok','markersize',5);
hold on
plot(real(Hk),-imag(Hk),'color','red','LineWidth',2)
axis equal; 
grid on;
xlabel('Re(Z) [\Omega]')
ylabel('-Im(Z) [\Omega]')
hold on
legend('Data','Loewner model','Location','NorthWest')
annotation(h,'textbox','String',{'(a)'},'FontSize',10,...
    'LineStyle','none',...
    'HorizontalAlignment','right','VerticalAlignment','top','Margin',10.0,...
    'Position',[0.16 0.61 0.3 0.3]);

% figure 1b: Residuals between the data and the computed transfer function 
subplot('Position',[0.57 0.6 0.35 0.3]);
plot(f,ReZ_rel,':+b','markersize',5)
hold on
plot(f,ImZ_rel,'--ob','markersize',5)
set(gca,'xscal','log')
grid on;
xlabel('Frequency [Hz]')
ylabel('Residual [%]')
legend('Real','Imag','Location','SouthWest')
annotation(h,'textbox','String',{'(b)'},'FontSize',10,...
    'LineStyle','none',...
    'HorizontalAlignment','right','VerticalAlignment','top','Margin',10.0,...
    'Position',[0.63 0.61 0.3 0.3]);

% figure 1c: Tempered distribution of relaxation times
subplot('Position',[0.165 0.2 0.75 0.28]);
stem(tau_i,R_i,'-.o','filled','markersize',5,'LineWidth',2) 
hold all
set(gca,'xscal','log')
%title('Gain vs Time')
xlabel('\tau [s]')
ylabel('R_i [\Omega]')
grid on; 
hold off;
%legend('Location','NorthWest')
annotation(h,'textbox','String',{'(c)'},'FontSize',10,...
    'LineStyle','none',...
    'HorizontalAlignment','right','VerticalAlignment','top','Margin',10.0,...
    'Position',[0.62 0.18 0.3 0.3]);   

% Figure 2 Singular values of the Loewner pencil [L Ls]
h=figure(2);hold on;
set(h,'Units','centimeters');
set(h,'Position',[5 0 20 20]);
subplot('Position',[0.165 0.3 0.56 0.28]);
semilogy(s_LLs1./s_LLs1(1),':*','MarkerSize',5,'DisplayName','[L Ls]')
xlabel('Order of the model [-]')
ylabel('Singular value of the Loewner pencil')
hold on
legend('Location','Best')
axis tight;
grid on;

% Figure 3 Comparison between analytical and numerical DRT
h=figure(3);hold on;
set(h,'Units','centimeters');
set(h,'Position',[5 0 20 20]);
subplot('Position',[0.165 0.2 0.75 0.28]);
stem(tau_i,R_i,'-.o','filled','markersize',5,'LineWidth',2,...
'DisplayName','DRT from Loewner model')
hold all
stem(tau(2:end-1),qq_c,'-.om','markersize',5,'LineWidth',2,...
'DisplayName','Analytical DRT','color','red')
hold all
set(gca,'xscal','log')
xlabel('\tau [s]')
ylabel('R_i [\Omega]')
grid on
hold off;
legend('Location','NorthEast')


