% Ex_Zarc.m             
%
% Extraction of the tempered distribution of relaxation times (DRT) from 
% impedance spectra related to an Zarc circuit. An analytical DRT is also 
% caluclated in order to compare it with an analytical one

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

R1 = 100; % resistance [Ohm]

Y1 = 1e-4; % capacitance [Ohm]
 
n1 = 0.85; % 

tau1 = (R1*Y1)^(1/n1); % time constant [s]

f = logspace(-1,5,61); % frequency range of the impedance spectra [Hz]

% Calculation of the impedance dataset 
for ii=1:length(f)
    
Q1(ii) = 1/(Y1*(1j*2*pi*f(ii))^n1);

Z(ii) = 1/((1/R1)+Q1(ii)^-1);

end


%% Construction of the Loewner framework matrices and vectors (Steps 2-4)
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*f,Z,'real');

  
%% Construction of the reduced state space model (Steps 5-6)
[Ek, Ak, Bk, Ck, s_LLs1] = state_space_mod(L,Ls,V,W);

% Transfer function interpolating the data
for ii=1:length(f)
    Hk(ii)=Ck*((1i*2*pi*f(ii)*Ek-Ak)\Bk);
end

% Residula between the data and the obtained impedance transfer function
ReZ_rel=abs(((real(Z)-real(Hk))./(abs(Z))))*100;
ImZ_rel=abs(((imag(Z)-imag(Hk))./(abs(Z))))*100;
   

%% Determination fo the resistance and time constants (Steps 7-9)

% Determination poles and residues of the impedance transfer function
[R_i,tau_i] = DRT_values(Ek,Ak,Bk,Ck);
 

%% Analytical DRT of Zarc circuit (see Section 3.3 in the paper)

% parameters 
tau = tau_i(1:end); % time vector [s]
x   = log(tau1./tau); % time vector in logarithmic scale
d_tau = 0.01 -tau_i(end); 
tau   = [tau;20];


% Analytical DRT function of Zarc circuit (see equation 8 in the
% Supplementary Information)
fun   = @(x) (R1./(2*pi)).*sin((n1)*pi)./(cosh(n1.*x) +cos((n1)*pi));

% Integration fo analytical DRT over log(tau) (Equation 19 in the paper)
yy= length(tau); 
for i=1:length(tau)
      
q = integral(fun,log(tau1/tau(1)),log(tau1/tau(yy)));

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
qq_c_end=qq_c(end-1)-(qq_c(end-1)-qq_c(end))*((tau(end-3)-tau(end-1))/(tau(end-3)-tau(end-2)));
qq_c=[qq_c;qq_c_end];    


%% Plotting diagrams 

h=figure(1);hold on;
set(h,'Units','centimeters');
set(h,'Position',[5 0 20 20]);

% figure 1a: Computed Nyquist plot vs Data
subplot('Position',[0.165 0.6 0.3 0.3]);
plot(real(Z),-imag(Z),'ok','markersize',5);
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
set(gca,'XLim',[1e-4 1e0]);
xlabel('\tau [s]')
ylabel('R_i [\Omega]')
grid on
hold off;
legend('Location','NorthEast')


