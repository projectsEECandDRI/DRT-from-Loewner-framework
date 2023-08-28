% Ex_RC.m             
%
% Extraction of the tempered distribution of relaxation times (DRT) from 
% impedance spectra related to an RC circuit. The circuit is constituted
% by two RC units and a resistor connected in parallel.

% The numbers of the algorithm steps and the equations mentioned in the 
% comments refer to the supplementary information document of the paper 
% "Sorrentino A., B. Patel,I. V. Gosea, A. Anthoulas, T. Vidakovic-Koch,
% Determination of the distribution function of relaxation times through
% Loewner framework: a versatile and regularization free approach".

% Main steps of the algorithm:
% -Creating an impedance data set related to the RC circuit <=> (Step 1)
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


%% Creating an impedance spectra data set of RC circuit

% Specifications
R1 = 200; % resistance [Ohm]
R2 = 100; % resistance [Ohm]
R0 = 70;  % resistance [Ohm]

C1 = 2.5e-3; % capacitance [F]
C2 = 1e-4;   % capacitance [F]

tau1 = R1*C1;  % time constants [s]
tau2 = R2*C2;  % time constants [s]

f = logspace(-2,2,41);  % frequency range of the impedance spectra [Hz]

% Calculation of the impedance dataset 
for i=1:length(f)
w=2*pi*f(i);

Z_RC(i)=1/((1/R1)+1j*w*C1) +1/((1/R2)+1j*w*C2) +R0;

end

   
%% Construction of the Loewner framework matrices and vectors (Steps 2-4)
% NOTE: Option `'real is set if one deals with data sets with imaginary
% entries
[L,Ls,V,W] = Loewner_Framework(1i*2*pi*f,Z_RC,'real');


%% Construction of the reduced state space model (Steps 5-6)
[Ek, Ak, Bk, Ck, s_LLs1] = state_space_mod(L,Ls,V,W);

% Transfer function interpolating the data
for ii=1:length(f)
    Hk(ii)=Ck*((1i*2*pi*f(ii)*Ek-Ak)\Bk);
end

% Residula between the data and the obtained impedance transfer function
ReZ_rel=abs(((real(Z_RC)-real(Hk))./(abs(Z_RC))))*100;
ImZ_rel=abs(((imag(Z_RC)-imag(Hk))./(abs(Z_RC))))*100;
   

%% Determination fo the resistance and time constants (Steps 7-9)

% Determination poles and residues of the impedance transfer function
[R_i,tau_i] = DRT_values(Ek,Ak,Bk,Ck);
 

%% Plotting diagrams 

h=figure(1);hold on;
set(h,'Units','centimeters');
set(h,'Position',[5 0 20 20]);

% figure 1a: Computed Nyquist plot vs Data
subplot('Position',[0.165 0.6 0.3 0.3]);
plot(real(Z_RC),-imag(Z_RC),'ok','markersize',5);
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
set(gca,'YLim',[0 300]);
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

