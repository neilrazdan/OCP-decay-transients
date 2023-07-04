clear all; clc
%% Import data
% CHOOSE SHEET;
% Sheet = '3 M NaClO4, Stirred'; 
Sheet = '3 M LiClO4, Stirred'; 
d = readmatrix('LargepHSwings.xlsx','Sheet',Sheet); % Import Excel Sheet
J = d(1,1); pH_initial = d(2,1); d = d(3:end,:); % Extract J, pH_initial
%% Process data
d = d(100:end,:); % cut off first 20 ms
time = d(:,1); time = time - time(1); E = d(:,2); % redefine time zero
time = time(400:end); E = E(400:end); % cut off 80 ms of eta_CT and eta_chem convolution
%% pH swing
KW = 10^-14; % water auto-ionization eqm
CH_bulk = 10^-1; % M; bulk hydronium concentration
COH_bulk = KW/CH_bulk; % M; bulk hydroxide concentration
COH_0 = 10^(-(14-pH_initial)); % M; steady state interfacial hydroxide concentration
CH_0 = KW/COH_0; % M; steady state interfacial hydronium concentration
F = 96485; % C/mol
%% Boundary Layer
alpha = (COH_0 - sqrt(KW))./(CH_bulk + COH_0 - 2*sqrt(KW)); % delta_n/delta_BL
DH = 9.312*10^-5; DOH = 5.26*10^-5; % cm^2/s
R = DOH/DH; % Ratio of diffusion coefficeints
krecomb = 10^9; % water recomb
k = KW*krecomb; % k dissociation
del = (COH_0 - sqrt(KW))*DOH/(J/F); % delta_n (called del herein)
delBL = del/alpha; % cm
eps = 1/(sqrt(KW).*DH/delBL^2/k); % small parameter
tau = delBL^2/DH; % timescale for non-dimensionalization
%% Dimensionless Parameters
T_bulk = CH_bulk/sqrt(KW);
T_0 = CH_0/sqrt(KW);
phi_bulk = 1./T_bulk;
%% Time and Space Arrays
x = linspace(0,1,500)*eps^(1/2); % dimensionless distance vector
t = logspace(-3,3,2000)*eps; % dimensionless time vector
%% PDEPE
P = pde(x,t,R,eps,T_0,T_bulk,phi_bulk,alpha);
%% Post Processing
CH = sqrt(KW)*P(:,:,1);
COH = sqrt(KW)*P(:,:,2);
VH = 59*log10(CH/CH_bulk);
VOH = -59*log10(COH/COH_bulk);
X = x/eps^0.5;
T = t/eps;
%% Potential Transient
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(time),E*1000,'k-','linewidth',1);
plot(log10(T*tau),VH(:,1),'k--','linewidth',3)
plot(log10(T*tau),VOH(:,1),'r-','linewidth',1)
xlabel('log_{10} t / s')
ylabel('\Delta E / mV')
title(Sheet)
legend(['j = ' num2str(J) ' mA cm^{-2}'],'pH-based OCP', 'pOH-based OCP')
legend boxoff
legend('location','southeast')

