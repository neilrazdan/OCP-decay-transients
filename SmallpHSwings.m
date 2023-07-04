clear all; close all; clc
%% Import data
% CHOOSE SHEET;
% Sheet = '3 M NaClO4, Stirred'; check = 1;
% Sheet = '3 M NaClO4, Quiescent'; check = 1;
Sheet = '3 M LiClO4, Stirred'; check = 1;
% Sheet = '3 M NaClO4, Nafion, Stirred'; check = 2;
M = readmatrix('SmallpHSwings.xlsx','Sheet',Sheet); % Import Excel Sheet
d = M(:,3:4); J = d(1,1); d = d(2:end,:); % Choose Data Set (1:2 = 10 mA cm2, 3:4 = 30 mA cm2, 5:6 = 53 mA cm2)
%% Process data
d = d(100:end,:); % cut off first 20 ms corresponding to capacitive decay
time = d(:,1); time = time - time(1); E = d(:,2); %  redefine time zero
%% Parameters
beta = 1 - 10.^((E(1)*1000)/59); % (C0 - Cinf)/Cinf
D = 9.312*10^-5; % cm^2 s^-1
C = 0.1; % mol/L; bulk hydronium concentration
C0 = C*(1-beta); % mol/L; steady state interfacial hydronium concentration
F = 96485; % C/mol
delBL = beta*(D*C/(J/F)); % cm
if check == 2
    tau = 2.4;
else
    tau = delBL^2/D;
end
X = linspace(0,3,10); % dimensionless distance vector
T = logspace(-4,4,100); % dimensionless time vector
omega = logspace(-10,10,10^4); % dummy variable for Fourier integral
n = [0:1:10^4]; m = n+0.5; % dummy variable for Fourier sum
%% Array Prep
C = zeros(length(X),length(T)); % empty matrix to hold finite solution
CC = zeros(length(X),length(T)); % empty matrix to hold infinite solution
%% Numerical Integration of Analytical Infinite Domain Solution
if check ~= 2
    for i = 1:length(X)
        for j = 1:length(T)
            x = X(i); t = T(j);
            u = 2*beta/pi.*(cos(omega) - 1)./omega.^2.*cos(omega.*x).*exp(-omega.^2.*t); % Fourier integral
            U = trapz(omega,u); % Integrate Fourier integral
            C(i,j) = U; % Allocate solution into concentration matrix
        end
    end
end
%% Numerical Integration of Analytical Finite Domain Solution
for i = 1:length(X)
    for j = 1:length(T)
        x = X(i); t = T(j);
        u = -2*beta.*cos(m*pi*x)./(m*pi).^2.*exp(-(m*pi).^2.*t); % Fourier integral
        U = sum(u); % Integrate Fourier integral
        CC(i,j) = U; % Allocate solution into concentration matrix
    end
end
%% Concentration Transient
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(time,10.^((E*1000)/59)-1,'k-','linewidth',1);
plot(T*tau,C(1,:),'k--','linewidth',2)
plot(T*tau,CC(1,:),'r-','linewidth',1)
xlabel('t / (\delta_{BL}^{2}/D)')
ylabel('(C_{0} - C_{\infty}) / C_{\infty}')
title(Sheet)
legend(['j = ' num2str(J) ' mA cm^{-2}'],'Quiescent Solution', 'Stirred Solution')
legend boxoff
legend('location','southeast')
xlim([0 4]);
%% Potential Transient
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(time),E*1000,'k-','linewidth',1);
plot(log10(T*tau),59*log10(C(1,:)+1),'k--','linewidth',2)
plot(log10(T*tau),59*log10(CC(1,:)+1),'r-','linewidth',1)
xlabel('log_{10} t / s')
ylabel('\Delta E / mV')
title(Sheet)
legend(['j = ' num2str(J) ' mA cm^{-2}'],'Quiescent Solution', 'Stirred Solution')
legend boxoff
legend('location','southeast')
