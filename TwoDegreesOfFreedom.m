% Two Degrees Of Freedom Model: State Variables and Finite Differences

clear;clc;close all;

% Parametros 
fn   = 0.15;
wn   = fn*2*pi;
m2   = 726*1e3;
m1   = m2*100/0.24;
c1   = 0;
c2   = 2*0.047*wn*m2;
c3   = 0;
k1   = m1*wn^2;
k2   = m2*wn^2;
k3   = 0;
Tsim = 1200;
dt   = 1e-2;
t = dt:dt:Tsim;

% Entrada
indini = 1;

%=================================%

% f1 = zeros(Tsim/dt,1);           % Delta
% f1((indini:(indini))/dt)=1000e9; % Delta

f1 = 1000e3*sin(2*pi*fn*0.7*t); % Sinusoide

f2 = 0;

%=================================%

% Condiciones Iniciales

v = [0;0;0;0];
V = zeros(Tsim/dt,4);

% Solución Numerica
B = [0,1,0,0;-(k1+k2)/m1,-(c1+c2)/m1,k2/m1,c2/m1;0,0,0,1;k2/m2,c2/m2,-(k2+k3)/m2,-(c2+c3)/m2];
Ms = (eye(4)-dt*B);
cont = 0;

for ii = dt:dt:Tsim
    cont = cont+1;
    F = [0;f1(cont)/m1;0;f2(min(cont,end))/2];
    v = Ms\(v+dt*F);
    V(cont,:) = v;
end

x1 = V((indini+100*0)/dt:end,1);
x2 = V((indini+100*0)/dt:end,3);

% ===================================================== %

% Plots

figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,2,1:2),plot(t,V(:,1),'k','LineWidth',1.2),hold on,plot(t,V(:,3),'r','LineWidth',1.2),xlim([0,500]),xlabel('Time [Seg]'),legend({'x_1','x_2'}),set(gca,'FontSize',12),grid on
title("Solución Numerica")

fs = 1/dt;
f  = (0:(length(x1)-1))*fs/length(x1);

subplot(323),plot(f,abs(fft(x1)),'k','LineWidth',1.2),hold on,plot(fn*[1,1],[0,max(max(abs(fft(x1))))],'--b','LineWidth',1.2)
xlim([0,2*fn]),ylim([0,max(max(abs(fft(x1))),max(abs(fft(x2))))]),grid on
xlabel('Frequency [Hz]'),ylabel('|DFT\{x1\}|'),set(gca,'FontSize',12)
subplot(325),plot(f,abs(fft(x2)),'r','LineWidth',1.2),hold on,plot(fn*[1,1],[0,max(max(abs(fft(x1))),max(abs(fft(x2))))],'--b','LineWidth',1.2)
xlim([0,2*fn]),ylim([0,max(max(abs(fft(x1))),max(abs(fft(x2))))]),grid on
xlabel('Frequency [Hz]'),ylabel('|DFT\{x2\}|'),set(gca,'FontSize',12)

subplot(324),plot(f,rad2deg(unwrap(angle(fft(x1)))),'k','LineWidth',1.2)
xlim([0,2*fn]),xlabel('Frequency [Hz]'),ylabel('\angle DFT\{x1\}'),set(gca,'FontSize',12),grid on
subplot(326),plot(f,rad2deg(unwrap(angle(fft(x2)))),'r','LineWidth',1.2)
xlim([0,2*fn]),xlabel('Frequency [Hz]'),ylabel('\angle DFT\{x1\}'),set(gca,'FontSize',12),grid on
