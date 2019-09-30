close all
clear all

T = 1000;
dt = 0.04;
tt = (1:1:round(T/dt))*dt/1000;

zx1 = sin(2*pi*1*tt);
zx2 = sqrt(tt); %sin(2*pi*f*tt1); %sqrt(tt2);
zx3 = tt.^2 - sqrt(tt);
zx4 = 1./exp(-tt) - tt.^3;

figure;
plot(tt, zx1, 'LineWidth', 2, 'Color' , 'k')
xlabel('Time(s)')
ylabel('x')
ax = gca;
ax.FontSize = 16;

figure;
plot(tt, zx2, 'LineWidth', 2, 'Color' , 'k')
xlabel('Time(s)')
ylabel('x')
ax = gca;
ax.FontSize = 16;

figure;
plot(tt, zx3, 'LineWidth', 2, 'Color' , 'k')
xlabel('Time(s)')
ylabel('x')
ax = gca;
ax.FontSize = 16;
