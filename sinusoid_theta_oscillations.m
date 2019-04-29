T = 5000;  %Total time in ms
dt = 0.04; %Integration time step in ms 
nt = round(T/dt); %Time steps

A = 5; % wave amplitude
f = 1; % oscillation frequency (Hz)
omega =  2*pi*f; % angular frequency, the rate of change of the function argument (radians per second)
phase = 0; % phase
 
x = (1:1:nt) * dt/1000;
y = A * sin(omega*x + phase);
 
plot(x,y)
 