clear all
close all
%%
T = 3000;  % Total time in ms
dt = 0.04; % Integration time step in ms 
nt = round(T/dt); % Time steps
imin = round(1000/dt);
icrit = round(2000/dt);
%%
A = 5; % wave amplitude
f = 3; % oscillation frequency (Hz)
omega =  2*pi*f; % angular frequency, the rate of change of the function argument (radians per second)
%%
step_per_second = round(1000/dt);
step_per_cycle = round(step_per_second / f); % number of iterations to make a cycle (2pi)
%omega_training = (f) *  omega;
train_every = round(1 * step_per_cycle); % train every 1/2 pi

phase = 1/2*pi; % phase
phase_training = round(1/4 * step_per_cycle);

%%
x = (1:1:nt) * dt/1000;
xx = (1:step_per_cycle:nt) * dt/1000;
xx2 = (1:step_per_cycle:nt) + phase_training; % ! 
xx3 = (1:train_every:nt);

y = A * sin(omega*x);
y2 = A * sin(omega*x + phase);

plot(x,y), hold on
plot(xx, 0*xx, 'x'), hold on
plot(xx2  * dt/1000, A*sin(omega*xx2 * dt/1000), 'o'), hold on
%plot(xx3 * dt/1000, 0*xx3 * dt/1000, 'o', 'LineWidth', 2)
legend('original', '2pi', '2pi+phase')

%plot(x,y2), hold on

j = [];
for i=1:1:nt
   if i > imin && i < icrit
       if mod(i-imin, step_per_cycle) == 1 + phase_training
          j(end+1) = i; 
       end
   end
end

plot(j*dt/1000, 0*j, 'o')
