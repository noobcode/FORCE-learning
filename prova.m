m=32;
dt=0.04;
target_length= 2301;

sigma = 0.1 * target_length/(m*1000); % standard deviation of gaussian, 
                                      % longer signals have wider pulses
a = 1;                                % amplitude
mu = (target_length/2)*dt/1000;       % mean, centred at half signal

t = (1:1:(target_length/dt))*dt;

%num = (t-mu).^2;
%den = 2*sigma^2;
%f = a * exp(-num/den);

gaussian_pulse = @(t,a,mu,sigma) a * exp(-(t-mu).^2/2*sigma^2);

y = pulstran(t,1:1/m:target_length,gaussian_pulse(t,a,mu,sigma)); %uses a sample rate of fs.

figure;
plot(t,y)