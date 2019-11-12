function [y] = gaussian_pulse_train(target_length, dt, m)
 %% 
 pulse_len = target_length/m/dt;
 t = (1:1:pulse_len)*dt/1000;
 
 sigma = 0.1 * target_length/(m*1000); % standard deviation of gaussian, 
                                   % longer signals have wider pulses
 a = 1;                            % amplitude
 mu = (pulse_len/2)*dt/1000;       % mean, centred at half signal
 
 num = (t-mu).^2;
 den = 2*sigma^2;
 y = a * exp(-num/den);
 
 size(y), target_length/dt, round(target_length/dt)
 
 temp = y;
 for i=1:m-1
    y = [y,temp]; 
 end
 
 %%
 do_plot = false;
 if do_plot
    tt = (1:1:(pulse_len*m))*dt;
    size(tt), size(y)
    figure;
    plot(tt,y)
 end
end