function [hdts1, hdts2, hdts3, hdts4] = load_four_hdts_same_length(m, target_length, compression, dt)
    %{
    creates 4 HDTS (sinusoid, rectangular, isoscele triangle, rectangular
    triangle), with 'm' pulses, and length 'target_length'.

    You can also create the compressed/dilated version by specifying
    'compression' to be different than 1.
    %}
    %%
    T = round(target_length/compression);
    tt = (1:1:round(T/dt))*dt;
    
     %% create pulse trains
    temp1 = abs(sin(m*pi*tt/T));
    
    width = (T/dt)/m; % pulse width
    temp2 = pulstran(tt, 1:width:T/dt, @rectpuls, width);
    
    x1 = sawtooth(2*pi*m*tt/T, 0.5);
    temp3 = (x1 - min(x1)) / ( max(x1) - min(x1) );
    
    %temp4 = gaussian_pulse_train(T4, dt, m);
    x2 = sawtooth(2*pi*m*tt/T);
    temp4 = (x2 - min(x2)) / ( max(x2) - min(x2) );
    
    %% create HDTS
    for qw = 1:1:m
        hdts1(qw,:) = temp1 .* (tt < qw*T/m) .* (tt > (qw-1)*T/m);
        hdts2(qw,:) = temp2 .* (tt < qw*T/m) .* (tt > (qw-1)*T/m);
        hdts3(qw,:) = temp3 .* (tt < qw*T/m) .* (tt > (qw-1)*T/m);
        hdts4(qw,:) = temp4 .* (tt < qw*T/m) .* (tt > (qw-1)*T/m);
    end
    
    
    %% save
    to_save = false;
    if to_save
       s = strcat('Data/forNRP/hdts_by_', num2str(compression), '.mat');
       save(s, 'hdts1', 'hdts2', 'hdts3', 'hdts4');
    end
    
    %% plot
    do_plot = false;
    if do_plot
       figure;
       plot(tt, temp1, 'LineWidth', 2, 'Color', 'k')
       figure;
       plot(tt, temp2, 'LineWidth', 2, 'Color', 'k')
       ylim([0, 1])
       figure;
       plot(tt, temp3, 'LineWidth', 2, 'Color', 'k')
       figure;
       plot(tt, temp4, 'LineWidth', 2, 'Color', 'k')
    end
    
end