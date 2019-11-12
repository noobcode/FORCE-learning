function hdts = generate_HDTS(m, length_signal, compression, dt, do_plot)

    T = round(length_signal/compression);
    tt = (1:1:round(T/dt))*dt;
    
    temp = abs(sin(m*pi*(tt/T)));
    %temp = gaussian_pulse_train(len_compr, dt, m); DOES NOT WORK WHEN
                                                    %target_length not divisible by m
    
    %x = sawtooth(2*pi*m*tt/T);
    %temp = (x - min(x)) / ( max(x) - min(x) );
    
    for qw = 1:1:m                     
        hdts(qw,:) = temp .* (tt < qw*T/m) ...
                          .* (tt > (qw-1)*T/m);
    end

    %% plot
    
    if do_plot
        figure;
        hold on
        for i=1:m
            plot(tt/1000, hdts(i,:)+i, 'LineWidth', 2)
        end
        title('HDTS');
        xlabel('Time (s)')
        yticks([])
        ax = gca;
        ax.FontSize = 16; 
        %xlabel('My Label','FontSize',10)
    end
    
end