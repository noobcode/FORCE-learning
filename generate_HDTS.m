function hdts = generate_HDTS(m, length_signal, compression, dt)
    
    len_compr = round(length_signal/compression); % length of compressed HDTS

    temp = abs(sin(m*pi*((1:1:(len_compr)/dt)*dt)/(len_compr)));

    for qw = 1:1:m                     
        hdts(qw,:) = temp .* ((1:1:(len_compr)/dt)*dt < qw*(len_compr)/m) ...
                           .* ((1:1:(len_compr)/dt)*dt > (qw-1)*(len_compr)/m);
    end

    %% plot
    plot = false;
    
    if plot
        figure;
        hold on
        for i=1:m
            plot((1:1:(len_compr)/dt)*dt/1000, hdts(i,:)+i, 'LineWidth', 2)
        end
        title('HDTS');
        xlabel('Time (s)')
        yticks([])
    end
end