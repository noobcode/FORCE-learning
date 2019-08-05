function [hdts, hdts_post] = generate_HDTS(m, length_signal, compression, dt)
    
    len_compr = round(length_signal/compression); % length of compressed HDTS

    temp1 = abs(sin(m*pi*((1:1:length_signal/dt)*dt)/length_signal));
    temp2 = abs(sin(m*pi*((1:1:(len_compr)/dt)*dt)/(len_compr)));

    for qw = 1:1:m
        hdts(qw,:) = temp1 .* ((1:1:length_signal/dt)*dt < qw*length_signal/m) ...
                           .* ((1:1:length_signal/dt)*dt > (qw-1)*length_signal/m);
                     
        hdts_post(qw,:) = temp2 .* ((1:1:(len_compr)/dt)*dt < qw*(len_compr)/m) ...
                                .* ((1:1:(len_compr)/dt)*dt > (qw-1)*(len_compr)/m);
    end

end