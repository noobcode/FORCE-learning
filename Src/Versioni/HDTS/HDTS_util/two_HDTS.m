function hdts_train = two_HDTS(m, length_signal_1, length_signal_2, dt)
    T1 = length_signal_1;
    T2 = length_signal_2;
    
    tt1 = (1:1:length_signal_1/dt)*dt;
    tt2 = (1:1:length_signal_2/dt)*dt;

    %% first HDTS sample
    temp1 = abs(sin(m*pi*tt1/T1));
    for qw = 1:1:m
        hdts(qw,:) = temp1 .* (tt1 < qw*T1/m) ...
                           .* (tt1 > (qw-1)*T1/m);
    end
    %% second HDTS sample
    
    % rectangular pulse
    width = (T2/dt)/m; % pulse width
    temp2 = pulstran(tt2, 1:T2/dt/m:T2/dt, @rectpuls, width);
    
    % triangular pulse
    %x = sawtooth(2*pi*m*tt2/T2, 0.5);
    %temp2 = (x - min(x)) / ( max(x) - min(x) );
    
    for qw = 1:1:m
        hdts2(qw,:) = temp2 .* (tt2 < qw*T2/m) ...
                            .* (tt2 > (qw-1)*T2/m);
    end
    
    %% third HDTS sample
    % triangular pulse
    x = sawtooth(2*pi*m*tt2/T2, 0.5);
    temp3 = (x - min(x)) / ( max(x) - min(x) );
    
    for qw = 1:1:m
        hdts3(qw,:) = temp3 .* (tt2 < qw*T2/m) ...
                            .* (tt2 > (qw-1)*T2/m);
    end
    
    %% fourth HDTS sample
    % different triangular pulse
    x = sawtooth(2*pi*m*tt2/T2);
    temp4 = (x - min(x)) / ( max(x) - min(x) );
    
    for qw = 1:1:m
        hdts4(qw,:) = temp4 .* (tt2 < qw*T2/m) ...
                            .* (tt2 > (qw-1)*T2/m);
    end
    
    %% try swap
    %tmp = hdts;     hdts = hdts2;   hdts2 = tmp;
    
    % try combination
    hdts_comb = 0.5 * hdts + 0.5 * hdts2;
    
    %% construct HDTS train
    imin = round(1000/dt);
    
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts2];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts2];
    
    % post-training
    hdts_train = [hdts_train, zeros(m, 2000/dt)]; % rest for 1 second
    % replay
    hdts_train = [hdts_train, hdts]; % hdts / hdts_comb
    hdts_train = [hdts_train, hdts]; % hdts / hdts_comb
    hdts_train = [hdts_train, hdts]; % hdts / hdts_comb
    hdts_train = [hdts_train, zeros(m, 2000/dt)]; % rest for 1 second
    hdts_train = [hdts_train, hdts2];   % hdts2 / -
    hdts_train = [hdts_train, hdts2];   % hdts2 / -
    hdts_train = [hdts_train, hdts2];   % hdts2 / -
    
    
    %% reverse replay (!!! overwrite hdts_train !!!)
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    % post training
    hdts_train = [hdts_train, flip(hdts,2)];
    hdts_train = [hdts_train, flip(hdts,2)];
    
    %% two target one HDTS (!!! overwrite hdts_train !!!)
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    % post training
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    hdts_train = [hdts_train, hdts];
    
    %% four HDTS interpolation
    
    hdts_12 = 0.5*hdts + 0.5*hdts2;
    hdts_13 = 0.5*hdts + 0.5*hdts3;
    hdts_23 = 0.5*hdts2 + 0.5*hdts3;
    hdts_123 = 0.34*hdts + 0.33*hdts2 + 0.33*hdts3;
    hdts_1234 = 0.25*hdts + 0.25*hdts2 + 0.25*hdts3 + 0.25*hdts4;
    
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts, hdts2, hdts3];
    hdts_train = [hdts_train, hdts, hdts2, hdts3];
    %hdts_train = [hdts_train, hdts4];
    % post training
    hdts_train = [hdts_train, hdts, hdts2, hdts3];
    hdts_train = [hdts_train, zeros(m, imin), hdts_12, hdts_12];
    hdts_train = [hdts_train, zeros(m, imin), hdts_13, hdts_13];
    hdts_train = [hdts_train, zeros(m, imin), hdts_23, hdts_23];
    hdts_train = [hdts_train, zeros(m, imin), hdts_123, hdts_123];
    
    %hdts_train = [hdts_train, hdts_12, hdts_12, hdts_12, hdts_12, hdts_12, hdts_12];
    
    %hdts_train = [hdts_train, hdts_123, hdts_123];
    %hdts_train = [hdts_train, hdts_1234, hdts_1234];
    
    %% plot the two HDTS sample 
    figure;
    plot(tt1, hdts);
    figure;
    plot(tt2, hdts2);
    figure;
    plot(tt2, hdts3);
    figure;
    plot(tt2, hdts4);
    
    figure;
    dims = size(hdts_train);
    tt = (1:1:dims(2))*dt/1000;
    plot(tt, hdts_train)
    
    figure;
    hold on
    for i=1:m
        plot(tt, hdts_train(i,:)+i, 'LineWidth', 2)
    end
    title('HDTS train');
    xlabel('Time (s)')
    yticks([])
    
    figure;
    hdts_comb = 0.5 * hdts + 0.5 * hdts2;
    plot(tt1, hdts_comb)
    
end