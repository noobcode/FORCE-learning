function [zx_train, hdts_train, target_lengths] = trajectory_train(m, dt)
    
    %% trajectories
    [zx1, T1] = load_trajectory('trajectory_1.csv', dt, false);
    [zx2, T2] = load_trajectory('trajectory_2.csv', dt, false);
    [zx3, T3] = load_trajectory('trajectory_3.csv', dt, false);
    [zx4, T4] = load_trajectory('trajectory_4.csv', dt, false);
    
    target_lengths = [T1 T2 T3 T4];
    
    %% tt
    tt1 = (1:1:round(T1/dt))*dt;
    tt2 = (1:1:round(T2/dt))*dt;
    tt3 = (1:1:round(T3/dt))*dt;
    tt4 = (1:1:round(T4/dt))*dt;

    %% create pulse trains
    temp1 = abs(sin(m*pi*tt1/T1));
    
    width = (T2/dt)/m; % pulse width
    temp2 = pulstran(tt2, 1:T2/dt/m:T2/dt, @rectpuls, width);
    
    x = sawtooth(2*pi*m*tt3/T3, 0.5);
    temp3 = (x - min(x)) / ( max(x) - min(x) );
    
    %temp4 = gaussian_pulse_train(T4, dt, m);
    x = sawtooth(2*pi*m*tt4/T4);
    temp4 = (x - min(x)) / ( max(x) - min(x) );
    
    %% create HDTS
    for qw = 1:1:m
        hdts1(qw,:) = temp1 .* (tt1 < qw*T1/m) .* (tt1 > (qw-1)*T1/m);
        hdts2(qw,:) = temp2 .* (tt2 < qw*T2/m) .* (tt2 > (qw-1)*T2/m);
        hdts3(qw,:) = temp3 .* (tt3 < qw*T3/m) .* (tt3 > (qw-1)*T3/m);
        hdts4(qw,:) = temp4 .* (tt4 < qw*T4/m) .* (tt4 > (qw-1)*T4/m);
    end
    
    %% create HDTS train
    imin = round(T1/dt);
    
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    % wait
    hdts_train = [hdts_train, zeros(m, imin)];
    % replay
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    
    %% create trajectory train
    [k, ~] =  size(zx1);

    zx_train = zeros(k, imin);
    % train
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zeros(k,imin)];
    % replay
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    
    %% plot
    %{
    [~, t_len] = size(zx_train);
    tt = (1:1:t_len)*dt/1000;
    
    size(tt), size(zx_train)
    figure;
    plot(tt, zx_train);
    figure;
    plot(tt, hdts_train);
    %}
    
end