function zx_train = two_target(target_length_1, target_length_2, dt)
    T1 = target_length_1;
    T2 = target_length_2;
    
    tt1 = (1:1:T1/dt)*dt/1000;
    tt2 = (1:1:T2/dt)*dt/1000;
    
    %% supervisor samples
    f = 1; % Hz, frequency sinusoid
    zx1 = sin(2*pi*f*tt1);
    %zx2 = %t2.^2 - sqrt(tt2);
    zx2 = sqrt(tt2); %sin(2*pi*f*tt1); %sqrt(tt2);
    
    %% try combination
    zx_comb = 0.5 * zx1 + 0.5 * zx2; %[0.5 0.5]
    %% create supervisor train
    imin = round(1000/dt);
    
    % pre-training
    zx_train = zeros(1, imin);
    
    % during training
    zx_train = [zx_train, zx1];
    zx_train = [zx_train, zx2];
    zx_train = [zx_train, zx1];
    zx_train = [zx_train, zx2];
    
    % post-training
    zx_train = [zx_train, zeros(1, round(2000/dt))];
    % replay
    zx_train = [zx_train, zx1]; % zx1 / zx_comb
    zx_train = [zx_train, zx1]; % zx1 / zx_comb
    zx_train = [zx_train, zx1]; % zx1 / zx_comb
    zx_train = [zx_train, zeros(1, round(2000/dt))];
    zx_train = [zx_train, zx2];    % zx2 / -
    zx_train = [zx_train, zx2];    % zx2 / -
    zx_train = [zx_train, zx2];    % zx2 / -
    
    %% reverse replay (!!! overwrite zx_train !!!)
    
    % pre-training
    zx_train = zeros(1, imin);
    % during training
    zx_train = [zx_train, zx2];
    zx_train = [zx_train, zx2];
    zx_train = [zx_train, zx2];
    % post-training - replay
    zx_train = [zx_train, zx2];
    zx_train = [zx_train, flip(zx2)];
    zx_train = [zx_train, flip(zx2)];
    
    %% plot
    dims = size(zx_train);
    tt = (1:1:dims(2))*dt/1000;
    figure;
    plot(tt, zx_train);
    xlabel('Time (s)');
    title('Supervisor');
    
    figure;
    plot(tt2, zx2);
end