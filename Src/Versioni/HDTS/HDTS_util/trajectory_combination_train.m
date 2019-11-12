function [zx_train, hdts_train, target_length] = trajectory_combination_train(m, dt)    
    [zx1, zx2, zx3, zx4, target_length] = load_four_trajectories_same_length(dt);
    T = target_length;
    [hdts1, hdts2, hdts3, hdts4] = load_four_hdts_same_length(m, T, 1, dt);
    
    %% save to file
    save('Data/forNRP/trajectories_and_hdts.mat', ...
          'zx1', 'zx2', 'zx3', 'zx4', ...
          'hdts1', 'hdts2', 'hdts3', 'hdts4', ...
          'target_length', 'dt', 'm');
      
   %% create HDTS train
    hdts_12 = 0.50*hdts1 + 0.50*hdts2;
    hdts_13 = 0.50*hdts1              + 0.50*hdts3;
    hdts_23 =              0.50*hdts2 + 0.50*hdts3;
    hdts_34 =                           0.50*hdts3 + 0.50*hdts4;
    hdts_14 = 0.50*hdts1                           + 0.50*hdts4;
    hdts_24 =              0.50*hdts2              + 0.50*hdts4;
    hdts_123 = 0.34*hdts1 + 0.33*hdts2 + 0.33*hdts3;
    hdts_124 = 0.34*hdts1 + 0.33*hdts2             + 0.33*hdts4;
    hdts_134 = 0.34*hdts1              + 0.33*hdts3 + 0.33*hdts4;
    hdts_1234 = 0.25*hdts1 + 0.25*hdts2 + 0.25*hdts3 + 0.25*hdts4;
    
    imin = round(T/dt);
    % pre-training
    hdts_train = zeros(m, imin);
    % during training
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    % post training
    %hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4];
    hdts_train = [hdts_train, zeros(m, imin), hdts_12, hdts_12, hdts_12];
    %hdts_train = [hdts_train, zeros(m, imin), hdts_13, hdts_13];
    %hdts_train = [hdts_train, zeros(m, imin), hdts_23, hdts_23];
    %hdts_train = [hdts_train, zeros(m, imin), hdts_123, hdts_123];
    
    %% create trajectory train
    zx_12 =    0.50*zx1 + 0.50*zx2;
    zx_13 =    0.50*zx1            + 0.50*zx3;
    zx_23 =               0.50*zx2 + 0.50*zx3;
    zx_34 =                        + 0.50*zx3 + 0.50*zx4;
    zx_14 =    0.50*zx1                       + 0.50*zx4;
    zx_24 =               0.50*zx2            + 0.50*zx4;
    zx_123 =   0.34*zx1 + 0.33*zx2 + 0.33*zx3;
    zx_124 =   0.34*zx1 + 0.33*zx2            + 0.33*zx4;
    zx_134 =   0.34*zx1            + 0.33*zx3 + 0.33*zx4;
    zx_1234 =  0.25*zx1 + 0.25*zx2 + 0.25*zx3 + 0.25*zx4;
    
    [k, ~] = size(zx1);
    
    % pre-training
    zx_train = zeros(k, imin);
    % during training
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zx1, zx2, zx3, zx4];
    % post-training - replay
    %zx_train = [zx_train, zx1, zx2, zx3, zx4];
    zx_train = [zx_train, zeros(k, imin), zx_12, zx_12, zx_12];
    %zx_train = [zx_train, zeros(k, imin), zx_13, zx_13];
    %zx_train = [zx_train, zeros(k, imin), zx_23, zx_23];
    %zx_train = [zx_train, zeros(k, imin), zx_123, zx_123];
    
end