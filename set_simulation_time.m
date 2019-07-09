function [imin, icrit, icrit_2, nt] = set_simulation_time(training_setting, dt, iteration_per_target_cycle)
    %% define simulation intervals
    % [0, imin] -- time before starting RLS, gets the network to chaotic attractor
    % (imin, icrit) -- training interval
    % [icrit, icrit_2) -- optional interval, only for training_setting == 2
    % [icrit/icrit_2, nt] -- evaluation interval, measure network
    % performance.
    % T -- total time in ms
    % ns -- total number of iterations
    
    if training_setting == 0
        %% Nicola & Clopath.
        % 5 sec pre-training, 5 sec training, 5 sec post-training
        imin = round(5000/dt);
        icrit = round(10000/dt);    icrit_2 = -1; % dummy value
        T = 15000;
        nt = round(T/dt);
    elseif training_setting == 1
        %% Train for certain number of cycles of the target signal.
        % 5 sec pre-training, 20 target cycles training, 10 target cycles post-training
        imin = round(5000/dt);
        icrit = imin + 20 * iteration_per_target_cycle;      icrit_2 = -1;
        nt = icrit + 10 * iteration_per_target_cycle; % evaluate for 10 cycles
    elseif training_setting == 2
        %% Compute error w.r.t past version of the target.
        % 3 sec pre-training, max(_ sec, _ cycles) training, stage1, stage2 
        imin = round(3000/dt);
        icrit = imin + max(round(5000/dt), 20 * iteration_per_target_cycle);
        icrit_2 = icrit + max(round(5000/dt), 10 * iteration_per_target_cycle);
        nt = icrit_2 + max(round(8000/dt), 10 * iteration_per_target_cycle);
    elseif training_setting == 3
        %% only simulate, no training
        imin = inf; % never train, (i > imin) will never be true
        icrit = 1; % so you compute the loss from the beginning to the end.
        icrit_2 = -1;
        T = 10000;
        nt = round(T/dt);
    end 
end