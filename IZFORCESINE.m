function [AverageFiringRate, AverageError, AverageErrorStage1, AverageErrorStage2, tspike, current] = ...
    IZFORCESINE(varargin)
%{
The function can be called in two ways:
    1) IZFORCESINE
    2) IZFORCESINE(seed, load_weights, training_setting, replay)
    - seed -- seed for random generator
    - load_weights -- true/false, load weights from file or initialize them randomly
    - training_setting -- [0,3], define simulation intervals
    - replay -- [0,2], replay experiment. 0: ignore,
                                          1: first run
                                          2: second run
%} 
    %% Force Method with Izhikevich Network 
    clearvars -except varargin
    close all
    clc
    
    % TODO: move into a separate function
    if isempty(varargin)
        seed = 1;
        load_weights = false;
        training_setting = 0;
        replay = 0;
    else
        seed = varargin{1};
        if length(varargin) > 1
            load_weights = varargin{2};
        end
        if length(varargin) > 2
           training_setting = varargin{3};
        end
        
        if length(varargin) > 3
            replay = varargin{4};
        end
    end

    plot_results = 1; % set to 1 to plot the results
    
    AverageErrorStage1 = -1; AverageErrorStage2 = -1; % empty variables 
    
    %% Main experimental parameters
    dt = 0.04; %0.2, 0.04; % Integration time step in ms (0.04)
    
    past_ms = 0;  % time-shift of the target (in ms) [keep it in 0-500]
    pEta = 1; % probability of a reservoir neuron to receive feedback from the readout.

    bias_OMEGA = 0.0; % 0 for Nicola&Clopath
 
    % input distributed onto reservoir neurons with probability pInput
    constant_BIAS = 0; pBias = 0.1; A_add_BIAS = 100;
    
    % external sinusiod params
    input_external_sin_to_network = 0;  pSin = 0.2;
    
    % target injection params
    input_target_to_network = 0; pInput = 0.25; A_it = 50;

    global_inhibition = 0; % set to 1 to implement a model of global inhibition
    exhaustion = 0;  % set to 1 to implement a model of synaptic fatigue / STD
    selective_feedback = 0; % set to 1 to project each dimension of 'z' to a different neurons cluster
    high_dimensional_temporal_signal = 1; 
    
    % Target discretization and smoothing
    n_slices = 0; % 0: mono-dimensional / otherwise discretize
    smoothing = 0; % 0: don't smooth (use one-hot) / 1: apply smoothing
    sigma_smoothing = 0.1; % standard deviation of gaussian filter
    
    % replay
    BIAS_LOW = 900;     BIAS_HIGH = 1100;   pBiasHigh = 0.1;
    replay = replay; % 0: ignore; 1: first run (training); 2: second run (replay)
    
    %% frequency and phase of target signal
    target_frequency = 1;       target_phase = 0;   target_length = 1000; % in milliseconds
    fprintf("target frequency / phase: %.2f / %.2f\n", target_frequency, target_phase);

    %% load trajectories
    % overwrites target_length
    %[original_zx, ~,~, target_length] = load_trajectory('trajectory_4.csv', dt, false);
    
    %[original_zx, hdts_train, target_length] = trajectory_train(32, dt); % train
    [original_zx, hdts_train, target_length] = trajectory_combination_train(32, dt);
    %%
    load_weights = load_weights; % 0: load weights from file / 1: initialize weights randomly
    training_setting = training_setting; % 0: Nicola&Clopath, time-based training
                          % 1: cycle-based training
                          % 2: max(5 seconds, 20 cycles) training
                          
    %fprintf("training setting: %d\n", training_setting);
    %% simulation timings (start/end training, total time)
    iterations_per_second = round(1000/dt); % number of loop iterations per second
    iterations_per_target_cycle = round(iterations_per_second / target_frequency);
    
    [imin, icrit, icrit_2, nt] = set_simulation_time(training_setting, dt, iterations_per_target_cycle, target_length);
    
    % TODO: put i_reingect_BIAS above, but dt not defined...
    if replay == 2
        i_reinject_BIAS = round(3000/dt);
    end
    
    fprintf("train start / train end / icrit2 / simulation end: %d / %d / %d / %d ms\n", imin*dt , icrit*dt, icrit_2*dt, nt*dt);
      
    %% Izhikevich Parameters
    N = 1000; %2000 % Number of neurons
    a = 0.01; % adaptation reciprocal time constant  (0.01 IZFORCESINE)
    b = -2; % resonance parameter (-2 IZFORCESICE)
    vreset = -65; % reset voltage 
    d = 0;   %adaptation jump current (200 IZFORCESINE) (0 or negative for low-frequency signals)
    C = 250;  %capacitance 
    vr = -60;   %resting membrane
    ff = 2.5;  %k parameter for Izhikevich, gain on v
    vpeak = 30;  % peak voltage
    vt = vr + 40 - (b/ff); %threshold (vr + 40 - (b/ff) IZFORCESINE)
    u = zeros(N,1);  %initialize adaptation 
    tr = 2;  %synaptic rise time 
    td = 20; %decay time 
    p = 0.1; % sparsity
    G = 5*10^3; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
    Q = 5*10^3; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA. Set Q=0 to remove feedback term
    WE2 = 4*10^3;
    
    if high_dimensional_temporal_signal
       a = 0.002;
       b = 0;
       vt = -40;
       Q = 4 * 10^2;
    end
        
    %% Storage variables for synapse integration  
    IPSC = zeros(N,1); % post synaptic current 
    h = zeros(N,1);     r = zeros(N,1);     hr = zeros(N,1);    JD = zeros(N,1);

    %% Membrane potential initialization
    rng(seed)
    v = vr + (vpeak-vr) * rand(N,1); %initial distribution 
    v_ = v; % These are just used for Euler integration, previous time step storage
    
    %% Target signal, discretization, smoothing, dimensionality
    tt = (1:1:target_length/dt)*dt; % x-axis of target signal
    
    %original_zx = sin(2*pi*target_frequency*tt + target_phase);
    %original_zx = sqrt((tt)) - 0.5 * (tt);
    %original_zx = sqrt(tt) - 0.5 * (tt).^2 +  0.2 * (tt).^3;
    %original_zx = (sin(2*pi*target_frequency*(1:1:nt)*dt/1000 + target_phase));
    %original_zx = tt.^2 - sqrt(tt);
    %original_zx = 1./exp(-tt) - tt.^3;
    %si=3; mu=1; original_zx = 1/(si*sqrt(2*pi)) * exp(- 0.5 * ((tt-mu)/si).^2);
    
    % comment / uncomment
    %{
    target_length_1 = 1000; target_length_2 = 1000; % ms, length of target signals
    original_zx = two_target(target_length_1, target_length_2, dt);
    %}
    
    %original_zx = original_zx/(max(max(original_zx)));  original_zx(isnan(original_zx) == 1) = 0;
    
    zx = original_zx;
    %plot(tt, zx)
    %reconstructed_zx = zeros(size(original_zx));
    
    if n_slices > 0
        [zx, bin_edges, bin_centers] = discretize_and_smooth(original_zx, n_slices, smoothing, sigma_smoothing);
    end
 
    dims_zx = size(zx);     k = dims_zx(1); % target signal dimensionality
    
    fprintf("target dimension: %d\n", k);
    %% generate HDTS
    m = 32; % number of upstates of the HDTS
    compression = 2; % post-training speedup of the replay. > 1 accelerate, < 1 decelerate
    
    hdts = generate_HDTS(m, target_length(1), 1, dt);
    
    hdts_post = generate_HDTS(m, target_length(1), 0.5, dt);
    hdts_post2 = generate_HDTS(m, target_length(1), compression, dt);
    
    dims_hdts_post = size(hdts_post);
    dims_hdts_post2 = size(hdts_post2);
    
    %hdts_train = two_HDTS(m, target_length_1, target_length_2 , dt);
    
    %% Reinjected target signal parameters
    input_target_amplitude = A_it; % amplitude of reinjected target signal
    input_target_recipients = (rand(N,1) < pInput);
    input_target_indiv_coeff = 0.5 + 0.5 * rand(N,1);
    input_target_coeff = input_target_amplitude .* input_target_recipients .* input_target_indiv_coeff;
    
    past_it = round(past_ms/dt); % past_ms in terms of iterations
    
    fprintf("past: %d ms\n", past_ms);
    
    %% General Network Activity and Global Inhibition Term parameters 
    y = zeros(nt,1); % general network activity measure
    A_inh = 100; %100 % coefficient of y
    tau_inh = 100/dt; % parameter of spike smoothing, frequency and/or duration of network burst (keep it to 100-150)
    
    %g = 1.6; % network high gain
    g = 1; % Clopath gain
    sigma = g / ( p*sqrt(N) ); % std. deviation of static weight matrix  
    
    fprintf("network gain: %.2f\n", g);
    %% Synaptic fatigue / Short-term Depression
    has_spiked = zeros(N,1); % binary vector that tells if a neuron has spiked
    index = []; % indeces of neurons that spiked at the current iteration

    y_ad = ones(N,1); % activity descriptor of each neuron
    tau_ad = 10/dt; % time constant of activity descriptor
    exhaustion_coeff = 0.1;
    
    %% oscillations by means of external sinusoidal input current
    A_ext_sin = 200; % wave amplitude of external sinusoid
    f_ext_sin = 4; % oscillation frequency (Hz) [theta oscillations 4-10 Hz]
    omega =  2*pi*f_ext_sin; % angular frequency (radians per second)
    omega = omega * (dt/1000);
    phase_ext_sin = 0; % phase -- e.g. 1/2*pi
    
    ext_sinusoid =  A_ext_sin * sin(omega * (1:1:nt) + phase_ext_sin); % external sinuoid
    if input_external_sin_to_network == 1
       ext_sinusoid_distrib = ones(N,1) .* (rand(N,1)<pSin);
       ext_sin = ext_sinusoid .* ext_sinusoid_distrib;
    end
    
    %% Gating variables (1 open, 0 closed)
    % if train_gate=0 then the weight update is done once at every cycle of
    % the target signal, otherwise every 20 milliseconds
    %train_gate = 1; % gate for training the weights. 0: phase training
                                                   % 1: continuous training
    %phase_percentage = 0; % 0 : 0 / 0.25 : 1/2 pi / 0.5 : pi ...
    %phase_training = round(phase_percentage * iteration_per_target_cycle); % at what percentage of the cycle to start training - [0,1]
    
    %fprintf("train gate (0 phase, 1 time): %d\n", train_gate);
    %fprintf("train at phase percentage %.2f\n", phase_percentage);
    
     %% load weights omega, phi and eta or initialize them randomly
    if load_weights
        fprintf('loading old weights...');
        weights_file = load('weights.mat');
        OMEGA = weights_file.OMEGA;
        BPhi = weights_file.BPhi;
        E = weights_file.E;
        E2 = weights_file.E2;
    else
        fprintf('initializing weights...');
        OMEGA = G * (bias_OMEGA + randn(N,N)) .* (rand(N,N) < p) * sigma; % static weight matrix.
        BPhi = zeros(N,k); %initial decoder.  Best to keep it at 0.
        E = (2*rand(N,k)-1)*Q .* (rand(N,k) < pEta); % Weight matrix is OMEGA0 + E*BPhi';
        E2 = (2*rand(N,m)-1) * WE2;
    end
    fprintf("done.\n")
    %% Some more initializations
    z = zeros(k,1);  % initial approximant
    tspike = zeros(20*nt,2);  % store spike times, 
    ns = 0; % count cumulative total number of spikes
    ns_t = 0; % number of spikes at time t
    BIAS = 1000; % Bias current. The Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding co tant currents to neurons  
    
    if constant_BIAS == 1
        BIAS = BIAS + A_add_BIAS * (rand(N,1) < pBias);  % add constant bias to some neurons
    end
    
    if replay == 1
        % first run
        BIAS = bias_replay(N, BIAS_LOW, BIAS_HIGH, pBiasHigh);
    elseif replay == 2
        % second run
        BIAS = BIAS_LOW;
        BIAS_replay = bias_replay(N, BIAS_LOW, BIAS_HIGH, pBiasHigh);
    end
    
    %% RLS parameters
    Pinv = 2*eye(N); %initial correlation matrix, coefficient is the regularization constant as well (2)
    step = 20; % optimize with RLS only every 20 steps  
    current = zeros(k,nt);  %store the approximant 
    RECB = zeros(nt,5); %store the decoders 
    REC = zeros(nt,10); %Store voltage and adaptation variables for plotting
    RECy_ad = zeros(nt,5); % stores the exhaustion state of synapses of a given neuron
    i=1;
    qq = 1; % to index the HDTS
    qq_post = 1; % to index the compressed version of the HTDS
    qq_post2 = 1;
    losses = zeros(nt, 1);
   
    %% SIMULATION
    tic
    ilast = i ;
    %icrit = ilast; %uncomment this, and restart cell if you want to test
    % performance before icrit.  
    for i = ilast:1:nt
        %% Elapsed time
        if mod(i,250/dt) == 0
            fprintf('Elapsed simulation time %dms\n', round((i-1)*dt));
        end
        % measure general network activity
        y(i+1) = y(i) + dt * (- y(i) + ns_t) / tau_inh;
        
        if exhaustion == 1
            y_ad = max(0, y_ad * (1-dt/tau_ad) + (dt / tau_ad) * ones(N,1)...
            - exhaustion_coeff * has_spiked); % activity descriptor of individual neuron 
        end
        
        %% current
        I = IPSC + E*z + BIAS;    % postsynaptic current (PSC)
        
        if selective_feedback == 1
            % overwrite previous value of I
            I = IPSC + selective_projection(E)*z + BIAS; % PSC + selective feedback of network output
        end
        
        if global_inhibition == 1
            I = I - A_inh*y(i+1);     % PSC + Global Inhibition
        end
        
        if input_external_sin_to_network == 1
            I = I + ext_sin(:,i);     % PSC + External Sinusoidal Input
        end
        
        if input_target_to_network == 1
            I = I + zx(i) .* input_target_coeff * (i > 1 & i < icrit_2); % PSC + target injection
        end
        
        if replay == 2 && i > i_reinject_BIAS
            % overwrite previous value of I
            I = IPSC + E*z + BIAS_replay; 
        end
        
        if high_dimensional_temporal_signal == 1
            if training_setting == 5
                % HDTS, compress
                %I = I + E2 * hdts(:,qq) * (i < icrit_2) + E2 * hdts_post(:,qq_post) * (i >= icrit_2);
                %{ 
                I = I + E2 * hdts(:,qq) * (i < icrit_2) ...
                      + E2 * hdts_post(:,qq_post) * (i >= icrit_2 && i < icrit_2 + 2*round(target_length/dt)) ...
                      + E2 * hdts_post2(:,qq_post2) * (i >= icrit_2 + 2*round(target_length/dt));
                %}
                
                % HDTS, trigger replay
                %I = I + E2*hdts(:,qq)*(i < icrit_2 || i > icrit_2 + 2*round(target_length/dt));
                I = I + E2*hdts(:,qq)*(i < icrit || i > icrit_2);
                %if i >= round(6000/dt) && i <= round(8000/dt)  BIAS = 0; else BIAS = 1000; end
            end
            if training_setting == 6 || training_setting == 8
                % HDTS, change activity
                I = I + E2*hdts_train(:,i);
                %if i >= round(5000/dt) && i <= round(7000/dt)   BIAS = 0; else BIAS = 1000; end
                %if i >= round(10000/dt) && i <= round(12000/dt)   BIAS = 0; else BIAS = 1000; end
            end
            
            if training_setting == 7
                % just HDTS to aid learning
                I = I + E2*hdts(:,qq);
            end
        end
        
        %% EULER INTEGRATE
        % the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
        v = v + dt * ((ff .* (v-vr) .* (v-vt) - u + I))/C ; % v(t) = v(t-1) + dt*v'(t-1)
        u = u + dt * (a*(b*(v_-vr) - u));                   % u(t) = u(t-1) + dt*u'(t-1).
        %%
        if (~isempty(index))
            has_spiked(index) = 0; % reset all values to 0 before recomputing 'index'   
        end
        
        index = find(v >= vpeak); % get index of neurons that fired
        has_spiked(index) = 1; % set binary values to 1 for neurons that have spiked
        ns_t = length(index); % number of spikes at time t
                
        if ~isempty(index)
            if exhaustion==1
                OMEGA_temp = OMEGA;
                for ii = index'
                    OMEGA_temp(:,ii) = y_ad(ii) .* OMEGA(:,ii);
                end
                JD = sum(OMEGA_temp(:,index),2); % compute the increase in current due to spiking
            else
                JD = sum(OMEGA(:,index),2); % compute the increase in current due to spiking  
            end
              
            %tspike(ns+1:ns+length(index),:) = [index,0*index + dt*i];  %  store spike times. Takes longer.
            ns = ns + ns_t; 
        end
        %% exponential filters
        if tr == 0
            %synapse for single exponential
            IPSC = IPSC*exp(-dt/td) + JD*(~isempty(index))/(td);
            r = r * exp(-dt/td) + (v>=vpeak)/td;
        else
            %synapse for double exponential
            IPSC = IPSC*exp(-dt/tr) + h*dt;
            h = h*exp(-dt/td) + JD*(~isempty(index))/(tr*td);  %Integrate the current

            r = r*exp(-dt/tr) + hr*dt; 
            hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
        end
        %% Compute approximant and error
        z = BPhi' * r; % approximant
        
        % to compute error with respect to time-shifted target, set past_ms > 0. past_ms = 0 to compute the "normal" error
        %if ( (i-past_it) > 1 )
        %    err = z - zx(: , i-past_it);
        %else
        %    err = zeros(k,1);
        %end

        if qq      >= dims_zx(2)         qq = 1;      end
        if qq_post >= dims_hdts_post(2)  qq_post = 1; end
        if qq_post2 >= dims_hdts_post2(2)  qq_post2 = 1; end
        
        err = z - zx(:, qq);
        qq = qq + 1;    qq_post = qq_post + 1; qq_post2 = qq_post2 + 1;
        
        %% RLS 
        % apply RLS only every 'step' steps, i.e every dt*step milliseconds
        if i > imin && i < icrit
            % double check --> if (~train_gate) * mod(i-imin, iteration_per_target_cycle) == 1 + phase_training || train_gate * mod(i,step)==1
            if mod(i,step)==1    
                cd = Pinv * r;
                BPhi = BPhi - (cd * err');
                Pinv = Pinv - ((cd)*(cd'))/( 1 + (r')*(cd));
            end 
        end

        %% End iteration...reset.         
        % if spike occurred, reset variables
        u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
        v = v + (vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
        v_ = v;  % sets v(t-1) = v for the next iteration of loop
        
        REC(i,:) = [v(1:5)',u(1:5)'];  
        current(:,i) = z;
        RECB(i,:) = BPhi(1:5);
        RECy_ad(i,:) = y_ad(1:5)'; 
        
        % store loss
        if n_slices > 0
            %% uncomment this for misclassification error
            %one_hot_z = round(z/sum(z)); % normalize to 0-1 then round
            %losses(i) = ~isequal(one_hot_z, zx(:,i));
            
            %% uncomment this for target reconstruction
            [argvalue, index_max] = max(z);
            %reconstructed_zx(:,i) = bin_centers(index_max);
            %reconstructed_zx(:,i) = z' * bin_centers;
        else
            % mean squared error
            losses(i) = mean(err.^2);
        end
        
        %% control frequency
        if training_setting == 4 && i >= icrit_2
            if mod(i, round(10000/dt)) == 1
                %a = a + 0.005
                %a = a - 0.00025
                %a = a + 0.005*(i < icrit_2 + round(60000/dt)) - 0.005*(i > icrit_2 + round(60000/dt))
                
                %a = a + a*0.25
                %a = a - a*0.05
                %a = a + a*0.25*(i < icrit_2 + round(60000/dt)) - a*0.25*(i > icrit_2 + round(60000/dt))
                
                %d = d + 20*(i < icrit_2 + round(60000/dt)) - 20*(i > icrit_2 + round(60000/dt))
                %d = d + 20
                %d = d - 5
                
                %d = d + abs(d)*0.25
                %d = d - abs(d)*0.05
                %d = d + d*0.25*(i < icrit_2 + round(60000/dt)) - d*0.25*(i > icrit_2 + round(60000/dt))
                
                %BPhi = BPhi*1.05
            end
           %b = 0.5 * b
           %d = 2 * d
        end
        
    end
    fprintf("Simulation finished\n");
    %% smooth output
    tau_smooth = 1/dt;
    
    for j=1:k
        current(j,:) = filter_output(current(j,:), dt, tau_smooth);
    end
    
    %% filtering the reconstructed target
    tau_rec = 10/dt; % low-pass filter time constant
    if n_slices > 0
        %reconstructed_zx = filter_output(reconstructed_zx, dt, tau_rec);
        %losses = (reconstructed_zx - original_zx).^2;
    end
    %% Save weights
    save('weights.mat', 'OMEGA', 'BPhi', 'E', 'E2');
    
    save('Data/forNRP/parameters.mat', 'N', 'a', 'b', 'vreset', 'd', 'C', 'vr', 'ff', ...
                'vpeak', 'vt', 'p', 'G', 'Q', 'WE2', 'BIAS');
    save('Data/forNRP/weights.mat', 'OMEGA', 'BPhi', 'E', 'E2', 'N', 'k', 'm');
    
    %% Average Firing Rate and Average Error after training
    % only consider the spikes from the moment RLS is turned off till the end of the simulation.
    % multiply by 1000 to convert milliseconds to seconds.
    tspike = tspike(tspike(:,2) ~= 0,:); 
    M = tspike(tspike(:,2) > dt*icrit); 
    AverageFiringRate = 1000*length(M)/(N*(nt-icrit)*dt);

    % Average Error after training
    errors_after_training = losses(icrit:end); % losses at each dt after RLS is turned off
    AverageError = mean(errors_after_training);
    
    % Average Error in stage-1 and stage-2
    if icrit_2 ~= -1
        errors_stage1 = losses(icrit:icrit_2);
        errors_stage2 = losses(icrit_2:end);
        AverageErrorStage1 = mean(errors_stage1);
        AverageErrorStage2 = mean(errors_stage2);
    end
    
    fprintf("Average Firing Rate: %f\n", AverageFiringRate);
    fprintf("Average Error: %f\n", AverageError);
    
    if icrit_2 ~= -1
        fprintf("Average Error Stage-1: %f\n", AverageErrorStage1);
        fprintf("Average Error Stage-2: %f\n", AverageErrorStage2);
    end
    %% Graphs
    if plot_results == 1
        if k == 1
            plot_approximant_vs_target(zx, current, dt, nt, imin, icrit, icrit_2);
            plot_approximant(current, dt, nt, imin, icrit, icrit_2)
        else
            % plots joint trajectories
            plot_multidim_approx_vs_target(zx, current, dt, nt, imin, icrit, icrit_2, target_length)
        end
            
        %plot_loss(losses, dt, nt, imin, icrit);
        %plot_raster(tspike, ns, dt, nt, imin, icrit);
        %plot_general_network_activity(y, dt, nt);
        %plot_target_vs_reconstruction(original_zx, reconstructed_zx, dt, nt, imin, icrit);
        %plot_membrane_pre_and_post_learning(REC, dt, nt, vpeak, vreset, imin, icrit);
        %plot_eigenvalues_pre_and_post_learning(OMEGA, E, BPhi);
        %plot_decoders(RECB, dt, nt);
        %plot_phase_portait_pre_and_post_training(REC, dt, imin, icrit, nt);
        %plot_spectrogram(current, dt)
        
        %{
        figure;
        plot(tt, zx)
        xlabel('Time (s)')
        title('supervisor')
        %}
        
    end  

%{            
    % plot activity descriptor of neuron 1
    figure(6)
    plot(dt*(1:1:nt)/1000, y_ad), 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('$y_{ad}(t)$','Interpreter','LaTeX')
    title('Activity Descriptor of neuron 1')
%}
end

