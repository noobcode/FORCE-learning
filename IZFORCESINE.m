function [] = IZFORCESINE(load_weights)
    %% Force Method with Izhikevich Network 
    clear all
    close all
    clc 
    
    if nargin < 1
        load_weights = 0;
    end

    T = 15000; %15000; %Total time in ms
    dt = 0.04; %Integration time step in ms 
    nt = round(T/dt); %Time steps
    N =  2000;  %Number of neurons
    
    losses = zeros(nt, 1);
    %% Target signal  COMMENT OUT TEACHER YOU DONT WANT, COMMENT IN TEACHER YOU WANT. 
    %zx = (sin(2*5*pi*(1:1:nt)*dt/1000)); original target
    zx = (sin(2*5*pi*(1:1:nt)*dt/1000));
    
    k = min(size(zx)); %used to get the dimensionality of the approximant correctly.  Typically will be 1 unless you specify a k-dimensional target function.  
    %% Izhikevich Parameters
    a = 0.01; %adaptation reciprocal time constant
    b = -2;  %resonance parameter 
    vreset = -65; % reset voltage 
    d = 200;   %adaptation jump current 
    C = 250;  %capacitance 
    vr = -60;   %resting membrane
    ff = 2.5;  %k parameter for Izhikevich, gain on v
    vpeak = 30;  % peak voltage
    vt = vr + 40 - (b/ff); %threshold
    u = zeros(N,1);  %initialize adaptation 
    tr = 2;  %synaptic rise time 
    td = 20; %decay time 
    p = 0.3; %0.1; %sparsity
    G = 5*10^3; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
    Q = 5*10^3; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA 

    %Storage variables for synapse integration  
    IPSC = zeros(N,1); %post synaptic current 
    h = zeros(N,1); 
    r = zeros(N,1);
    hr = zeros(N,1);
    JD = zeros(N,1);

    %-----Initialization---------------------------------------------
    v = vr + (vpeak-vr) * rand(N,1); %initial distribution 
    v_ = v; %These are just used for Euler integration, previous time step storage
    rng(1)
     %% General Network Activity and Global Inhibition Term parameters 
     y = zeros(nt,1); % general network activity measure
     gamma = 0; %100 % coefficient of y
     tau = 100/dt; % parameter of spike smoothing, frequency and/or duration of network burst (keep it to 100-150)
     g = 1.6; %1.6 % network gain
     sigma = g / sqrt(p * N); % std. deviation of static weight matrix

     %% Synaptic fatigue / Short-term Depression
     has_spiked = zeros(N,1); % binary vector that tells if a neuron has spiked
     index = []; % indeces of neurons that spiked at the current iteration
     
     % activity descriptor of neuron 1
     y_ad = zeros(nt,1);
     tau_ad = 10/dt; % time constant of activity descriptor
     
     %% oscillations by means of external sinusoidal input current
     A = 25; % wave amplitude
     f = 4; % oscillation frequency (Hz) [theta oscillations 4-10 Hz]
     omega =  2*pi*f; % angular frequency (radians per second)
     omega = omega * (dt/1000);
     phase = 0; % phase -- e.g. 1/2*pi
     
     %% Gating variables (1 open, 0 closed)
     feedback_gate = 0; % gate for the feedback current {0,1}, how about [0,1]?
     train_gate = 0; % gate for training the weights
     
     iteration_per_second = round(1000/dt); % number of loop iterations per second
     iteration_per_cycle = round(iteration_per_second / f); % number of iterations to make a wave cycle (2pi)
     phase_training = round(0 * iteration_per_cycle); % at what percentage of the cycle to start training - [0,1]
     
    %%
    % load weights omega, phi and eta or initialize them randomly
    if load_weights == 1
        fprintf('loading old weights\n');
        weights_file = load('weights.mat');
        OMEGA = weights_file.OMEGA;
        BPhi = weights_file.BPhi;
        E = weights_file.E;
    else
        fprintf('initializing weights\n');
        OMEGA = G * (randn(N,N)) .* (rand(N,N) < p) * sigma; % static weight matrix.
        BPhi = zeros(N,k); %initial decoder.  Best to keep it at 0.
        E = (2*rand(N,k)-1)*Q;  %Weight matrix is OMEGA0 + E*BPhi';
    end
    
    z = zeros(k,1);  %initial approximant
    tspike = zeros(5*nt,2);  %If you want to store spike times, 
    ns = 0; %count total number of spikes
    ns_t = 0; % number of spikes at time t
    BIAS = 1000; % Bias current, note that the Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding co tant currents to neurons  
    %% 
     Pinv = eye(N)*2; %initial correlation matrix, coefficient is the regularization constant as well 
     step = 20; %optimize with RLS only every 20 steps 
     imin = round(5000/dt); %time before starting RLS, gets the network to chaotic attractor 
     icrit = round(10000/dt); %end simulation at this time step 
     current = zeros(nt,k);  %store the approximant 
     RECB = zeros(nt,5); %store the decoders 
     REC = zeros(nt,10); %Store voltage and adaptation variables for plotting 
     i=1;

    %% SIMULATION
    tic
    ilast = i ;
    %icrit = ilast; %uncomment this, and restart cell if you want to test
    % performance before icrit.  
    for i = ilast:1:nt
        % measure general network activity
        y(i+1) = y(i) + dt * (- y(i) + ns_t) / tau;
        % activity descriptor of neuron 1
        y_ad(i+1) = y_ad(i) + dt * (- y_ad(i) - has_spiked(1) + 1) / tau_ad;
        
        %% EULER INTEGRATE
        % uncomment the type of current you want to use
        %I = IPSC + E*z + BIAS;                    % postsynaptic current (PSC)
        % I = IPSC + E*z + BIAS - gamma*y(i+1);     % PSC + Global Inhibition
        % I = IPSC + E*z + BIAS + OMEGA*has_spiked; % PSC + Short-term Depression
        I = IPSC + feedback_gate * E*z + BIAS + A * sin(omega * i + phase); % PSC + External Sinusoidal Input
        
        % the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
        v = v + dt * ((ff .* (v-vr) .* (v-vt) - u + I))/C ; % v(t) = v(t-1) + dt*v'(t-1)
        u = u + dt * (a*(b*(v_-vr) - u)); % u(t) = u(t-1) + dt*u'(t-1).
        %%
        % reset all values to 0 before recomputing 'index'
        has_spiked(index) = 0;    
        index = find(v >= vpeak);
        has_spiked(index) = 1; % set binary values to 1 for neurons that have spiked
        ns_t = length(index); % number of spikes at time t
        
        if length(index) > 0
            JD = sum(OMEGA(:,index), 2); %compute the increase in current due to spiking  
            tspike(ns+1:ns+length(index),:) = [index,0*index + dt*i];  %uncomment this
            %if you want to store spike times.  Takes longer.
            ns = ns + ns_t; 
        end

        if tr == 0
            %synapse for single exponential
            IPSC = IPSC*exp(-dt/td) + JD*(length(index)>0)/(td);
            r = r * exp(-dt/td) + (v>=vpeak)/td;
        else
            %synapse for double exponential
            IPSC = IPSC*exp(-dt/tr) + h*dt;
            h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);  %Integrate the current

            r = r*exp(-dt/tr) + hr*dt; 
            hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
        end

        z = BPhi' * r; % approximant 
        err = z - zx(:,i); % error 
        %% RLS 
        % apply RLS only every 'step' steps, i.e every dt*step milliseconds
        if i > imin && i < icrit
            if (~train_gate) * mod(i-imin, iteration_per_cycle) == 1 + phase_training || train_gate * mod(i,step)==1 
                cd = Pinv * r;
                BPhi = BPhi - (cd * err');
                Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
            end 
         end
        %% End iteration...store and reset.
        % if spike occurred, reset variables
        u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
        v = v + (vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
        v_ = v;  % sets v(t-1) = v for the next iteration of loop
        
        REC(i,:) = [v(1:5)',u(1:5)'];  
        current(i,:) = z';
        RECB(i,:) = BPhi(1:5);
        losses(i) = err^2;
        %% plot
        if mod(i,round(100/dt))==1 
            drawnow
            gg = max(1,i - round(3000/dt));  %only plot for last 3 seconds

            % plot approximant and target
            figure(2)
            plot(dt*(1:1:i)/1000, zx(:,1:1:i),'k','LineWidth',2), hold on
            plot(dt*(1:1:i)/1000, current(1:1:i,:),'b--','LineWidth',2), hold off
            xlabel('Time (s)')
            ylabel('$\hat{x}(t)$','Interpreter','LaTeX')
            legend('Target Signal', 'Approximant')
            %xlim([dt*i-3000,dt*i]/1000)

            % plot decoders
            figure(3)
            plot((1:1:i)*dt/1000, RECB(1:1:i,:))
            xlabel('Time (s)')
            ylabel('Decoders $\phi(t)$', 'Interpreter', 'LaTeX')
            title('Decoders')
            
            % plot squared error
            figure(5)
            plot((1:1:i)*dt/1000, losses(1:1:i,:))
            xlabel('Time (s)')
            ylabel('Error $e(t)$', 'Interpreter', 'LaTeX')
            set(gca, 'YScale', 'log')
            title('Error curve')

            % plot population activity
            figure(14)
            plot(tspike(1:ns,2), tspike(1:ns,1),'k.')
            xlabel('Time (ms)')
            ylabel('Neuron Index')
            title('Raster plot')
            %ylim([0,100])
            
            % plot general network activity
            figure(4)
            plot(dt*(1:1:i)/1000, y(1:1:i), 'LineWidth', 2)
            xlabel('Time (s)')
            ylabel('$y(t)$','Interpreter','LaTeX')
            title('General Network Activity')
            
            % plot activity descriptor of neuron 1
            figure(6)
            plot(dt*(1:1:i)/1000, y_ad(1:1:i), 'LineWidth', 2)
            xlabel('Time (s)')
            ylabel('$y_{ad}(t)$','Interpreter','LaTeX')
            title('Activity Descriptor of neuron 1')
        end   
    end
    % end simulation
    %% saving weights
    save('weights.mat', 'OMEGA', 'BPhi', 'E');

    %% Average Firing Rate after training
    % only consider the spikes from the moment RLS is turned off till the end of the simulation.
    % multiply by 1000 to convert milliseconds to seconds.
    tspike = tspike(tspike(:,2)~=0,:); 
    M = tspike(tspike(:,2)>dt*icrit); 
    AverageFiringRate = 1000*length(M)/(N*(T-dt*icrit)) 

    % Average Error after training
    errors_after_training = losses(icrit*dt:end); % losses at each dt after RLS is turned off
    AverageError = mean(errors_after_training)
    %% Plotting neurons before and after learning
    % "normalize" membrane potential,
    % add 'j' to the membrane potential of neuron j so that the curves do not
    % overlap

    % post learning
    figure(30)
    for j = 1:1:5
        plot((1:1:i)*dt/1000, REC(1:1:i,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([T/1000-2,T/1000]) %  from 13 seconds to 15 seconds
    xlabel('Time (s)')
    ylabel('Neuron Index') 
    title('Post Learning')

    % before learning
    figure(31)
    for j = 1:1:5
        plot((1:1:i)*dt/1000, REC(1:1:i,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([0,imin*dt/1000]) % from 0 to 5 seconds
    xlabel('Time (s)')
    ylabel('Neuron Index') 
    title('Pre-Learning')

    %% plot eigenvalues
    figure(40)
    Z = eig(OMEGA + E*BPhi'); % eigenvalues after learning 
    Z2 = eig(OMEGA); % eigenvalues before learning 

    plot(Z2,'r.'), hold on 
    plot(Z,'k.') 
    legend('Pre-Learning','Post-Learning')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title('Eigenvalues')

end

