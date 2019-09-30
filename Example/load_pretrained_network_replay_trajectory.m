clear all
close all
clc

%% load trajectories and hdts
trajectories_and_hdts = load('../Data/forNRP/trajectories_and_hdts.mat');

zx1 = trajectories_and_hdts.zx1;
zx2 = trajectories_and_hdts.zx2;
zx3 = trajectories_and_hdts.zx3;
zx4 = trajectories_and_hdts.zx4;

hdts1 = trajectories_and_hdts.hdts1;
hdts2 = trajectories_and_hdts.hdts2;
hdts3 = trajectories_and_hdts.hdts3;
hdts4 = trajectories_and_hdts.hdts4;

dt = trajectories_and_hdts.dt;
target_length = trajectories_and_hdts.target_length;
%%
zx = zx1;
hdts = hdts1;

[n_joints, n_samples] = size(zx);
k = n_joints;

%%
nt = n_samples * 6;

imin = n_samples; % start HDTS
icrit = imin + 3 * n_samples; % remove HDTS

%% Load Izhikevich Parameters
params = load('../Data/forNRP/parameters.mat');

N = params.N; %2000 % Number of neurons
a = params.a; % adaptation reciprocal time constant  (0.01 IZFORCESINE)
b = params.b; % resonance parameter (-2 IZFORCESICE)
vreset = params.vreset; % reset voltage 
d = params.d;   %adaptation jump current (200 IZFORCESINE) (0 or negative for low-frequency signals)
C = params.C;  %capacitance 
vr = params.vr;   %resting membrane
ff = params.ff;  %k parameter for Izhikevich, gain on v
vpeak = params.vpeak;  % peak voltage
vt = params.vt; %threshold (vr + 40 - (b/ff) IZFORCESINE)
u = zeros(N,1);  %initialize adaptation 
tr = 2;  %synaptic rise time 
td = 20; %decay time 
p = params.p; % sparsity
G = params.G; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
Q = params.Q; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA. Set Q=0 to remove feedback term
WE2 = params.WE2;

%% load weights
weights_file = load('../Data/forNRP/weights_4.mat');
OMEGA = weights_file.OMEGA;
BPhi = weights_file.BPhi;
E = weights_file.E;
E2 = weights_file.E2;

%% Storage variables for synapse integration  
IPSC = zeros(N,1); % post synaptic current 
h = zeros(N,1);     r = zeros(N,1);     hr = zeros(N,1);    JD = zeros(N,1);

%%
rng('shuffle')
v = vr + (vpeak-vr) * rand(N,1); %initial distribution 
v_ = v; % These are just used for Euler integration, previous time step storage

%%
z = zeros(k,1);  % initial approximant
BIAS = 1000; % Bias current. The Rheobase is around 950 

%% RLS parameters
Pinv = 2*eye(N); % initial correlation matrix, coefficient is the regularization constant as well (2)
step = 20; % optimize with RLS only every 20 steps  
current = zeros(k,nt);  % store the approximant 
RECB = zeros(nt,5); % store the decoders 
REC = zeros(nt,10); % store voltage and adaptation variables for plotting

qq = 1; % to index the HDTS
losses = zeros(nt, 1);

%% simulation

for i = 1:1:nt
    %% Elapsed time
    if mod(i,250/dt) == 0
        fprintf('Elapsed simulation time %dms\n', round((i-1)*dt));
    end

    %% current
    I = IPSC + E*z + BIAS + E2*hdts(:,qq) .* (i >= imin && i <= icrit);

    %% EULER INTEGRATE
    % the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
    v = v + dt * ((ff .* (v-vr) .* (v-vt) - u + I))/C ; % v(t) = v(t-1) + dt*v'(t-1)
    u = u + dt * (a*(b*(v_-vr) - u));                   % u(t) = u(t-1) + dt*u'(t-1).
    %%
    index = find(v >= vpeak); % get index of neurons that fired

    if ~isempty(index)
        JD = sum(OMEGA(:,index),2); % compute the increase in current due to spiking  
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

    if qq >= n_samples        
        qq = 1;      
    end
   
    err = z - zx(:, qq);
    qq = qq + 1;

    %% End iteration...reset.         
    % if spike occurred, reset variables
    u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    v = v + (vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    v_ = v;  % sets v(t-1) = v for the next iteration of loop
    current(:,i) = z;

    losses(i) = mean(err.^2);

end

plot_multidim_approx_vs_target(zx, current, dt, nt, -1, -1, -1, target_length)

fprintf("Simulation finished\n");