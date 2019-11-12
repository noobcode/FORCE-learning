function [AverageFiringRate, AverageError, AverageErrorStage1,  AverageErrorStage2] = HDTS_IZFORCE
%% Force Method with Izhikevich Network 
clear all
close all
clc 

dt = 0.04; %0.04;%0.04; %Integration time step in ms 
N =  1000;  %Number of neurons 

%% High dimensional temporal signals
high_dimensional_temporal_signal = 1;

%% load trajectories

%% EXP 1 simply learn with/without HDTS and replay
%{
training_setting = 7;
m = 32; % number of pulses of the HDTS
compression = 1;
[original_zx, ~,~, target_length] = load_trajectory('trajectory_4.csv', dt, compression, false);
hdts = generate_HDTS(m, target_length, compression, dt, false);
%}

%% EXP 2 learn with HDTS then replay at half then double velocity
%{
training_setting = 9;
m = 32;
compression_1 = 1; compression_2 = 0.5; compression_3 = 2;
[original_zx_1, ~,~, target_length] = load_trajectory('trajectory_4.csv', compression_1, dt, false);
[original_zx_2, ~,~, target_length_2] = load_trajectory('trajectory_4.csv', compression_2, dt, false);
[original_zx_3, ~,~, target_length_3] = load_trajectory('trajectory_4.csv', compression_3, dt, false);

hdts_1 = generate_HDTS(m, target_length, compression_1, dt, false);
hdts_2 = generate_HDTS(m, target_length_2, compression_2, dt, false);
hdts_3 = generate_HDTS(m, target_length_3, compression_3, dt, false);

original_zx = [original_zx_1, original_zx_1, original_zx_1, original_zx_1, ...
               original_zx_2, original_zx_3, original_zx_3];

hdts_train = zeros(m, round(target_length/dt));
hdts_train = [hdts_train, hdts_1, hdts_1, hdts_1, hdts_2, hdts_3, hdts_3];
%}

%% EXP 3 learn with HDTS then trigger the replay
%{
training_setting = 5;
m = 32;
compression = 1;
[original_zx, ~,~, target_length] = load_trajectory('trajectory_4.csv', dt, compression, false);
hdts = generate_HDTS(m, target_length, compression, dt, false);
%}
%% EXP 4 learn with HDTS, replay forward then backward
%{
training_setting = 10;
m = 32;
compression = 1;
[original_zx_1, ~,~, target_length] = load_trajectory('trajectory_4.csv', dt, compression, false);
hdts = generate_HDTS(m, target_length, compression, dt, false);

original_zx = [original_zx_1, original_zx_1, original_zx_1, original_zx_1, flip(original_zx_1,2)];

hdts_train = zeros(m, round(target_length/dt));
hdts_train = [hdts_train, hdts, hdts, hdts, flip(hdts,2)];
%}

%% EXP 5 learn all trajectories with HDTS then replay them

training_setting = 11;
m = 32;
compression = 1;
[zx1, zx2, zx3, zx4, target_length] = load_four_trajectories_same_length(dt);
[hdts1, hdts2, hdts3, hdts4] = load_four_hdts_same_length(m, target_length, compression, dt);

original_zx = zeros(7,round(target_length/dt));
original_zx = [original_zx, zx1, zx2, zx3, zx4, ...
                            zx1, zx2, zx3, zx4, ...
                            zx1, zx2, zx3, zx4, ...
                            zeros(7,round(target_length/dt)), ...
                            zx1, zx2, zx3, zx4];
%{
original_zx = [original_zx, zx1, zx2, zx3, zx4, ...
                            zx1, zx2, zx3, zx4, ...
                            zx1, zx2, zx3, zx4, ...
                            zx1, zx2, zx3, zx4, ...
                            zeros(7,round(target_length/dt)), ...
                            zx1, zx2, zx3, zx4];
%}
                        
hdts_train = zeros(m, round(target_length/dt));
hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4, ...
                          hdts1, hdts2, hdts3, hdts4, ...
                          hdts1, hdts2, hdts3, hdts4, ...
                          zeros(m, round(target_length/dt)), ...
                          hdts1, hdts2, hdts3, hdts4];
%{
hdts_train = [hdts_train, hdts1, hdts2, hdts3, hdts4, ...
                          hdts1, hdts2, hdts3, hdts4, ...
                          hdts1, hdts2, hdts3, hdts4, ...
                          hdts1, hdts2, hdts3, hdts4, ...
                          zeros(m, round(target_length/dt)), ...
                          hdts1, hdts2, hdts3, hdts4];
%}

%%
%[original_zx, hdts_train, target_length] = trajectory_train(m, dt); % train

%[original_zx, hdts_train, target_length] = trajectory_combination_train(m, dt);   


%% generate HDTS
%compression = 2; % post-training speedup of the replay. > 1 accelerate, < 1 decelerate
%[hdts1, hdts2, hdts3, hdts4] = load_four_hdts_same_length(m, target_length, compression, dt)


%hdts = generate_HDTS(m, target_length(1), 1, dt);
    
%hdts_post = generate_HDTS(m, target_length(1), 0.5, dt);
%hdts_post2 = generate_HDTS(m, target_length(1), compression, dt);
    
%dims_hdts_post = size(hdts_post);
%dims_hdts_post2 = size(hdts_post2);
    
%hdts_train = two_HDTS(m, target_length_1, target_length_2 , dt);

%% training settting
%training_setting = 7;
[imin, icrit, icrit_2, nt] = set_simulation_time(training_setting, dt, -1, target_length);

fprintf("imin: %d ms\n", imin*dt);
fprintf("icrit: %d ms\n", icrit*dt);
fprintf("icrit2: %d ms\n", icrit_2*dt);
fprintf("nt: %d ms\n", nt*dt);

%% Izhikevich Parameters
C = 250;  %capacitance
vr = -60;   %resting membrane 
b = 0;  %b = -2;  %resonance parameter 
ff = 2.5;  %k parameter for Izhikevich, gain on v 
vpeak = 30;  % peak voltage
vreset = -65; % reset voltage 
vt = -40;   %vt = vr+40-(b/ff); %threshold  %threshold 
u = zeros(N,1);  %initialize adaptation 
a = 0.002; %a = 0.01; %adaptation reciprocal time constant 
d = 0; % oppure 0 d = 200; %adaptation jump current 
tr = 2;  %synaptic rise time 
td = 20; %decay time 
p = 0.1; %sparsity 
G = 5*10^3; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
Q = 4 * 10^2;   %Q =5*10^3; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA 
WE2 = 4*10^3;

%Storage variables for synapse integration  
IPSC = zeros(N,1); %post synaptic current 
h = zeros(N,1); 
r = zeros(N,1);
hr = zeros(N,1);
JD = zeros(N,1);

%-----Initialization---------------------------------------------
v = vr+(vpeak-vr)*rand(N,1); %initial distribution 
v_ = v; %These are just used for Euler integration, previous time step storage
rng('shuffle')

%% Target signal  COMMENT OUT TEACHER YOU DONT WANT, COMMENT IN TEACHER YOU WANT. 
zx = original_zx;

%%
dims_zx = size(zx);    
k = dims_zx(1);

OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(p*sqrt(N)); %Static weight matrix.  
z = zeros(k,1);  %initial approximant
BPhi = zeros(N,k); %initial decoder.  Best to keep it at 0.  
tspike = zeros(25*nt,2);  %If you want to store spike times, 
ns = 0; %count toal number of spikes
BIAS = 1000; %Bias current, note that the Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding constant currents to neurons 
E = (2*rand(N,k)-1)*Q;  %Weight matrix is OMEGA0 + E*BPhi'; 
E2 = (2*rand(N,m)-1) * WE2;
%% 
 Pinv = eye(N)*2; %initial correlation matrix, coefficient is the regularization constant as well 
 step = 20; %optimize with RLS only every 50 steps 
 current = zeros(k,nt); %store the approximant 
 RECB = zeros(nt,5); %store the decoders 
 REC = zeros(nt,10); %Store voltage and adaptation variables for plotting 
 i=1;
 
 qq = 1; % to index the HDTS
 qq_post = 1; % to index the compressed version of the HTDS
 qq_post2 = 1;
 
losses = zeros(nt,1);

%% SIMULATION
tic
ilast = i ;
%icrit = ilast; %uncomment this, and restart cell if you want to test
% performance before icrit.  
for i = ilast:1:nt; 
    if mod(i,250/dt) == 0
            fprintf('Elapsed time %dms\n', round((i-1)*dt));
    end
%% EULER INTEGRATE
I = IPSC + E*z  + BIAS;  %postsynaptic current 

if high_dimensional_temporal_signal == 1
    if training_setting == 5
        % HDTS, compress
        %I = I + E2 * hdts(:,qq) * (i < icrit_2) + E2 * hdts_post(:,qq_post) * (i >= icrit_2);
        %{ 
        I = I + E2 * hdts(:,qq) * (i < icrit_2) ...
              + E2 * hdts_post(:,qq_post) * (i >= icrit_2 && i < icrit_2 + 2*round(target_length/dt)) ...
              + E2 * hdts_post2(:,qq_post2) * (i >= icrit_2 + 2*round(target_length/dt));
        %}

        %% HDTS, trigger replay
        %I = I + E2*hdts(:,qq)*(i < icrit_2 || i > icrit_2 + 2*round(target_length/dt));
        I = I + E2*hdts(:,qq)*(i < icrit || i >= icrit_2);
        %if i >= round(6000/dt) && i <= round(8000/dt)  BIAS = 0; else BIAS = 1000; end
    end
    if training_setting == 6 || training_setting == 8 || training_setting == 10
        %% EXP4 inverse replay; change activity
        I = I + E2*hdts_train(:,i);
        %if i >= round(5000/dt) && i <= round(7000/dt)   BIAS = 0; else BIAS = 1000; end
        %if i >= round(10000/dt) && i <= round(12000/dt)   BIAS = 0; else BIAS = 1000; end
    end

    if training_setting == 7
        %% EXP1 just HDTS to aid learning
        I = I + E2*hdts(:,qq);
    end
    
    if training_setting == 9 || training_setting == 11
        %% EXP 2 decelerate (compress = 0.5) then accelerate (compress = 2)
        %% EXP 5 learn all then replay
        I = I + E2*hdts_train(:,i);
    end
    
end

v = v + dt*(( ff.*(v-vr).*(v-vt) - u + I))/C ; % v(t) = v(t-1)+dt*v'(t-1)
u = u + dt*(a*(b*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
%% 
index = find(v>=vpeak);
if length(index)>0
JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking  
tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i];  %uncomment this
%if you want to store spike times.  Takes longer.  
ns = ns + length(index); 
end

%synapse for single exponential 
if tr == 0 
    IPSC = IPSC*exp(-dt/td)+   JD*(length(index)>0)/(td);
    r = r *exp(-dt/td) + (v>=vpeak)/td;
else
    
%synapse for double exponential
IPSC = IPSC*exp(-dt/tr) + h*dt;
h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td);  %Integrate the current

r = r*exp(-dt/tr) + hr*dt; 
hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
end


 z = BPhi'*r; %approximant 
 if qq      >= dims_zx(2)         qq = 1;      end
%if qq_post >= dims_hdts_post(2)  qq_post = 1; end
%if qq_post2 >= dims_hdts_post2(2)  qq_post2 = 1; end
        
err = z - zx(:, qq);
qq = qq + 1;    %qq_post = qq_post + 1; qq_post2 = qq_post2 + 1;

%% RLS 
 if mod(i,step)==1 
if i > imin 
 if i < icrit 
   cd = Pinv*r;
   BPhi = BPhi - (cd*err');
   Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
 end 
end 
 end

% mean squared error
losses(i) = mean(err.^2);


%% Store, and plot.  
u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
v = v+(vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
v_ = v;  % sets v(t-1) = v for the next itteration of loop
REC(i,:) = [v(1:5)',u(1:5)'];  
current(:,i) = z;
RECB(i,:)=BPhi(1:5);
   

end
%%
tspike = tspike(tspike(:,2)~=0,:); 
M = tspike(tspike(:,2)>dt*icrit); 
%AverageFiringRate = 1000*length(M)/(N*(T-dt*icrit))
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

%%
%plot_raster(tspike, N, ns, dt, nt, imin, icrit);
%plot_multidim_approx_vs_target(zx, current, dt, nt, imin, icrit, icrit_2, target_length);
end
