function [AverageFiringRate, AverageError] = STD_IZFORCESINE
%% Force Method with Izhikevich Network 
clear all
close all
clc 

T = 15000; %Total time in ms
dt = 0.1;%0.04; %0.04; %Integration time step in ms 
nt = round(T/dt); %Time steps
N =  2000;  %Number of neurons 
%% Izhikevich Parameters
C = 250;  %capacitance
vr = -60;   %resting membrane 
b = -2;  %resonance parameter 
ff = 2.5;  %k parameter for Izhikevich, gain on v 
vpeak = 30;  % peak voltage
vreset = -65; % reset voltage 
vt = vr+40-(b/ff); %threshold  %threshold 
u = zeros(N,1);  %initialize adaptation 
a = 0.01; %adaptation reciprocal time constant 
d = 200; %adaptation jump current 
tr = 2;  %synaptic rise time 
td = 20; %decay time 
p = 0.1; %sparsity 
G =5*10^3; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
Q =5*10^3; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA 
Irh = 0.25*ff*(vt-vr)^2; 

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
%rng(1)
%% Target signal  COMMENT OUT TEACHER YOU DONT WANT, COMMENT IN TEACHER YOU WANT. 
f_sin = 1;
zx = (sin(2*pi*f_sin*(1:1:nt)*dt/1000));

%%
k = min(size(zx)); %used to get the dimensionality of the approximant correctly.  Typically will be 1 unless you specify a k-dimensional target function.  

OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(p*sqrt(N)); %Static weight matrix.  
z = zeros(k,1);  %initial approximant
BPhi = zeros(N,k); %initial decoder.  Best to keep it at 0.  
tspike = zeros(5*nt,2);  %If you want to store spike times, 
ns = 0; %count toal number of spikes
nt_t = 0;
BIAS = 1000; %Bias current, note that the Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding constant currents to neurons 
E = (2*rand(N,k)-1)*Q;  %Weight matrix is OMEGA0 + E*BPhi'; 
%% 
 Pinv = eye(N)*2; %initial correlation matrix, coefficient is the regularization constant as well 
 step = 20; %optimize with RLS only every 50 steps 
 imin = round(5000/dt); %time before starting RLS, gets the network to chaotic attractor 
 icrit = round(10000/dt); %end simulation at this time step 
 current = zeros(nt,k);  %store the approximant 
 RECB = zeros(nt,5); %store the decoders 
 REC = zeros(nt,10); %Store voltage and adaptation variables for plotting 
 RECy_ad = zeros(nt,5); % stores the exhaustion state of synapses of a given neuron
 i=1;
 
losses = zeros(nt,1);

%% Short term depression / exhaustion / synaptic fatigue
exhaustion = 1;  % set to 1 to implement a model of synaptic fatigue / STD
has_spiked = zeros(N,1); % binary vector that tells if a neuron has spiked
index = []; % indeces of neurons that spiked at the current iteration

y_ad = ones(N,1); % activity descriptor of each neuron
tau_ad = 100/dt; % 10, time constant of activity descriptor
exhaustion_coeff = 0.1; % 0.1 (nella nuova versione questo deve essere alto 
                        % tipo 1000 per influenzare un minimo)

fprintf("tau_ad = %.2f\n", tau_ad*dt);
fprintf("exhaustion_coeff = %.2f\n", exhaustion_coeff);

%% SIMULATION
tic
ilast = i ;
%icrit = ilast; %uncomment this, and restart cell if you want to test
% performance before icrit.  
for i = ilast:1:nt; 
    if mod(i,250/dt) == 0
            fprintf('Elapsed time %dms\n', round((i-1)*dt));
    end   

if exhaustion == 1
    % activity descriptor of neurons
    %y_ad = max(0, y_ad * (1-dt/tau_ad) + (dt / tau_ad) * (ones(N,1) - exhaustion_coeff * has_spiked)); 
    y_ad = max(0, y_ad * (1-dt/tau_ad) + (dt / tau_ad) * ones(N,1) - exhaustion_coeff * has_spiked);  
end

%% EULER INTEGRATE
I = IPSC + E*z  + BIAS;  %postsynaptic current 
v = v + dt*(( ff.*(v-vr).*(v-vt) - u + I))/C ; % v(t) = v(t-1)+dt*v'(t-1)
u = u + dt*(a*(b*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
%% 

if (~isempty(index))
    has_spiked(index) = 0; % reset all values to 0 before recomputing 'index'   
end

index = find(v >= vpeak); % get index of neurons that fired
has_spiked(index) = 1; % set binary values to 1 for neurons that have spiked
ns_t = length(index); % number of spikes at time t

if ~isempty(index)
    if exhaustion==1
        %OMEGA_temp = OMEGA;
        %for ii = index'
        %    OMEGA_temp(:,ii) = y_ad(ii) .* OMEGA(:,ii);
        %end
        %JD = sum(OMEGA_temp(:,index),2); % compute the increase in current due to spiking
        % oppure equivalentemente
        tmp = y_ad(index)' .* ones(N, length(index)) .* OMEGA(:, index); % compute the increase in current due to spiking
        JD = sum(tmp,2);
    else
        JD = sum(OMEGA(:,index),2); % compute the increase in current due to spiking  
    end

    tspike(ns+1:ns+length(index),:) = [index,0*index + dt*i];  %  store spike times. Takes longer.
    ns = ns + ns_t; 
end

%% filters
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
 err = z - zx(:,i); %error 
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
current(i,:) = z';
RECB(i,:)=BPhi(1:5);
RECy_ad(i,:) = y_ad(1:5)';

end
%%
tspike = tspike(tspike(:,2)~=0,:); 
M = tspike(tspike(:,2)>dt*icrit); 
AverageFiringRate = 1000*length(M)/(N*(T-dt*icrit))

errors_after_training = losses(icrit:end); % losses at each dt after RLS is turned off
AverageError = mean(errors_after_training)

%%
plot_raster(tspike, ns, dt, nt, imin, icrit);
%figname = strcat("Images/pre-learning/short_term_depression/",...
%    "raster_tau_", num2str(tau_ad*dt), "_coeff_", num2str(exhaustion_coeff), ".png");
%saveas(gcf,figname);
plot_approximant_vs_target(zx, current, dt, nt, imin, icrit, -1);
%plot_general_network_activity(y, dt, nt);

% plot activity descriptor of neuron 1
figure;
plot(dt*(1:1:nt)/1000, RECy_ad, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('$y_{ad}(t)$','Interpreter','LaTeX')
title('Activity Descriptors Sample')

%figname = strcat("Images/pre-learning/short_term_depression/",...
%    "descriptor_tau_", num2str(tau_ad*dt), "_coeff_", num2str(exhaustion_coeff), ".png");
%saveas(gcf,figname);
end
