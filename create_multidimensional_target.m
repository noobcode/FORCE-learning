clear all
close all

%% define target
T = 1000;  % Total time in ms
dt = 0.04; % Integration time step in ms 
nt = round(T/dt); % Time steps

target_frequency = 1;
target_phase = 0;
zx = (sin(2*pi*target_frequency*(1:1:nt)*dt/1000 + target_phase));

%% one-hot encode
n_slices = 10; % number of classes
[zx_bin, bin_edges] = discretize(zx, n_slices);
md_zx_one_hot = full(ind2vec(zx_bin)); % one-hot

%% plot discretized signal and bin histogram
figure(1)
hold on
plot((1:1:nt)*dt/1000, zx, "LineWidth", 2)
stairs((1:1:nt)*dt/1000, rescale(zx_bin,-1,1), "LineWidth", 2)
hold off
xlabel("Time (s)")
ylabel('$x(t)$','Interpreter','LaTeX')
legend("Original signal", "Discretized signal")
title("Signal discretization")

figure(2)
histogram(zx_bin)
xlabel("Bin")
ylabel("Frequency")
title("Histogram")

%% bin centers
bin_centers = zeros(n_slices,1);
for i=1:n_slices
    bin_centers(i) = (bin_edges(i) + bin_edges(i+1))/2;
end

%% smooth 1
sigma_smoothing = 0.001;
zx_smoothed = zeros(size(md_zx_one_hot));

for i=1:nt
   zx_smoothed(:,i) = exp(-(md_zx_one_hot(:,i) - bin_centers)/sigma_smoothing);
end

%% smooth 2
zx_smoothed2 = zeros(size(md_zx_one_hot));
for i=1:nt
   x0 = bin_centers(zx_bin(i));
   zx_smoothed2(:,i) = md_zx_one_hot(:,i) *  exp(-(zx(:,i) - x0)/sigma_smoothing);
end

%% smooth 3 - OK
zx_smoothed3 = zeros(size(md_zx_one_hot));
for i=1:nt
   zx_smoothed3(:,i) = exp(-(zx(:,i) - bin_centers).^2/sigma_smoothing);
end

figure
plot((1:1:nt)*dt/1000, zx_smoothed3, "LineWidth", 2)
dim_labels = "dim " + ["1" "2" "3" "4" "5"];
xlabel("Time (s)")
ylabel('$x_i(t)$','Interpreter','LaTeX')
legend(dim_labels, "Location", "best")
title("Smooothed dimensions")

