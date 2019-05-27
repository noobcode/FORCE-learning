clear all
close all

T = 5;  % Total time in ms
dt = 0.04; % Integration time step in ms 
nt = round(T/dt); % Time steps

n_slices = 3;

target_frequency = 1;
target_phase = 0;
zx = (sin(2*pi*target_frequency*(1:1:nt)*dt/1000 + target_phase));

[md_zx, bin_edges] = discretize(zx, n_slices);
%md_zx = full(ind2vec(md_zx));

figure(1)
plot((1:1:nt)*dt/1000, zx)
hold on
stairs((1:1:nt)*dt/1000, rescale(md_zx,-1,1))
hold off

figure(2)
histogram(md_zx)

