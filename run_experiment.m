clear all 
close all 
clc

trials = 10;
FiringRates = zeros(trials, 1);
Errors = zeros(trials, 1);

tic
for i=1:1:trials
    [avgFR, avgE] = IZFORCESINE(i);
    FiringRates(i) = avgFR;
    Errors(i) = avgE;
end

MeanFiringRate = mean(FiringRates);
MeanError = mean(Errors);

stdFiringRate = std(FiringRates);
stdError = std(Errors);

fprintf("Mean Firing Rate: %.4f +/- %.4f\n", MeanFiringRate, stdFiringRate);
fprintf("Mean Error: %.4f +/- %.4f\n", MeanError, stdError);
fprintf("Elapsed time: %d\n", toc);