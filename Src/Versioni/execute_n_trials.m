clear all 
close all 
clc

trials = 5;
FiringRates = zeros(trials, 1);
Errors = zeros(trials, 1);
ErrorsStage1 = zeros(trials, 1);
ErrorsStage2 = zeros(trials, 1);

for i=1:1:trials
    %[avgFR, avgE] = ORIGINAL_IZFORCESINE;
    %[avgFR, avgE] = external_sinusoid_IZFORCESINE;
    %[avgFR, avgE] = global_inhibition_IZFORCESINE;
    [avgFR, avgE, avgES1, avgES2] = HDTS_IZFORCE;
    FiringRates(i) = avgFR;
    Errors(i) = avgE;
    ErrorsStage1(i) = avgES1;
    ErrorsStage2(i) = avgES2;
end

MeanFiringRate = mean(FiringRates);
MeanError = mean(Errors);
MeanErrorS1 = mean(ErrorsStage1);
MeanErrorS2 = mean(ErrorsStage2);

stdFiringRate = std(FiringRates);
stdError = std(Errors);
stdErrorS1 = std(ErrorsStage1);
stdErrorS2 = std(ErrorsStage2);

fprintf("Mean Firing Rate: %.4f +/- %.4f\n", MeanFiringRate, stdFiringRate);
fprintf("Mean Error: %.4f +/- %.4f\n", MeanError, stdError);
fprintf("Mean Error stage 1: %.4f +/- %.4f\n", MeanErrorS1, stdErrorS1);
fprintf("Mean Error stage 2: %.4f +/- %.4f\n", MeanErrorS2, stdErrorS2);