function [BIAS_replay] = bias_replay(N, BIAS_LOW, BIAS_HIGH, pBiasHigh)
    %% Bias vector
    % Create a vector of size N, where 'pBiasHigh' percent of the values
    % are set to BIAS_HIGH and the rest are set to BIAS_LOW.
    
    fixed_seed = 1; % to keep the same identity of HIGH neurons and LOW neurons
    rng(fixed_seed);
    
    BIAS_distr = rand(N,1) < pBiasHigh;
    
    BIAS_replay = BIAS_LOW * ones(N,1);
    BIAS_replay(BIAS_distr) = BIAS_HIGH;
end