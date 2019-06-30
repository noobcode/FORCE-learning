function [x_filtered] = lowpass(x, dt, tau)
    %% Return RC low-pass filter output samples, given input samples,
    % time interval dt, and time constant tau
    % y[i] := α * x[i] + (1-α) * y[i-1]
    
    dims = size(x);
    x_filtered = zeros(dims);
    alpha = dt / (tau + dt);
    
    x_filtered(1) = alpha * x(1);
    n = dims(2);
    for i=2:n
        x_filtered(i) = x_filtered(i-1) + alpha * (x(i) - x_filtered(i-1));
    end
    
end