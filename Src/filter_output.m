function [x_filtered] = filter_output(x, dt, tau)
    %% Low-pass filter output signal for target reconstruction
    % apply a low-pass filter to input samples 'x', using integration time
    % step 'dt', and time constant 'tau'
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