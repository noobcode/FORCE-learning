function [zx_new, bin_edges, bin_centers] = discretize_and_smooth(zx, n_slices, smoothing, sigma_smoothing)
    %% input target signal discretization and smoothing
    % discretize target signal 'zx' in 'n_slices' classes.
    % if 'smoothing' = 1, then apply gaussian smoothing using standard deviation 'sigma_smoothing',
    % if 'smoothing' = 0, then use the one-hot encoded target.
    
    [zx_bin, bin_edges] = discretize(zx, n_slices); % discretize target signal
 
    if smoothing
        %% smooth multi-dimensional target
        bin_centers = zeros(n_slices,1);
        
        % compute bin centers
        for i=1:n_slices
            bin_centers(i) = (bin_edges(i) + bin_edges(i+1))/2;
        end
        
        % smooth the one-hot encoded target
        zx_new = exp(-(zx - bin_centers).^2/sigma_smoothing);
    else
        %% one-hot multi-dimensional target
        zx_new = full(ind2vec(zx_bin)); % one-hot encoding
    end

end