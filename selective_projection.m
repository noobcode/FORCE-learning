function [E] = selective_projection(E)
    %% Selective feedback of the network approximant 'z' weighed by static matrix 'E' (eta).
    % Inject the network approximant dimension 'z_i' only to the i-th
    % neuron cluster. The neurons clusters are consecutive with respect to
    % the neurons indeces.
    
    dims = size(E);
    N = dims(1); % number of neurons
    k = dims(2); % number of bins

    cluster_size = round(N/k); % numbers of neurons for each bin
    
    for i=1:k
        % columns to zero
        columns_to_zero = setdiff(1:k, i);
        start_row = (i-1)*cluster_size+1;
        end_row = start_row + cluster_size -1;
        if i == k
            end_row = N;
        end
        % rows to zero
        rows_to_zero = start_row:end_row;
        % set to zero corresponding rows and columns
        E(rows_to_zero, columns_to_zero) = 0; 
    end

end