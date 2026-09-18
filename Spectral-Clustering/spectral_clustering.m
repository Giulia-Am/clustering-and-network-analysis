function [labels, Y_norm] = spectral_clustering(L, k, m)
% SPECTRAL_CLUSTERING_FULL  spectral clustering implementation
%
%   INPUTS:
%   L : (N x N)   Laplacian matrix
%   k : number of clusters
%   m : dimension of the Krylov subspace used in the Arnoldi method
%
%   OUTPUTS:
%   labels : (N x 1) cluster assignment for each node
%   Y_norm : (N x k) normalized spectral embedding

    %% 1. Extraction of Ritz vectors via Arnoldi shift-invert method
    % Compute the k eigenvectors associated with the smallest eigenvalues of L
    [V_ritz, ~] = arnoldi_shift_invert(L, k, m);
    
    %% 2. Row-wise normalization of the spectral embedding
    % Each row is normalized to unit Euclidean norm
    % This projects the nodes onto the unit hypersphere
    Y_norm = V_ritz ./ (sqrt(sum(V_ritz.^2, 2)) + eps);
    
    %% 3. Manual k-means clustering with multiple replicates
    % Multiple runs are used to reduce sensitivity to local minima
    n = size(Y_norm, 1);
    best_total_dist = inf;
    labels = zeros(n, 1);
    
    for rep = 1:20
        
        % Random initialization of centroids (sampled from data points)
        centroids = Y_norm(randperm(n, k), :);
        curr_labels = zeros(n, 1);
        
        for iter = 1:150
            
            % Compute squared Euclidean distances to each centroid
            Dist = zeros(n, k);
            for c = 1:k
                Dist(:, c) = sum((Y_norm - centroids(c, :)).^2, 2);
            end
            
            % Assign each point to the nearest centroid
            [min_dists, curr_labels] = min(Dist, [], 2);
            
            % Update centroids as the mean of assigned points
            new_centroids = zeros(k, size(Y_norm, 2));
            for c = 1:k
                idx_c = (curr_labels == c);
                if any(idx_c)
                    new_centroids(c, :) = mean(Y_norm(idx_c, :), 1);
                end
            end
            
            % Convergence check based on centroid displacement
            if norm(centroids - new_centroids, 'fro') < 1e-9
                break;
            end
            
            centroids = new_centroids;
        end
        
        % Evaluate the quality of the current replicate
        % using the total within-cluster distance
        total_dist = sum(min_dists);
        
        if total_dist < best_total_dist
            best_total_dist = total_dist;
            labels = curr_labels;
        end
    end
end
