function [eigvec_norm, lambda1]=validateWithEigs(M, name)
% VALIDATEWITHEIGS Compute dominant eigenvector with eigs
% Inputs:
%   M    : matrix to analyze (can be adjacency or Google matrix)
%   name : string, description of the matrix (for display)
%
% This function computes:
%   - Dominant eigenvector (PageRank)
%   - Dominant eigenvalue (lambda1)

fprintf('\n--- Validation with eigs: %s ---\n', name);

% Compute dominant eigenvector and eigenvalue
[eigvec, eigval] = eigs(M, 1, 'largestabs');

% Normalize eigenvector (PageRank)
eigvec_norm = eigvec / sum(eigvec);

% Extract dominant eigenvalue
lambda1 = eigval(1,1);

end
