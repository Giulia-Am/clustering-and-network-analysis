function [V_ritz, d_ritz] = arnoldi_shift_invert(L, k, m)
% arnoldi_shif_invert computes k smallest eigenvalues and eigenvectors
% of a symmetric Laplacian matrix L using the shift-and-invert Arnoldi method.
%
% Inputs:
%   L - Unnormalized Laplacian matrix (n x n)
%   k - Number of eigenvalues/eigenvectors to compute
%   m - Dimension of the Krylov subspace (m >= k)
%
% Outputs:
%   V_ritz - Approximated eigenvectors (n x k)
%   d_ritz - Approximated eigenvalues (k x 1)
%
% This implementation uses a shift sigma close to zero to improve
% convergence for the smallest eigenvalues.

    n = size(L, 1);                 % Size of the matrix
    Q = zeros(n, m + 1);            % Orthonormal basis of Krylov subspace
    H = zeros(m + 1, m);            % Upper Hessenberg matrix

    % Shift for shift-and-invert; small positive value to avoid singularity
    sigma = 1e-6;
    A_shifted = L - sigma * eye(n);

    % Initial random vector normalized
    q1 = rand(n, 1);
    Q(:, 1) = q1 / norm(q1);

    % Arnoldi iteration to construct Krylov subspace
    for j = 1:m
        % v=(L - sigma*I)^(-1)*Q(:,j)
        % Solve shifted system: (L - sigma*I)v = Q(:,j)
        % Using manual_palu function
        [L, U, P] = manual_palu(A_shifted);
        v = U \ (L \ (P * Q(:,j)));

        % Modified Gram-Schmidt orthogonalization
        for i = 1:j
            H(i, j) = Q(:, i)' * v;
            v = v - H(i, j) * Q(:, i);
        end

        H(j+1, j) = norm(v);
        if H(j+1, j) > 1e-12
            Q(:, j+1) = v / H(j+1, j);
        else
            % Breakdown: Krylov subspace has converged
            break;
        end
    end

    % Truncate Hessenberg to m x m
    Hm = H(1:m, 1:m);

    % Compute eigenvalues and eigenvectors of Hm
    [V_h, D_h] = eig(Hm);

    % Shift-and-invert relation: approximate eigenvalues of L
    theta = diag(D_h);
    lambda_approx = (1 ./ theta) + sigma;

    % Sort eigenvalues in ascending order and select the first k
    [d_sorted, sort_idx] = sort(real(lambda_approx), 'ascend');
    d_ritz = d_sorted(1:k);

    % Compute the Ritz vectors corresponding to the k smallest eigenvalues
    V_ritz = Q(:, 1:m) * V_h(:, sort_idx(1:k));
end
