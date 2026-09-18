function [L, U, P] = manual_palu(A)
% MANUAL_PALU Computes the PA = LU factorization of a matrix A
% using Gaussian elimination with partial pivoting.
%
% INPUT:
%   A - square matrix (n x n)
%
% OUTPUT:
%   L - lower triangular matrix with unit diagonal
%   U - upper triangular matrix
%   P - permutation matrix such that P*A = L*U

% Get matrix size
n = size(A,1);

% Initialize U as A, L as identity, P as identity
U = A;
L = eye(n);
P = eye(n);

% Loop over each column (except the last)
for k = 1:n-1

    % -------------------------
    % Partial pivoting
    % -------------------------
    % Find the index of the maximum absolute value in column k
    % from row k to n
    [~, pivot] = max(abs(U(k:n, k)));
    pivot = pivot + k - 1;

    % If the pivot row is not the current row, swap rows
    if pivot ~= k
        % Swap rows in U
        U([k pivot], :) = U([pivot k], :);

        % Swap rows in P to keep track of permutations
        P([k pivot], :) = P([pivot k], :);

        % Swap the previously computed multipliers in L
        % (only columns 1 to k-1 are affected)
        if k > 1
            L([k pivot], 1:k-1) = L([pivot k], 1:k-1);
        end
    end

    % -------------------------
    % Gaussian elimination
    % -------------------------
    % Eliminate entries below the pivot
    for i = k+1:n
        % Compute the multiplier
        L(i,k) = U(i,k) / U(k,k);

        % Update the i-th row of U
        U(i,:) = U(i,:) - L(i,k) * U(k,:);
    end
end

end
