function [v, k, lambda1] = pagerank_power(A, tol, maxiter)
% PAGERANK_POWER Applies the power method to compute the PageRank vector
% and returns the dominant eigenvalue.
%
% Inputs:
%   A        : column-stochastic matrix 
%   tol      : convergence tolerance
%   maxiter  : maximum number of iterations
%
% Outputs:
%   v        : PageRank vector (sum(v)=1)
%   k        : number of iterations performed
%   lambda1  : estimated dominant eigenvalue of A

[n,~] = size(A);

% Initialize vector uniformly
v = ones(n,1) / n;

% Initialize Rayleigh quotient 
lambda_prev=Inf;

for k = 1:maxiter

    v_next = A * v;
    
    % Normalize to be a probability vector
    v_next = v_next / sum(v_next);

    % Estimate dominant eigenvalue using Rayleigh quotient
    lambda1 = (v' * (A * v)) / (v' * v);
    
    % Check convergence (relative change in lambda)
    if abs(lambda1-lambda_prev)<tol*abs(lambda1)
        v=v_next; 
        break; 
    end
   
   % Prepare for next iteration 
    v = v_next;
    lambda_prev=lambda1; 
end

if k == maxiter
    warning('PageRank power method did not converge.');
end
end
