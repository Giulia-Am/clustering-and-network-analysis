function A_stoch = makeStochastic(A)
%MAKESTOCHASTIC Converts a matrix A into a column-stochastic matrix 
% Input: A- square matrix (n x n)
% Output: A_stoch - matrix with columns summing to 1 


%Ensure the matrix is valid: replace any negative entries with zero
A(A<0)=0; % safety check 
[n,~] = size(A);
A_stoch=zeros(n); 

for j=1:n 
    col_sum=sum(A(:,j)); 
    if col_sum ~= 0
        A_stoch(:,j)=A(:,j)/col_sum; 
    else 
        % Handle dangling nodes: evenly distribute 
        A_stoch(:,j)=1/n; 
    end
end

end