clc; clear;
%% PageRank Main Script 
%---------------------------------------------------------
% This script performs PageRank computation across three scenarios: 
% 1. Example from Exercises 1 and 11 (simple web with few pages)
% 2. Example corresponding to Figure 2.2 (textbook small web)
% 3. General case using a real dataset of links loaded from file 


% Each section is clearly separated and follows a logical workflow:
% -Load or define the adjacency structure
% -Build adjacency matrix
% -Convert to stochastic matrix 
% -Apply PageRank via the Power Method 
% -Validate results using eigs
% -Print formatted Page
% 
% 1Rank scores 
% -Visuaize results 

% Functions:
% -makeStochastic.m
% -pagerank_power.m
% -validateWithEigs.m
% -dispFormattedPageRank.m

%% ===============================================================
% SECTION 1 — PAGE RANK FOR EXERCISES 1 and 11
% ===============================================================
% EXERCISE 1: Extend the 4-page web (Figure 1) by adding Page 5 
% linking to Page 3, and create a mutual link to check if Page 3's 
% rank surpasses Page 1, using the adjacency matrix.

% Step 1: Define the adjacency matrix for a 4-page web graph (Figure 1)
n1 = 4;
A = zeros(n1);
A(2,1)=1; A(3,1)=1; A(4,1)=1;     % page 1 --> 2,3,4
A(3,2)=1; A(4,2)=1;               % page 2 --> 3,4
A(1,3)=1;                         % page 3 --> 1
A(1,4)=1; A(3,4)=1;               % page 4 --> 1,3

disp('Original adjacency matrix A (4x4):');
disp(A);

% Step 2: Convert the columns into a stochastic form (including handling
% dangling nodes)
A = makeStochastic(A);
disp('Column-stochastic matrix A:');
disp(A);

% Step 3: Compute PageRank vector using the Power Method
tol = 1e-12;
maxiter = 500;
[rank_vector, iters, lambda1] = pagerank_power(A, tol, maxiter);

disp('PageRank vector (Power Method, d=1):');
dispFormattedPageRank(rank_vector);
fprintf('Number of iterations: %d\n', iters);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1);

% Step 4: Validate the PageRank vector using MATLAB's eigs function
[eigvec1,eig1]=validateWithEigs(A, 'Original web (4 pages)');

disp('PageRank vector (from eigs):');
dispFormattedPageRank(eigvec1);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', eig1);


%Step 5: Add Page 5 linking to Page 3, and update Page 3 linking to Page 5
n2 = 5;
A2 = zeros(n2);
A2(2,1)=1; A2(3,1)=1; A2(4,1)=1;    % page 1 --> 2,3,4
A2(3,2)=1; A2(4,2)=1;               % page 2 --> 3,4
A2(1,3)=1; A2(5,3)=1;               % page 3 --> 1,5
A2(1,4)=1; A2(3,4)=1;               % page 4 --> 1,3
A2(3,5)=1;                          % page 5 --> 3

disp('New adjacency matrix A2 (5x5):');
disp(A2);

% Step 6: Convert the columns into a stochastic form (including handling
% dangling nodes)
A2 = makeStochastic(A2);
disp('Column-stochastic matrix A2:');
disp(A2);

%Step 7: Compute PageRank vector using the Power Method
[rank_vector2, iters2, lambda1_2] = pagerank_power(A2, tol, maxiter);
disp('PageRank vector (Power Method, new web, d=1):');
dispFormattedPageRank(rank_vector2);
fprintf('Number of iterations: %d\n', iters2);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1_2);

% Step 8: Validation with MATLAB eigs
[eigvec2, eig2]=validateWithEigs(A2,'New web (5 pages)');

% Display results
disp('PageRank vector (from eigs):');
dispFormattedPageRank(eigvec2);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', eig2);

% CONCLUSIONS
% The PageRank vector was successfully computed and validated. Initially, Page 1
% had the highest score. After adding Page 5 with mutual links to Page 3, Page 3
% became the most important page, achieving a higher score as intended by its owners.
% This demonstrates how adding strategic links can influence PageRank in a web network.

%% EXERCISE 11: Add Page 5 with mutual links to Page 3 and 
% compute the new PageRank vector from the Google matrix (d=0.85).

% Step 9: Construct the Google matrix with damping factor d=0.85 for the
% 4-page web
d = 0.85;
M1 = d*A + (1-d)/n1 * ones(n1);

disp('Google matrix M1 (d=0.85):');
disp(M1);

%Step 10: Compute PageRank vector using the Power Method
[rank_vector_m1, iters_m1, lambda1_m1] = pagerank_power(M1, tol, maxiter);

disp('PageRank vector (Power Method, M1):');
dispFormattedPageRank(rank_vector_m1);
fprintf('Number of iterations: %d\n', iters_m1);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1_m1);

% Step 11: Validation with MATLab eig
[eigvec_m1,eig_m1]= validateWithEigs(M1, 'Google matrix M1');

disp('PageRank vector (from eigs):');
dispFormattedPageRank(eigvec_m1);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', eig_m1);

% Step 12:Construct the Google matrix with damping factor d=0.85 for 5-page
% web
M2 = d*A2 + (1-d)/n2 * ones(n2);

disp('Google matrix M2 (d=0.85):');
disp(M2);

% Step 13: Compute PageRank vector using the Power Mehod for the 5-page web
[rank_vector_m2, iters_m2, lambda1_m2] = pagerank_power(M2, tol, maxiter);
disp('PageRank vector (Power Method, M2):');
dispFormattedPageRank(rank_vector_m2);
fprintf('Number of iterations: %d\n', iters_m2);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1_m2);

% Step 14: Validationi with MATLab eigs
[eigvec_m2,eig_m2]=validateWithEigs(M2, 'Google matrix M2');

disp('PageRank vector (from eigs):');
dispFormattedPageRank(eigvec_m2);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', eig_m2);

% Step 15: Comparison for page 3 vs page 1
fprintf('\nComparison after adding page 5:\n');
fprintf('Page 1 PageRank: %.6f\n', rank_vector2(1));
fprintf('Page 3 PageRank: %.6f\n', rank_vector2(3));

% CONCLUSIONS
% The PageRank vector was efficiently computed, and the obtained values agree
% with the validation results and those found in Case A. Since there are no closed
% sub-webs, the differences are minimal, indicating that the problem remains well-posed.
% Additionally, using the Google matrix M improves the convergence speed of the algorithm.

%% ===============================================================
% SECTION 2 — PAGE RANK FOR FIGURE 2.2
% ===============================================================
% This section replicates the standard PageRank example 
% illustrated in Figure 2.2

% Step 1: Define adjacency matrix for 4 pages
n3=5;
A3 = zeros(n3);
A3(2,1)=1;               % page 1 --> 2
A3(1,2)=1;               % page 2 --> 1
A3(4,3)=1;               % page 3 --> 4
A3(3,4)=1;               % page 4 --> 3
A3(3,5)=1; A3(4,5)=1;    % page 5 --> 3,4

disp('Fig 2.2 djacency matrix A3 (5x5):');
disp(A3);

% Step 2: Convert the columns into a stochastic form (including handling
% dangling nodes)
A3 = makeStochastic(A3);
disp('Column-stochastic matrix A3:');
disp(A3);

% Step 3: Construct the Google matrix with damping factor d=0.85 for web
% Fig 2.2
d = 0.85;
M3 = d*A3 + (1-d)/n3 * ones(n3);

disp('Google matrix M3 (d=0.85):');
disp(M3);

% Step 4: Compute PageRank vector using the Power Mehod 
[rank_vector_m3, iters_m3, lambda1_m3] = pagerank_power(M3, tol, maxiter);
disp('PageRank vector (Power Method, M3):');
dispFormattedPageRank(rank_vector_m3);
fprintf('Number of iterations: %d\n', iters_m3);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1_m3);

% Step 5: Validate with MATLab eigs 
[eigvec_m3, eig_m3]=validateWithEigs(M3, 'Google matrix M3');

disp('PageRank vector (from eigs):');
dispFormattedPageRank(eigvec_m3);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', eig_m3);

% CONCLUSIONS
% The PageRank values were successfully computed and are consistent with the validation using eigs.
% Pages 4 and 3 have the highest scores, sharing the top rank, while Pages 2 and 1 are ranked second.

%% ===============================================================
% SECTION 3 — GENERAL CASE (DATASET OF LINKS)
% ===============================================================
% In this section we load a dataset describing a large set of web pages.
% The dataset must contain pairs (fromPage, toPage) representing hyperlinks.

% Example format:
% 1 5 
% 1 7
% 2 4
% meaning: page 1 links to pages 5 and 7; page 2 links to page 4.

% Step 1: Load dataset
file_dataset=fopen("hollins.dat");

% Step 2: Read the number of web pages and links from the dataset
n_WebPages=fscanf(file_dataset,"%d",1);
n_links=fscanf(file_dataset,"%d",1);

% Check for excessive links
if(n_links>n_WebPages^2)
    fprintf("Too many links in the dataset");
end

% Load the list of URLs
lista_url = strings(n_WebPages, 1);
for i=1:n_WebPages
    fscanf(file_dataset, "%d", 1);
    lista_url(i) = fscanf(file_dataset, "%s", 1);

end

% Step 3: Define the adjacency matrix
A4=zeros(n_WebPages,n_WebPages);
for k = 1:n_links
    link = fscanf(file_dataset, "%d %d", 2);
    
    if length(link) < 2
    error("Error reading link %d", k);
    end

    A4(link(2), link(1)) = 1;
end

% Step 4: Convert the columns into a stochastic form (including handling
% dangling nodes)
A4 = makeStochastic(A4);

% Step 5: Construct the Google matrix with damping factor d=0.85
M4 = d*A4 + (1-d)/n_WebPages * ones(n_WebPages);

% Step 6: Compute PageRank vector using the Power Method 
[rank_vector_m4, iters_m4, lambda1_m4] = pagerank_power(M4, tol, maxiter);

top=5;

dispFormattedPageRank(rank_vector_m4, top);
fprintf('Number of iterations: %d\n', iters_m4);
fprintf('Dominant eigenvalue (lambda1): %.12f\n', lambda1_m4);

% Step 8: Validate with MATLAB eigs
[eigvec_m4, eig_m4]=validateWithEigs(M4, 'Google matrix M4');

dispFormattedPageRank(eigvec_m4, top);
fprintf('Dominant eigenvalue (lambda1):%.12f\n', eig_m4);

[sorted_scores_m, sorted_indices_m] = sort(rank_vector_m4, 'descend');
fprintf('\nTop %d pages (Google Matrix d=0.85):\n', top);
for i = 1:min(top, length(rank_vector_m4))
    idx = sorted_indices_m(i);
    fprintf('%d. URL: %s\tScore: %.4f\n', i, lista_url(idx), sorted_scores_m(i));
end


% CONCLUSIONS
% The Power Method converged in 134 iterations. The PageRank 
% values were successfully validated using MATLAB's eigs 
% function. Page 2 is identified as the most important page.
% URL: http://www1.hollins.edu/Docs/CompTech/Network/webmail_faq.htm