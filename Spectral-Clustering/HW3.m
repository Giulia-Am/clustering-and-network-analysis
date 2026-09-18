%% Spectral Clustering Analysis on Les Misérables Network
% The network represents character co-occurrences in Victor Hugo's novel
% "Les Misérables". Edge weights indicate the number of shared scenes.
% The goal is to identify communities using spectral clustering techniques.

clear; clc; close all;

%% --- 1. Read GML File ---
filename = 'lesmis.gml';
fid = fopen(filename,'r');
if fid == -1
    error('File not found.');
end

nodes = {};
edges = [];
weights = [];

line = fgetl(fid);
while ischar(line)
    line = strtrim(line);

    if startsWith(line,'node')
        node_id = -1;
        while true
            line = strtrim(fgetl(fid));
            if contains(line,'id')
                node_id = str2double(regexp(line,'\d+','match','once'));
            elseif contains(line,'label')
                tokens = regexp(line,'"([^"]+)"','tokens');
                node_label = tokens{1}{1};
            elseif strcmp(line,']')
                break;
            end
        end
        nodes{node_id+1} = node_label;
    end

    if startsWith(line,'edge')
        source_id = -1;
        target_id = -1;
        value = 1;
        while true
            line = strtrim(fgetl(fid));
            if contains(line,'source')
                source_id = str2double(regexp(line,'\d+','match','once'));
            elseif contains(line,'target')
                target_id = str2double(regexp(line,'\d+','match','once'));
            elseif contains(line,'value')
                value = str2double(regexp(line,'\d+','match','once'));
            elseif strcmp(line,']')
                break;
            end
        end
        edges = [edges; source_id+1, target_id+1];
        weights = [weights; value];
    end

    line = fgetl(fid);
end
fclose(fid);

%% --- 2. Adjacency Matrix and Graph Laplacians ---
N = length(nodes);
Adj = zeros(N);

for i = 1:size(edges,1)
    Adj(edges(i,1), edges(i,2)) = weights(i);
    Adj(edges(i,2), edges(i,1)) = weights(i);
end

D = diag(sum(Adj,2));
L_un = D - Adj;                          % Unnormalized Laplacian
L_sym = eye(N) - D^(-0.5) * Adj * D^(-0.5);  % Symmetric normalized Laplacian
L_rw  = eye(N) - D^(-1) * Adj;                % Random Walk Laplacian

G = graph(Adj, nodes);

% Visualize graph
figure;
plot(G);
title('Les Misérables Network');

%% --- 3. Spectral Analysis and Eigengap Computation and Validation ---
k_max = 10; 
m_krylov = 25;

% --- Unnormalized Laplacian ---
[V_un, d_un] =  arnoldi_shift_invert(L_un, k_max,m_krylov);
d_un = sort(d_un);
diff_un = diff(d_un);
[~, idx_gap_un] = max(diff_un(2:end));
k_opt_un = idx_gap_un + 1;

 % validation MATLab eigs
[~, D_ref_un] = eigs(L_un, k_max, 'sa');      % MATLAB reference
d_ref_un = sort(diag(D_ref_un));
abs_error_un = abs(d_ref_un - d_un);

T_eig_un = table((1:k_max)', d_ref_un, d_un, abs_error_un, ...
    'VariableNames', {'Index', 'eigs_System', 'Arnoldi_Manual', 'Absolute_Error'});

fprintf('\nEigenvalue comparison for Unnormalized Laplacian L:\n');
disp(T_eig_un);
fprintf('Optimal k for L (unnormalized): %d\n', k_opt_un);


% --- Symmetric Normalized Laplacian ---
[V_sym, d_sym] = arnoldi_shift_invert(L_sym, k_max, m_krylov);
d_sym = sort(d_sym);
diff_sym = diff(d_sym);
[~, idx_gap_sym] = max(diff_sym(2:end));
k_opt_sym = idx_gap_sym + 1;

% Validation MATLab eigs 
[~, D_ref_sym] = eigs(L_sym, k_max, 'sa');      % MATLAB reference
d_ref_sym = sort(diag(D_ref_sym));
abs_error_sym = abs(d_ref_sym - d_sym);

T_eig_sym = table((1:k_max)', d_ref_sym, d_sym, abs_error_sym, ...
    'VariableNames', {'Index', 'eigs_System', 'Arnoldi_Manual', 'Absolute_Error'});

fprintf('\nEigenvalue comparison for Symmetric Normalized Laplacian L_sym:\n');
disp(T_eig_sym);
fprintf('Optimal k for L_sym (symmetric normalized): %d\n', k_opt_sym);

% --- Random Walk Laplacian ---
[V_rw, d_rw] = arnoldi_shift_invert(L_rw, k_max, m_krylov);
d_rw = sort(d_rw);
diff_rw = diff(d_rw);
[~, idx_gap_rw] = max(diff_rw(2:end));
k_opt_rw = idx_gap_rw + 1;

% Validation MATLAb eigs 
[~, D_ref_rw] = eigs(L_rw, k_max, 'sa');      % MATLAB reference
d_ref_rw = sort(diag(D_ref_rw));
abs_error_rw = abs(d_ref_rw - d_rw);

T_eig_rw = table((1:k_max)', d_ref_rw, d_rw, abs_error_rw, ...
    'VariableNames', {'Index', 'eigs_System', 'Arnoldi_Manual', 'Absolute_Error'});

fprintf('\nEigenvalue comparison for Random Walk Laplacian L_rw:\n');
disp(T_eig_rw);
fprintf('Optimal k for L_rw (random walk): %d\n', k_opt_rw);
fprintf('Optimal k values:\n  L (unnormalized) = %d, L_sym = %d, L_rw = %d\n', k_opt_un, k_opt_sym, k_opt_rw);


%% --- 4. Eigengap Visualization with Highlighted k_opt ---
% --- Unnormalized Laplacian ---
[~, D_ref_un] = eigs(L_un, k_max, 'sa');      % MATLAB reference
d_ref_un = sort(diag(D_ref_un));

figure('Color','w');
plot(1:k_max, d_ref_un, 'bo-', 'LineWidth',1.5, 'MarkerSize',6); hold on;
plot(1:k_max, d_un, 'rs--', 'LineWidth',1.5, 'MarkerSize',6);
plot(k_opt_un, d_un(k_opt_un), 'r*', 'MarkerSize',12,'LineWidth',2);
xlabel('Index'); ylabel('Eigenvalue');
title('Spectrum and Eigengap for L (Unnormalized)');
grid on;
legend('eigs System','Arnoldi Manual','Maximum eigengap (k_{opt})','Location','best');

% --- Symmetric Normalized Laplacian ---
[~, D_ref_sym] = eigs(L_sym, k_max, 'sa');    % MATLAB reference
d_ref_sym = sort(diag(D_ref_sym));

figure('Color','w');
plot(1:k_max, d_ref_sym, 'bo-', 'LineWidth',1.5, 'MarkerSize',6); hold on;
plot(1:k_max, d_sym, 'rs--', 'LineWidth',1.5, 'MarkerSize',6);
plot(k_opt_sym, d_sym(k_opt_sym), 'r*', 'MarkerSize',12,'LineWidth',2);
xlabel('Index'); ylabel('Eigenvalue');
title('Spectrum and Eigengap for L_{sym} (Symmetric Normalized)');
grid on;
legend('eigs System','Arnoldi Manual','Maximum eigengap (k_{opt})','Location','best');

% --- Random Walk Laplacian ---
[~, D_ref_rw] = eigs(L_rw, k_max, 'sa');      % MATLAB reference
d_ref_rw = sort(diag(D_ref_rw));

figure('Color','w');
plot(1:k_max, d_ref_rw, 'bo-', 'LineWidth',1.5, 'MarkerSize',6); hold on;
plot(1:k_max, d_rw, 'rs--', 'LineWidth',1.5, 'MarkerSize',6);
plot(k_opt_rw, d_rw(k_opt_rw), 'r*', 'MarkerSize',12,'LineWidth',2);
xlabel('Index'); ylabel('Eigenvalue');
title('Spectrum and Eigengap for L_{rw} (Random Walk)');
grid on;
legend('eigs System','Arnoldi Manual','Maximum eigengap (k_{opt})','Location','best');
% Based on the computed spectra:
% 1. Eigengap Consistency: While the unnormalized Laplacian (L) suggests 
%    a coarse partition (k=3), both normalized formulations (L_sym, L_rw) 
%    reveal a sharper and more significant eigengap at k=9. This indicates 
%    that normalization effectively uncovers finer cluster structures by 
%    mitigating the influence of node degree heterogeneities.
%
% 2. Spectral Properties: The eigenvalues of L_sym and L_rw are identical 
%    and well-bounded, confirming the theoretical consistency of the 
%    normalized approach. This leads to a more stable embedding in the 
%    spectral domain.

% Key observations from the spectral plots:
% 1. Spectral Stability: For the normalized formulations (L_sym and L_rw), 
%    the manual Arnoldi method shows perfect agreement with the system's 
%    reference values across the entire computed spectrum (k=1 to 10).
%
% 2. Numerical Divergence in Unnormalized L: The unnormalized Laplacian (L) 
%    exhibits a visible divergence between Arnoldi and 'eigs' for higher-index 
%    eigenvalues (index > 7). This is theoretically expected due to the 
%    larger spectral radius and poorer conditioning of the unnormalized 
%    Laplacian, which accelerates the accumulation of rounding errors during 
%    the orthogonalization of the Krylov subspace.
%
% 3. Eigengap Reliability: The sharp eigengap at k=9 in the normalized 
%    plots—upheld by both solvers—confirms that normalization not only 
%    improves cluster resolution but also ensures the numerical robustness 
%    of the eigenvector extraction process.

%% --- 5. Automatic Connected Components Check ---
tolerance = 1e-10;  % tolerance for zero eigenvalue

% --- L_un ---
num_zero_un = sum(abs(d_un) < tolerance);
fprintf('\nNumber of connected components (L_un) = %d\n', num_zero_un);

% --- L_sym ---
num_zero_sym = sum(abs(d_sym) < tolerance);
fprintf('Number of connected components (L_sym) = %d\n', num_zero_sym);

% --- L_rw ---
num_zero_rw = sum(abs(d_rw) < tolerance);
fprintf('Number of connected components (L_rw) = %d\n', num_zero_rw);

% Commento:
% For all three Laplacians, the number of eigenvalues close to zero indicates
% the number of connected components in the graph. If num_zero = 1, the graph
% is fully connected, as expected. This also explains why the first eigenvalue
% is unique and the clustering methods must partition nodes within a single 
% connected component, without any naturally disconnected subgraphs.


%% --- 6. Spectral Clustering (Fixed k=3 and Eigengap-based k) ---
k_fixed = 3;
% --- L (unormalized) ---
% Compute spectral embedding for k=3
[V_un_fixed, ~] = eigs(L_un, k_fixed, 'sa');
Y_un_fixed = V_un_fixed ./ sqrt(sum(V_un_fixed.^2, 2) + eps);
% Clustering with built-in eigs + kmeans
labels_un_builtin = kmeans(Y_un_fixed, k_fixed, 'Replicates', 15);

% Clustering with my_spectral_clustering_full
[labels_un_manual_fixed, Y_un_manual_fixed] = spectral_clustering(L_un, k_fixed, m_krylov);

%% --- L_sym --- 
% Eigengap-based optimal L_sym (k=9)
[V_sym_opt, ~] = eigs(L_sym, k_opt_sym, 'sa');
Y_sym_opt = V_sym_opt ./ sqrt(sum(V_sym_opt.^2, 2) + eps);
labels_sym_opt = kmeans(Y_sym_opt, k_opt_sym, 'Replicates', 15);

% My full spectral clustering function 
[labels_sym_manual, Y_sym_manual] = spectral_clustering(L_sym, k_opt_sym, m_krylov);

% Brute-force kmeans for k=3
[V_sym_fixed, ~] = eigs(L_sym, k_fixed, 'sa');
Y_sym_fixed = V_sym_fixed ./ sqrt(sum(V_sym_fixed.^2,2) + eps);
labels_sym_fixed = kmeans(Y_sym_fixed, k_fixed, 'Replicates', 15);


%% --- L_rw ---
% Eigengap-based optimal k
[V_rw_opt, ~] = eigs(L_rw, k_opt_rw, 'sa');
Y_rw_opt = V_rw_opt ./ sqrt(sum(V_rw_opt.^2, 2) + eps);
labels_rw_opt = kmeans(Y_rw_opt, k_opt_rw, 'Replicates', 15);

% My full spectral clustering function for random walk Laplacian
[labels_rw_manual, Y_rw_manual] = spectral_clustering(L_rw, k_opt_rw, m_krylov);

% Brute-force kmeans for k=3
[V_rw_fixed, ~] = eigs(L_rw, k_fixed, 'sa');
Y_rw_fixed = V_rw_fixed ./ sqrt(sum(V_rw_fixed.^2,2) + eps);
labels_rw_fixed = kmeans(Y_rw_fixed, k_fixed, 'Replicates', 15);

%% --- 7. Plot Spectral Clustering Comparison with Fixed Colors --- 

% Define fixed color map
k_fixed = 3;
cluster_colors = lines(k_fixed);  % same colors for all plots

% --- L_un (k=3) ---
labels_un_builtin_aligned = matchClusters(labels_un_manual_fixed, labels_un_builtin);

figure('Color','w');
subplot(1,2,1);
p1 = plot(G,'Layout','force'); hold on;
for c = 1:k_fixed
    highlight(p1, find(labels_un_manual_fixed==c), 'NodeColor', cluster_colors(c,:), 'MarkerSize',6);
end
title('L Unnormalized: Manual Spectral Clustering');

subplot(1,2,2);
p2 = plot(G,'Layout','force'); hold on;
for c = 1:k_fixed
    highlight(p2, find(labels_un_builtin_aligned==c), 'NodeColor', cluster_colors(c,:), 'MarkerSize',6);
end
title('L Unnormalized: Built-in eigs + kmeans');
% 1. VALIDATION: The 'Manual' implementation perfectly matches the
% 'Built-in' software results, confirming the correctness of the
% eigenvector decomposition and k-means siggnment.
% 2. TOPOLOGY: the algorithm effectively identifies three distinct
% communities: 
%   - BLUE: the bishop Myriel social circle 
%   - RED: small peripheral components/outliers 
%   - YELLOW: the dense core og the network containing the main
%   protagonists. 
% --- L_sym (k_opt = 9) ---
k_opt = k_opt_sym;  % use the optimal k
cluster_colors_opt = lines(k_opt);  % fixed color map for optimal clusters

labels_sym_builtin_aligned = matchClusters(labels_sym_manual, labels_sym_opt);

figure('Color','w');
subplot(1,2,1);
p1 = plot(G,'Layout','force'); hold on;
for c = 1:k_opt
    highlight(p1, find(labels_sym_manual==c), 'NodeColor', cluster_colors_opt(c,:), 'MarkerSize',6);
end
title(sprintf('L_{sym} (k=%d): Manual Spectral Clustering', k_opt));

subplot(1,2,2);
p2 = plot(G,'Layout','force'); hold on;
for c = 1:k_opt
    highlight(p2, find(labels_sym_builtin_aligned==c), 'NodeColor', cluster_colors_opt(c,:), 'MarkerSize',6);
end
title(sprintf('L_{sym} (k=%d): Built-in eigs + kmeans', k_opt));
% Spectral Clustering Comparison (k=9)
% The results show that the manual implementation of the symmetric Laplacian (L_sym) 
% is consistent with the built-in solvers. The clustering effectively identifies 
% narrative communities, such as the Bishop's circle (top) and the revolutionary 
% student groups, proving that the algorithm captures the graph's topology correctly.

% --- L_rw (k_opt = 9) ---
k_opt = k_opt_rw;
cluster_colors_opt_rw = lines(k_opt);

labels_rw_builtin_aligned = matchClusters(labels_rw_manual, labels_rw_opt);

figure('Color','w');
subplot(1,2,1);
p1 = plot(G,'Layout','force'); hold on;
for c = 1:k_opt
    highlight(p1, find(labels_rw_manual==c), 'NodeColor', cluster_colors_opt_rw(c,:), 'MarkerSize',6);
end
title(sprintf('L_{rw} (k=%d): Manual Spectral Clustering', k_opt));

subplot(1,2,2);
p2 = plot(G,'Layout','force'); hold on;
for c = 1:k_opt
    highlight(p2, find(labels_rw_builtin_aligned==c), 'NodeColor', cluster_colors_opt_rw(c,:), 'MarkerSize',6);
end
title(sprintf('L_{rw} (k=%d): Built-in eigs + kmeans', k_opt));
% Spectral Clustering Comparison (k=9)
% The results show that the manual implementation of the symmetric Laplacian (L_sym) 
% is consistent with the built-in solvers. The clustering effectively identifies 
% narrative communities, such as the Bishop's circle (top) and the revolutionary 
% student groups, proving that the algorithm captures the graph's topology correctly.

%% --- 8. Comparative Visualization of All Laplacians (k=3, Aligned Clusters) ---
k_fixed = 3;
cluster_colors_fixed = [1 0 0; 0 0.7 0; 0 0 1];  % Red, Green, Blue for clusters 1-3

% Compute manual spectral clustering for k=3 if not done
labels_sym_manual_fixed_aligned = matchClusters(labels_un_manual_fixed, labels_sym_fixed);
labels_rw_manual_fixed_aligned  = matchClusters(labels_un_manual_fixed, labels_rw_fixed);

% Create a figure with 3 subplots side by side for L_un, L_sym, L_rw
figure('Color','w','Name','Comparison of Laplacians with k=3 (Aligned)');

% --- 1. L_un ---
subplot(1,3,1);
p_un = plot(G,'Layout','force'); hold on;
for c = 1:k_fixed
    highlight(p_un, find(labels_un_manual_fixed==c), 'NodeColor', cluster_colors_fixed(c,:), 'MarkerSize',6);
end
title('L_{un} (Unnormalized)');

% --- 2. L_sym ---
subplot(1,3,2);
p_sym = plot(G,'Layout','force'); hold on;
for c = 1:k_fixed
    highlight(p_sym, find(labels_sym_manual_fixed_aligned==c), 'NodeColor', cluster_colors_fixed(c,:), 'MarkerSize',6);
end
title('L_{sym} (Symmetric Normalized)');

% --- 3. L_rw ---
subplot(1,3,3);
p_rw = plot(G,'Layout','force'); hold on;
for c = 1:k_fixed
    highlight(p_rw, find(labels_rw_manual_fixed_aligned==c), 'NodeColor', cluster_colors_fixed(c,:), 'MarkerSize',6);
end
title('L_{rw} (Random Walk)');

% Overall figure title
sgtitle('Comparison of Spectral Clustering on Different Laplacians (k=3, Aligned Clusters)');

% To evaluate the impact of normalization, we compare the different Laplacian
% formulations. The unnormalized Laplacian (L_un) performs poorly, producing
% a highly unbalanced clustering dominated by the central hub (Jean Valjean),
% whose high degree strongly influences the spectral embedding. This behavior
% is expected, as L_un does not account for degree heterogeneity and is known
% to be sensitive to high-degree nodes in non-regular graphs.
%
% In contrast, the normalized Laplacians (L_sym and L_rw) yield much more
% stable and balanced results, effectively partitioning the network into
% three well-defined social and narrative macro-areas. From a theoretical
% standpoint, this consistency is expected since L_sym and L_rw are
% algebraically similar and therefore share the same eigenvalues, leading to
% equivalent spectral embeddings up to a degree-based scaling. Moreover, both
% normalized formulations are directly related to the optimization of the
% Normalized Cut objective, which explains their robustness and the close
% agreement observed in the resulting clusterings.

%The unnormalized Laplacian emphasizes the central role of Jean Valjean and 
% produces a dominant cluster that merges multiple narrative threads. This 
% representation flattens the novels structure and obscures the distinction 
% between different storylines. In contrast, the normalized Laplacians 
% preserve the plurality of perspectives and reflect the novel as a coherent, 
% multi-voiced narrative, clearly separating the main storyline around Valjean 
% and his family, the revolutionary student movement, and the marginal world 
% of crime and social exclusion.

%% --- 9. Spectral Embedding Visualization for L_un, L_sym, L_rw (Fixed k=3) ---

% Fixed number of clusters
k_fixed = 3;

% Define consistent colors for clusters: Red, Green, Blue
cluster_colors_fixed = [1 0 0; 0 0.7 0; 0 0 1];

% Align labels for comparison if needed
labels_sym_aligned_fixed = matchClusters(labels_un_manual_fixed, labels_sym_manual);
labels_rw_aligned_fixed  = matchClusters(labels_un_manual_fixed, labels_rw_manual);

% Create figure
figure('Color','w','Name','Spectral Embeddings Comparison');

%% --- 1. L_un (Unnormalized Laplacian) ---
subplot(1,3,1);
hold on;
for c = 1:k_fixed
    idx = labels_un_manual_fixed == c;
    scatter(Y_un_manual_fixed(idx,2), Y_un_manual_fixed(idx,3), 50, ...
        'MarkerFaceColor', cluster_colors_fixed(c,:), 'MarkerEdgeColor', 'k');
end
grid on;
xlabel('Dimension 2'); ylabel('Dimension 3');
title(sprintf('Spectral Embedding L_{un} (k=%d)', k_fixed));

%% --- 2. L_sym (Symmetric Normalized Laplacian) ---
subplot(1,3,2);
hold on;
for c = 1:k_fixed
    idx = labels_sym_aligned_fixed == c;
    scatter(Y_sym_manual(idx,2), Y_sym_manual(idx,3), 50, ...
        'MarkerFaceColor', cluster_colors_fixed(c,:), 'MarkerEdgeColor', 'k');
end
grid on;
xlabel('Dimension 2'); ylabel('Dimension 3');
title(sprintf('Spectral Embedding L_{sym} (k=%d)', k_fixed));

%% --- 3. L_rw (Random Walk Laplacian) ---
subplot(1,3,3);
hold on;
for c = 1:k_fixed
    idx = labels_rw_aligned_fixed == c;
    scatter(Y_rw_manual(idx,2), Y_rw_manual(idx,3), 50, ...
        'MarkerFaceColor', cluster_colors_fixed(c,:), 'MarkerEdgeColor', 'k');
end
grid on;
xlabel('Dimension 2'); ylabel('Dimension 3');
title(sprintf('Spectral Embedding L_{rw} (k=%d)', k_fixed));

% Overall figure title
sgtitle('Comparison of Spectral Embeddings for L_{un}, L_{sym} & L_{rw} (k=3, Aligned Clusters)');

% Although L_sym and L_rw share the same eigenvalue spectrum, their embeddings 
% differ due to the relationship between their eigenvectors: u_sym = D^(1/2) * u_rw.
% 
% Key observation: The scaling factor (sqrt of node degree 'd') is not global.
% 1. Hub Nodes (e.g., Valjean): High degree 'd' causes a significant displacement 
%    in the L_sym space compared to L_rw.
% 2. Peripheral Nodes: Low degree 'd' results in minimal scaling, keeping them 
%    closer to their relative L_rw positions.
%
% Conclusion: This node-dependent scaling explains why the plots are not simple 
% translations or global scalings of each other. L_rw preserves the transition 
% probability geometry (Random Walk), while L_sym distorts the embedding 
% based on local connectivity density, affecting the visual "spread" of clusters.

%% ---10. Post-Clustering Validation with Silhouette Score (k = 1:10) ---

k_range = 1:10;  % range of k to test
avg_sil_un = zeros(length(k_range),1);
avg_sil_sym = zeros(length(k_range),1);
avg_sil_rw = zeros(length(k_range),1);

for i = 1:length(k_range)
    k_test = k_range(i);
    
    if k_test == 1
        % silhouette non definito per 1 cluster, mettiamo NaN
        avg_sil_un(i) = NaN;
        avg_sil_sym(i) = NaN;
        avg_sil_rw(i) = NaN;
        continue;
    end
    
    % --- L_un ---
    [labels_un_k, ~] = kmeans(Y_un_manual_fixed(:,2:3), k_test, 'Replicates',5);
    s_un = silhouette(Y_un_manual_fixed(:,2:3), labels_un_k);
    avg_sil_un(i) = mean(s_un);
    
    % --- L_sym ---
    [labels_sym_k, ~] = kmeans(Y_sym_manual(:,2:3), k_test, 'Replicates',5);
    s_sym = silhouette(Y_sym_manual(:,2:3), labels_sym_k);
    avg_sil_sym(i) = mean(s_sym);
    
    % --- L_rw ---
    [labels_rw_k, ~] = kmeans(Y_rw_manual(:,2:3), k_test, 'Replicates',5);
    s_rw = silhouette(Y_rw_manual(:,2:3), labels_rw_k);
    avg_sil_rw(i) = mean(s_rw);
end

% --- Plot silhouette score vs k ---
figure('Color','w');
plot(k_range, avg_sil_un, '-o','LineWidth',2); hold on;
plot(k_range, avg_sil_sym, '-s','LineWidth',2);
plot(k_range, avg_sil_rw, '-^','LineWidth',2);
xlabel('Number of Clusters k'); ylabel('Average Silhouette Score');
legend('L_{un}','L_{sym}','L_{rw}','Location','best');
grid on;
title('Silhouette Score Validation for Different Laplacians (k = 1:10)');

% --- Optimal k for each Laplacian ---
[~, idx_un] = max(avg_sil_un);
[~, idx_sym] = max(avg_sil_sym);
[~, idx_rw] = max(avg_sil_rw);

fprintf('Optimal k (L_un) = %d, silhouette = %.3f\n', k_range(idx_un), avg_sil_un(idx_un));
fprintf('Optimal k (L_sym) = %d, silhouette = %.3f\n', k_range(idx_sym), avg_sil_sym(idx_sym));
fprintf('Optimal k (L_rw) = %d, silhouette = %.3f\n', k_range(idx_rw), avg_sil_rw(idx_rw));

% Post-Clustering Interpretation:
%
% The difference between the number of clusters suggested by the silhouette
% score and the number chosen in the spectral embedding is likely due to
% the nature of the silhouette metric. The silhouette score aims to 
% maximize intra-cluster cohesion and inter-cluster separation, so clusters
% that are very close in the embedding space tend to be merged. This explains
% why L_sym and L_rw show an optimal k of 6-7, whereas in the spectral
% embedding k=9 was used to capture finer substructures.
%
% For the unnormalized Laplacian (L_un), the silhouette score indicates k=3
% as optimal, matching the embedding observation. This occurs because L_un
% tends to collapse nodes around high-degree hubs, producing a few large 
% clusters. In this case, intra-cluster cohesion appears high and 
% inter-cluster separation is clear, simply because L_un differentiates 
% less at a finer scale.
