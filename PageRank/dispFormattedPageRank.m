function dispFormattedPageRank(v, k_vals)
    % Se k_vals non specificato, stampa tutto
    if nargin < 2
        k_vals = length(v);
    end

    % Ordina i valori di PageRank in ordine decrescente
    [v_sorted, idx] = sort(v, 'descend');
    
    % Limita il numero di valori da stampare
    k_vals = min(k_vals, length(v));

    fprintf('--- Top %d PageRank Scores (from highest to lowest) ---\n', k_vals);
    for k = 1:k_vals
        fprintf('%d) Page %d  -->  %.4f\n', k, idx(k), v_sorted(k));
    end
end
