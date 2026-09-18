function new_labels = matchClusters(labels_ref, labels_to_align)
% MATCHCLUSTERS Aligns labels_to_align to labels_ref based on maximum overlap
%
% INPUTS:
% labels_ref      : reference labels (manual)
% labels_to_align : labels to be reindexed
%
% OUTPUT:
% new_labels      : labels_to_align reordered to match labels_ref

k = max(labels_ref);
new_labels = zeros(size(labels_to_align));

used = false(1,k); % Keep track of already matched clusters

for i = 1:k
    counts = zeros(k,1);
    for j = 1:k
        if ~used(j)
            counts(j) = sum(labels_ref==i & labels_to_align==j);
        else
            counts(j) = -1; % Ignore already matched clusters
        end
    end
    [~, best_j] = max(counts);
    new_labels(labels_to_align==best_j) = i;
    used(best_j) = true;
end
end