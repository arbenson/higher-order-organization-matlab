function score = Purity(labels, clusters)
%PURITY - calculates purity to evaluate clustering
% score=Purity(labels, clusters)  where labels assigns the
% ground truth and clusters is the clustering assignment.
  

assert(length(labels) == length(clusters));
overlap = 0;
u_clusters = unique(clusters);
for i = 1:length(u_clusters)
    k = u_clusters(i);
    % Find best cluster for this label
    assignments = labels(clusters == k);
    overlap = overlap + sum(assignments == mode(assignments));
end

score = overlap / length(labels);
