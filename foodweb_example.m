% Generate partial results from Section S7.1 in
% Higher-order organization of complex networks. 
% Austin R. Benson, David F. Gleich, and Jure Leskovec.
% Science, 353.6295 (2016): 163--166.

load foodweb_data.mat;

W = MotifAdjacency(A, 'M6');
[LCC, lcc_inds, ci, sizes] = LargestConnectedComponent(W);
% Get the other non-trivial component
other_cc = find(ci == find(sizes == 12));
% Partition the largest connected component into three groups
X = KmeansCluster(LCC, 3);

% Extract names of organisms in clusters
cluster0 = labels(other_cc);
cluster1 = labels(lcc_inds(X == 1));
cluster2 = labels(lcc_inds(X == 2));
cluster3 = labels(lcc_inds(X == 3));
