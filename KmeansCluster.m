function [X, T] = KmeansCluster(A, k)
% KMEANSCLUSTER performs a spectral clustering of the (weighted) adjacency
% matrix A into k clusters using the algorithm in
% "On Spectral Clustering: Analysis and an algorithm" by Ng et al.
%
%  [X, T] = KmeansCluster(A, k) returns
%    X: the assignment of nodes to the k clusters
%    T: the embedding of the nodes into R^k used by the algorithm

Ln = nlaplacian(A);
eig_opts = {};
eig_opts.tol = 1e-6;
eig_opts.isreal = true;
eig_opts.issym = true;
[U, ~] = eigs(Ln, k, 'sa', eig_opts);
D = spdiag(1 ./ sqrt(sum(abs(U).^2, 2)));
T = D * U;

rep = 200;  % Number of replications for kmeans
X = kmeans(T, k, 'Replicates', rep, 'EmptyAction', 'singleton');
