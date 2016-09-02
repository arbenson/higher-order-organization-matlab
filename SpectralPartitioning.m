function [cluster, condv, condc, order] = SpectralPartitioning(A)
% SPECTRALPARTITIONING performs a spectral partitioning of A and
% returns several relevant quantities.  It assumes that A is undirected
% and connected.
%
%  [order, condv, comm, condc] = SpectralPartitioning(A) returns
%   cluster: vector of nodes in the smaller side of the partition
%   condv: the sweep conductance vector
%   condc: the conductance of the cluster
%   order: the spectral ordering

part_vec = nfiedler(A);
[~, order] = sort(part_vec);

% Compute the conductance values (vectorized)
B = A(order, order);
B_lower = tril(B);
B_sums = full(sum(B, 2));
B_lower_sums = sum(B_lower, 2);
volumes = cumsum(B_sums);
num_cut = cumsum(B_sums - 2 * B_lower_sums);
total_vol = full(sum(sum(A)));
volumes_other = total_vol * ones(length(order), 1) - volumes;
vols = min(volumes, volumes_other);
scores = num_cut ./ vols;
scores = full(scores(1:(end-1)));

[condc, min_ind] = min(scores);
  
% The output cluster is the smaller of the two sides of the partition
n = size(A, 1);
if min_ind <= floor(n / 2)
  cluster = order(1:min_ind);
else
  cluster = order((min_ind+1):n);
end

condv = min(scores, scores(end:-1:1));
condv = condv(1:ceil(length(scores) / 2));
