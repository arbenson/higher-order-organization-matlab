function [LCC, lcc_inds, ci, sizes] = LargestConnectedComponent(A)
% LARGESTCONNECTEDCOMPONENT gets the largest connected component of A
% A is assumed to be undirected.
%
% [LCC, lcc_inds, ci, sizes] = LargestConnectedComponent(A) returns
%   LCC: the largest connected component
%   lcc_inds: the indices in A corresponding to the connected component
%   ci: the component indices of the nodes in A
%   sizes: the sizes of the connected components

[ci, sizes] = components(A);
[~, max_ind] = max(sizes);
lcc_inds = find(ci == max_ind);
LCC = A(lcc_inds, lcc_inds);
