function [x, lambda] = nfiedler(A, tol)
% NFIEDLER returns the fiedler vector of the normalized laplacian of A.

if nargin < 2
    tol = 1e-12;
end

L = nlaplacian(A);
n = size(A, 1);
[V, lambdas] = eigs(L + speye(n), 2 , 'sa', struct('tol', tol));
[~, eig_order] = sort(diag(lambdas));
ind = eig_order(end);
x = V(:, ind);
x = x ./ sqrt(sum(A, 2));
lambda = lambdas(ind, ind) - 1;
