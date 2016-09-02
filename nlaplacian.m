function L = nlaplacian(A)
% NLAPLACIAN returns the normalized laplacian of A

d = sum(A,2);
d = full(d); 
d(d ~= 0) = 1 ./ sqrt(d(d ~= 0));
[i, j, v] = find(A);
[m,n] = size(A);
L = sparse(i,j,-v.*(d(i).*d(j)), m, n);
L = L + speye(n); % note the negative above
