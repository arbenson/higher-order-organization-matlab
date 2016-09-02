function W = MotifAdjacency(A, motif)
% MOTIFADJACENCY forms the motif adjacency matrix for the adjacency
% matrix A and the specified motif.
% 'motif' is one of:
% M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M13
% bifan
% edge
% 
% See http://snap.stanford.edu/higher-order/code.html for
% the naming conventions.
%  
%  W = MotifAdjacency(A, motif) returns
%    W: the motif adjacency matrix

% Ignore diagonals and weights
A = A - diag(diag(A));
A = min(A, 1);

lmotif = lower(motif);
if     strcmp(lmotif, 'm1')
  W = M1(A);
elseif strcmp(lmotif, 'm2')
  W = M2(A);
elseif strcmp(lmotif, 'm3')
  W = M3(A);
elseif strcmp(lmotif, 'm4')
  W = M4(A);
elseif strcmp(lmotif, 'm5')
  W = M5(A);
elseif strcmp(lmotif, 'm6')
  W = M6(A);
elseif strcmp(lmotif, 'm7')
  W = M7(A);
elseif strcmp(lmotif, 'm8')
  W = M8(A);
elseif strcmp(lmotif, 'm9')
  W = M9(A);
elseif strcmp(lmotif, 'm10')
  W = M10(A);
elseif strcmp(lmotif, 'm11')
  W = M11(A);
elseif strcmp(lmotif, 'm12')
  W = M12(A);
elseif strcmp(lmotif, 'm13')
  W = M13(A);
elseif strcmp(lmotif, 'bifan')
  W = Bifan(A);
elseif strcmp(lmotif, 'edge')
  [~, ~, W] = DirectionalBreakup(A);  
else
  error('Error.\nUnknown motif %s', motif);
end
end  

function W = M1(A)
  [~, U, ~] = DirectionalBreakup(A);
  W = (U * U) .* U;
end

function W = M2(A)
  [B, U, ~] = DirectionalBreakup(A);
  C = (B * U) .* U' + (U * B) .* U' + (U * U) .* B;
  W = C + C';
end

function W = M3(A)
  [B, U, ~] = DirectionalBreakup(A);
  C = (B * B) .* U + (B * U) .* B + (U * B) .* B;
  W = C+ C';
end

function W = M4(A)
  [B, ~, ~] = DirectionalBreakup(A);
  W = (B * B) .* B;
end

function W = M5(A)
  [~, U, ~] = DirectionalBreakup(A);
  T1 = (U  * U ) .* U;
  T2 = (U' * U ) .* U;
  T3 = (U  * U') .* U;
  C = T1 + T2 + T3;
  W = C + C';
end

function W = M6(A)
  [B, U, ~] = DirectionalBreakup(A);
  C1 = (U * B) .* U;
  C1 = C1 + C1';
  C2 = (U' * U) .* B;
  W = C1 + C2;
end

function W = M7(A)
  [B, U, ~] = DirectionalBreakup(A);
  C1 = (U' * B) .* U';
  C1 = C1 + C1';
  C2 = (U * U') .* B;
  W = C1 + C2;
end

function W = M8(A)
  [~, U, G] = DirectionalBreakup(A);
  W = zeros(size(G));
  N = size(G, 1);
  for i = 1:N
    J = find(U(i, :));
    for j1 = 1:length(J)
        for j2 = (j1+1):length(J)
            k1 = J(j1);
            k2 = J(j2);
            if A(k1, k2) == 0 && A(k2, k1) == 0
               W(i, k1)  = W(i, k1) + 1; 
               W(i, k2)  = W(i, k2) + 1;
               W(k1, k2) = W(k1, k2) + 1;
            end
        end
    end
  end
  W = sparse(W + W');
end

function W = M9(A)
  [~, U, G] = DirectionalBreakup(A);
  W = zeros(size(G));
  N = size(G, 1);
  for i = 1:N
    J1 = find(U(i, :));
    J2 = find(U(:, i));
    for j1 = 1:length(J1)
      for j2 = 1:length(J2)
        k1 = J1(j1);
        k2 = J2(j2);
        if A(k1, k2) == 0 && A(k2, k1) == 0
          W(i, k1)  = W(i, k1) + 1; 
          W(i, k2)  = W(i, k2) + 1;
          W(k1, k2) = W(k1, k2) + 1;
        end
      end
    end
  end
  W = sparse(W + W');
end

function W = M10(A)
  W = M8(A');
end

function W = M11(A)
  [B, U, G] = DirectionalBreakup(A);
  W = zeros(size(G));
  N = size(G, 1);
  for i = 1:N
    J1 = find(B(i, :));
    J2 = find(U(i, :));
    for j1 = 1:length(J1)
      for j2 = 1:length(J2)
        k1 = J1(j1);
        k2 = J2(j2);
        if A(k1, k2) == 0 && A(k2, k1) == 0
          W(i, k1)  = W(i, k1) + 1; 
          W(i, k2)  = W(i, k2) + 1;
          W(k1, k2) = W(k1, k2) + 1;
        end
      end
    end
  end
  W = sparse(W + W');
end

function W = M12(A)
  W = M11(A');
end  

function W = M13(A)
  [B, ~, G] = DirectionalBreakup(A);
  W = zeros(size(G));
  N = size(G, 1);
  for i = 1:N
    J = find(B(i, :));
    for j1 = 1:length(J)
        for j2 = (j1+1):length(J)
            k1 = J(j1);
            k2 = J(j2);
            if A(k1, k2) == 0 && A(k2, k1) == 0
               W(i, k1)  = W(i, k1) + 1; 
               W(i, k2)  = W(i, k2) + 1;
               W(k1, k2) = W(k1, k2) + 1;
            end
        end
    end
  end
  W = sparse(W + W');
end

function W = Bifan(A)
  [~, U, G] = DirectionalBreakup(A);
  NA = ~A & ~A';
  W = zeros(size(G));
  [ai, aj] = find(triu(NA, 1));
  for ind = 1:length(ai)
    x = ai(ind);
    y = aj(ind);
    xout = find(U(x, :));
    yout = find(U(y, :));
    common = intersect(xout, yout);
    nc = length(common);
    for i = 1:nc
      for j = (i+1):nc
        w = common(i);
        v = common(j);
        if NA(w, v) == 1
          W(x, y) = W(x, y) + 1;
          W(x, w) = W(x, w) + 1;
          W(x, v) = W(x, v) + 1;
          W(y, w) = W(y, w) + 1;
          W(y, v) = W(y, v) + 1;
          W(w, v) = W(w, v) + 1;
        end
      end
    end
  end
  W = sparse(W + W');
end

