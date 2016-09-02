function D = spdiag(v)
  D = spdiags(v, 0, length(v), length(v));
