function L = poisson2d(n)
  L = sparse(n*n);
  for i = 1:n
    for j = 1:n
      k = i + n*(j-1);
      L(k,k) = 4;
      if i > 1, L(k,k-1) = -1; end
      if i < n, L(k,k+1) = -1; end
      if j > 1, L(k,k-n) = -1; end
      if j < n, L(k,k+n) = -1; end
    end
  end
