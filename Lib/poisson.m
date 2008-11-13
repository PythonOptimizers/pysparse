function poisson(n)
%fprintf('poisson2d      ');
%tic;
%  poisson2d(n);
%toc;
%fprintf('poisson2d_blk  ');
%tic;  
%  poisson2d_blk(n);
%toc;
fprintf('poisson2d_kron ');
tic;
  poisson2d_kron(n);
toc;
return;


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
return;

function L = poisson2d_blk(n)
e = ones(n,1);
P = spdiags([-e 4*e -e], [-1 0 1], n, n);
I = -speye(n);
L = sparse(n*n);
for i = 1:n:n*n
  L(i:i+n-1,i:i+n-1) = P;
  if i > 1, L(i:i+n-1,i-n:i-1) = I; end
  if i < n*n - n, L(i:i+n-1,i+n:i+2*n-1) = I; end
end
return;

function L = poisson2d_kron(n)
e = ones(n,1);
P = spdiags([-e 2*e -e], [-1 0 1], n, n);
L = kron(P, speye(n)) + kron(speye(n), P);
return;


