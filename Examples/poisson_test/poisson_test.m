n = 500;
tic
L = poisson2d_kron(n);
toc
tic
[x,flag,relres,iter] = pcg(L, ones(n*n,1), 1e-12, 2000, [], [], ...
			   zeros(n*n,1));
toc
