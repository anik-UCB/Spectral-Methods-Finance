function r = lagroots(N);

%  The function r = lagroots(N) computes the roots of the 
%  Laguerre polynomial of degree N.


J = diag([1:2:2*N-1])-diag([1:N-1],1)-diag([1:N-1],-1);  % Jacobi matrix
r = sort(eig(sparse(J)));                                % Compute eigenvalues