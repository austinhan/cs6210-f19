% Differentiate the Cholesky factor of M
%
function [dR] = p2_deriv_chol(M, R, dM)
  % Complete this code
    A = triu(R'\(dM/R));
    A = A - 1/2*diag(diag(A));
    dR=A*R;