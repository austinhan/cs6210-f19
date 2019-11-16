% Differentiate the squared M-norm of x
%
function [dnormx] = p2_deriv_mnorm(x, M, dx, dM)
  % Complete this code
dnormx = dx'*M*x+x'*dM*x+x'*M*dx;