% Compute x = A\b quickly
%
function [x] = p1_fast_solve(Q, R, b)

  n = length(b);
  k = size(R,1);
  Z = Q*R;
  A = eye(n) + Z*Z';
%   x = A\b;
  x = b-Z*((eye(k)+R'*R)\(Z'*b));
