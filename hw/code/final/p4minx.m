% Minimize x'*M*x s.t. Cx = b
%
function [xbar] = p4minx(M, C, b)
m = size(C,1);
n = size(C,2);

A = [2*M -C';
     C   zeros(m,m)];
b = [zeros(n,1);b];

sol = A\b;
xbar = sol(1:n,1);