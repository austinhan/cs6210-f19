% p4tester

m = 3;
n = 5;

A = triu(randn(n));
A = A+A'-diag(diag(A));
M = randn(n);
M = M'*M;
% M = eye(n);
rankc = 0;
while rankc ~= m
    C = randn(m,n);
    rankc = rank(C);
end
b = randn(m,1);
xbar = p4minx(M,C,b);
V = p4null(M,C);
[x,feasible]=p4optimize(A,M,C,b);
success = abs(x'*M*x-1)<1e-10;