% p4tester

m = 3;
n = 5;

A = rand(n);
M = A'*A;
rankc = 0;
while rankc ~= m
    C = rand(m,n);
    rankc = rank(C);
end
b = rand(m,1);
xbar = p4minx(M,C,b);
V = p4null(M,C);