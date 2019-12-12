% Return an M-orthonoormal basis V for the null space of C, i.e.
%   V'*M*V = I
%   C*V = 0
%
function [V] = p4null(M, C)

m = size(C,1);
n = size(M,1);

V = null(C); % or SVD
V(:,1) = V(:,1)/sqrt(V(:,1)'*M*V(:,1));
for i = 2:n-m
    v = V(:,i);
    v = v - V(:,1:i-1)*V(:,1:i-1)'*M*v;
    V(:,i) = v/sqrt(v'*M*v);
end