% Return an approximate eigenvector v s.t. r = A*v-mu*v is minimal (in norm)
% where v is from the space spanned by V (V'*V = I).  Also return the
% norm of the min backward error E s.t. (A+E)*v = mu*v.
%
function [v, normE] = p3eigen(A, mu, V)
n=size(A,1);

Abar=(A-mu*eye(n))*V;
[~,~,V1]=svd(Abar,'econ');
y=abs(V1(:,end));

% Equivalent:
% [y,~]=eigs(Abar'*Abar,1,'smallestabs');
v=V*y;

% normE <= largest singular value of (mu*I-A)
[~,S,~]=svd(A-mu*eye(n));
normE=S(1,1);