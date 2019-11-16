% Sanity checks -- I recommend checking that
%  a.  The multipliers correspond to stationary points
%  b.  The solution returned by qcqp_solve is a constrained minimizer
%  c.  For some reference problem, qcqp_solve returns the global min

n = 10;

% A = eye(n);
% b = zeros(n,1);
% [x,phi] = hw8solve(A,b);

A = 2*rand(n) - 1;
A = triu(A);
A = A + A' - diag(diag(A));
b = 2*rand(n,1)-1;
[x,mu,phi] = hw8solve(A,b);

% a.
x_test = (A-mu*eye(n))\b;

% b.
eigvals = eig(A);
fprintf('mu < lambda_2(A): %i\n',mu<eigvals(2));
% c.
% phi = x^2+1/3*x*y+10*y^2-2*x
% minimum from Wolfram Alpha
Ac = 2*[1 1/6;1/6 10];
bc = [2;0];
[xc,muc,phic] = hw8solve(Ac,bc);
fprintf('reference min: %d\n',-1.00278);
fprintf('calculated min: %d\n',phic);

