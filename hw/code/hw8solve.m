% [x] = hw8solve(A, b)
%
% Return a global minimizer for the constrained optimization problem
%   minimize x'Ax/2 - x'b s.t. x'x = 1
% Your code should run in an overall time of O(n^3).
%
function [x,mu_x,phi] = hw8solve(A, b)

n = length(b);
% Find all possible multipliers (O(n^3) time)
mu = polyeig(A^2-b*b',-2*A,eye(n));
mu = sort(real(mu(imag(mu)==0)));

x = (A-mu(1)*eye(n))\b;
phi = 1/2*x'*A*x-x'*b;
mu_x = mu(1);

for i = 2:length(mu)
    x_test = (A-mu(i)*eye(n))\b;
    phi_test = 1/2*x_test'*A*x_test-x_test'*b;
    if phi_test < phi
        phi = phi_test;
        x = x_test;
        mu_x = mu(i);
    end
end

% for i = 1:length(mu)
%     x(:,i) = (A-mu(i)*eye(n))\b;
%     phi(i) = 1/2*x(:,i)'*A*x(:,i)-x(:,i)'*b;
% end
% i=1;