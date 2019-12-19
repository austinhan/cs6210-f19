% Minimize x'*A*X s.t. x'*M*x = 1 and C*x = b.
% Return feasible=0 if no points satisfy the constraints.
% (if feasible=0, return x=xbar)
function [x, feasible] = p4optimize(A, M, C, b)

xbar = p4minx(M,C,b);

if xbar'*M*xbar > 1
    feasible = 0;
    x = xbar;
    return
else
    feasible = 1;
    V = p4null(M,C);
    k = size(V,2); % dim of null space
    a = V'*A*V;
    b = V'*A*xbar;
    d = 1/(1-xbar'*M*xbar);
    mu = polyeig(a^2-d*b*b',-2*a,eye(k));
    mu = sort(real(mu(imag(mu)==0)));
    lambda = mu(1);
    y = -(V'*A*V-lambda*eye(k))\(V'*A*xbar);
    x = xbar + V*y;
end

