% Solve the shifted Hessenberg system (H-s*I)*x = b in O(n^2)
%
function [x] = hw7hess(H, sigma, b)

% Code adapted from 10/28 notes

n = length(H);
H = H-sigma*eye(n);
Q = eye(n);
V = zeros(2,n-1);
x = zeros(n,1);

% Compute the QR factorization with Householder
for j = 1:n-1
    
    % -- Find W_j = I-2vv' to put zero into H(j+1,j)
    u      = H(j:j+1,j);
    u(1)   = u(1) + sign(u(1))*norm(u);
    v      = u/norm(u);
    V(:,j) = v;
    
    % -- H := W_j H
    H(j:j+1,:) = H(j:j+1,:)-2*v*(v'*H(j:j+1,:));
    
end

% Compute Q
for j = 1:n-1
    
    % -- H := WHW', Q := QW
    v = V(:,j);
%     H(:,j:j+1) = H(:,j:j+1)-(H(:,j:j+1)*(2*v))*v';
    Q(:,j:j+1) = Q(:,j:j+1)-(Q(:,j:j+1)*(2*v))*v';

    
end

% Backwards solve
y = Q'*b;
x(n) = y(n)/H(n,n);
for j = n-1:-1:1
    x(j) = (y(j)-H(j,j+1:end)*x(j+1:end))/H(j,j);
end
