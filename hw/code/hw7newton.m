function [sigma, gs] = hw7newton(A, b, c, sigma)

  n = length(b);
  I = eye(n);
  en = zeros(n+1,1);
  en(end) = 1;

  gs = zeros(4,1);
  dM = [-I, 0*b;0*c',0];
  for k = 1:4

    M = [A-sigma*I, b; c', 0];
%     [L, U] = lu(M);
    luM = decomposition(M,'lu');

    % TODO:
    %   Compute g(sigma) and save to gs(k), then 
    %   compute a Newton step to update the
    %   eigenvalue estimate sigma.
%     y = L\en;
%     x = U\y;
    x = luM\en;
    gs(k)=x(end);

    % Solve for g'(sigma)
%     y = L\(dM*x);
%     dx = -U\y;
    dx = -luM\(dM*x);
    dgs = dx(end);
    % Update
    sigma = sigma - gs(k)/dgs;


  end
