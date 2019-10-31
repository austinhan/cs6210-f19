function [sigma, gs] = hw7newton(A, b, c, sigma)

  n = length(b);
  I = eye(n);
  en = zeros(n+1,1);
  en(end) = 1;

  gs = zeros(4,1);
  for k = 1:4

    M = [A-sigma*I, b; c', 0];
    [L, U] = lu(M);

    % TODO:
    %   Compute g(sigma) and save to gs(k), then 
    %   compute a Newton step to update the
    %   eigenvalue estimate sigma.
    x = M\en;
    gs(k)=x(end);
%     ddet = trace(adjoint(A-sigma*I));
%     dg = ddet/det(-b*c');
%     sigma=sigma-gs(k)/dg;
    M = [A-(sigma+1e-6)*I, b; c', 0];
    x2 = M\en;
    M = [A-(sigma-1e-6)*I, b; c', 0];
    x1 = M\en;
    dgs = (x2(end)-x1(end))/2e-6;
    sigma = sigma - gs(k)/dgs;


  end
