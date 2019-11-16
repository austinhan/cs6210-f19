% Differentiate the solution to (I+ZZ')x = b
%
function [dx] = p2_deriv_solve(Q, R, x, b, dZ, db)
  % Complete this code
  k = size(R,1);
  db1 = db-dZ*(R'*(Q'*x))-Q*R*(dZ'*x);
  dx =db1-Q*R*((eye(k)+R'*R)\(R'*(Q'*db1)));