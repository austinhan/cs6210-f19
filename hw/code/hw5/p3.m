z = randn(5e6,1);
X = [z.^0,z.^1,z.^2,z.^3,z.^4,z.^5,z.^6,z.^7,z.^8];
y = cos(z);
global a1 
% a1 = lscov(X,y);
R =chol(X'*X, 'upper');
a1 = R\(R'\(X'*y));
% a1 = a1+1e-4*randn(9,1);
% a1 = [1 0 -1/factorial(2) 0 1/factorial(4) 0 -1/factorial(6) 0 1/factorial(8)]';% Taylor series
efx = gausshq(@f);

function fx = f(x)
global a1
fx = ([x.^0,x,x.^2,x.^3,x.^4,x.^5,x.^6,x.^7,x.^8]*a1-cos(x)).^2;
end