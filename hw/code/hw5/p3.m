z = randn(10000,1);
X = [z.^0,z.^1,z.^2,z.^3,z.^4,z.^5,z.^6,z.^7,z.^8];
y = cos(z);
global a1 
a1 = lscov(X,y);
efx = gausshq(@f);

function fx = f(x)
global a1
fx = ([x.^0,x,x.^2,x.^3,x.^4,x.^5,x.^6,x.^7,x.^8]*a1+cos(x)).^2;
end