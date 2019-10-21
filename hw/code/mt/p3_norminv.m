% Compute the norm of the tridiagonal inverse
%
function [normTinv] = p3_norminv(alpha, beta)

n = length(alpha);
eye_1 = [1 zeros(1,n-1)]';
eye_n = [zeros(1,n-1) 1]';
firstcol = zeros(n,1);
lastcol = zeros(n,1);
diaginv = zeros(n,1);

%% Compute first and last columns of inverse (Tridiagonal matrix algo)
% Copy original alpha
alpha1 = alpha;
for i = 1:n-1
    w = beta(i)/alpha1(i);
    alpha1(i+1) = alpha1(i+1)-w*beta(i);
    eye_1(i+1) = eye_1(i+1)-w*eye_1(i);
    eye_n(i+1) = eye_n(i+1)-w*eye_n(i);
end

firstcol(n) = eye_1(n)/alpha1(n);
lastcol(n) = eye_n(n)/alpha1(n);

for i = n-1:-1:1
    firstcol(i) = (eye_1(i)-beta(i)*firstcol(i+1))/alpha1(i);
    lastcol(i) = (eye_n(i)-beta(i)*lastcol(i+1))/alpha1(i);
end

%% Compute diagonal of inverse
diaginv(1) = firstcol(1);
diaginv(n) = lastcol(n);

for i = 2:n-1
    leftel = lastcol(i-1)*firstcol(i)/firstcol(n);
    rightel = firstcol(i+1)*lastcol(i)/lastcol(1);
    diaginv(i) = (1-beta(i-1)*leftel-beta(i)*rightel)/alpha(i);
end

%% Compute column sums of abs(S)
firstcol = abs(firstcol);
diaginv = abs(diaginv);
lastcol = abs(lastcol);

colsum=zeros(n,1);
supcolsum = colsum;
% column sum of subdiagonal
colsum(n-1)=firstcol(n);
for i = n-1:-1:2 
    colsum(i-1)=colsum(i)+firstcol(i);
end
colsum=colsum.*lastcol/lastcol(1);
% do likewise for superdiag
supcolsum(2)=lastcol(1);
for i = 2:n-1
    supcolsum(i+1)=supcolsum(i)+lastcol(i);
end
supcolsum=supcolsum.*firstcol/firstcol(n);

%% Sum entire column and find max
normTinv = max(colsum+diaginv+supcolsum);