function [xp,fxp] = goldensection(f,a1,b1,delta)
rho = (3-sqrt(5))/2;
L1 = b1-a1;
N = ceil(log(delta/L1)/log(0.61803))

a = NaN*zeros(N+1,1);
b = NaN*zeros(N+1,1);
alpha = NaN*zeros(N+1,1);
beta = NaN*zeros(N+1,1);
falpha = NaN*zeros(N+1,1);
fbeta = NaN*zeros(N+1,1);
L = NaN*zeros(N+1,1);

a(1) = a1;
b(1) = b1;
L(1) = L1;

for k=1:N
    alpha(k) = a(k)+rho*(b(k)-a(k));
    beta(k) = a(k)+(1-rho)*(b(k)-a(k));
    falpha(k) = f(alpha(k));
    fbeta(k) = f(beta(k));
    if falpha(k)<fbeta(k)
        a(k+1) = a(k);
        b(k+1) = beta(k);
    else
    a(k+1) = alpha(k);
    b(k+1) = b(k);
    end
    L(k+1) = b(k+1)-a(k+1);
end
xp = a(N+1)+(b(N+1)-a(N+1))/2;
fxp = f(xp);
[a,b,alpha,beta,falpha,fbeta,L]
end