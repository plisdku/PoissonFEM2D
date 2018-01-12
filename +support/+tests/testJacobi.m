import support.*

isClose_tol = @(a,b,tol) norm(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Verify that alpha = beta = 0 produces Legendre polynomials
% These coefficients actually diverge more and more for higher-order
% polynomials, but the divergence is slow.  Probably not important for my
% purposes.  (Someday project: learn how to calculate the coefficients more
% accurately.)

ll = legendrePolynomialsNormalized(5);
jj = jacobiPolynomialsNormalized(5, 0, 0);

assert(isClose_tol(ll,jj,1e-5));

%% Test orthonormality

alpha = 0.2;
beta = -0.3;

%alpha = 0.0;
%beta = 1.0;

jj = jacobiPolynomialsNormalized(5, alpha, beta);

innerProds = zeros(5);
for nn = 1:5
    for mm = 1:5
        weightFunc = @(x) (1-x).^alpha .* (1+x).^beta;
        integrand = @(x) weightFunc(x) .* polyval(jj(nn,:), x) .* ...
            polyval(jj(mm,:), x);
        
        innerProds(nn,mm) = quadgk(integrand, -1, 1);
    end
end

disp(innerProds)

%% Done!

disp('Jacobi polynomial test PASSED');
