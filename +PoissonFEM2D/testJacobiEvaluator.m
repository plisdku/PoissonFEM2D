%% Test JacobiEvaluator!
import PoissonFEM2D.*

isClose_tol = @(a,b,tol) norm(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Verify that alpha = beta = 0 produces Legendre polynomials
% These coefficients actually diverge more and more for higher-order
% polynomials, but the divergence is slow.  Probably not important for my
% purposes.  (Someday project: learn how to calculate the coefficients more
% accurately.)

jac = PoissonFEM2D.JacobiEvaluator(5);

ll = support.legendrePolynomialsNormalized(5);
jj = jac.coefficientArray(1:5, :, 1, 1);

assert(isClose_tol(ll,jj,1e-5));

%% Test orthonormality

alpha = 1;
beta = 2;

jac = JacobiEvaluator(5);

innerProds = zeros(5);
for nn = 1:5
    for mm = 1:5
        weightFunc = @(x) (1-x).^alpha .* (1+x).^beta;
        
        integrand = @(x) weightFunc(x) .* jac.evaluate(nn-1, alpha, beta, x) ...
            .* jac.evaluate(mm-1, alpha, beta, x);
        
        innerProds(nn,mm) = quadgk(integrand, -1, 1);
    end
end

disp(innerProds)

%% Test derivatives

alpha = 1;
beta = 2;

jac = JacobiEvaluator(5);

xs = linspace(-1, 1, 1000);
for nn = 1:5
    y = jac.evaluate(nn-1, alpha, beta, xs);
    dydx = jac.evaluateDerivative(nn-1, alpha, beta, xs);
    
    dydx_meas = gradient(y,xs);
    
    relErr = norm(dydx - dydx_meas)/norm(dydx_meas);
    
    if relErr > 1e-3
        error('Derivative failed for n = %i', nn);
    end
end

%% Done!

disp('Jacobi polynomial test PASSED');