import PoissonFEM2D.support.*

isClose_tol = @(a,b,tol) abs(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Test P0, the first Legendre polynomial

l = legendrePolynomials(1);

assert(isequal(l, 1));

%% Test P0 and P1

l = legendrePolynomials(2);
assert(isequal(l, [0 1; 1 0]));

%% Test P3 and P4 too

l = legendrePolynomials(4);

expectation = [0  0  0  1
               0  0  1  0
               0  3/2  0  -1/2
               5/2  0  -3/2  0];

assert(isequal(l, expectation));

%% Test normalized polynomials
% Make sure they integrate to 1.  Use Matlab's built-in adaptive quadrature
% integration.  It should be exact on polynomials, right?

l = legendrePolynomialsNormalized(10);
legendreSquaredFunc = @(p) @(x) polyval(l(p,:), x).^2;

assert(isClose(quad(legendreSquaredFunc(1), -1, 1), 1));
assert(isClose(quad(legendreSquaredFunc(2), -1, 1), 1));
assert(isClose(quad(legendreSquaredFunc(3), -1, 1), 1));

%% Done!

disp('Legendre polynomial test PASSED');
