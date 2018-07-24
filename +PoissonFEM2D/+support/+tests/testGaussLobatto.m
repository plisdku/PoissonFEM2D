import PoissonFEM2D.support.*

isClose_tol = @(a,b,tol) norm(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Test interval

[xx ww] = gaussLobatto(2);
assert(isequal(xx, [-1 1]));
assert(isequal(ww, [1 1]));

%% Test three nodes vs. tabulated answer

[xx ww] = gaussLobatto(3);
assert(isequal(xx, [-1 0 1]));
assert(isClose(ww, [1/3, 4/3, 1/3]));

%% Test four nodes vs. tabulated answer

[xx ww] = gaussLobatto(4);
assert(isClose(xx, [-1, -sqrt(1/5), sqrt(1/5), 1]));
assert(isClose(ww, [1/6 5/6 5/6 1/6]));

%% Test integration rules
% Gauss-Lobatto quadrature is accurate for polynomials up to degree 2n-3

poly3 = [1 -1 2 4];
poly7 = 1:8;

%% The 3-node rule is exact up to degree 3

[xx ww] = gaussLobatto(3);
glIntegral = sum(ww .* polyval(poly3, xx));
assert(isClose(glIntegral, intPoly(poly3, -1, 1)));

%% The 5-node rule is exact up to degree 7
% Test as well on poly3

[xx ww] = gaussLobatto(5);

glIntegral = sum(ww .* polyval(poly3, xx));
assert(isClose(glIntegral, intPoly(poly3, -1, 1)));

glIntegral = sum(ww .* polyval(poly7, xx));
assert(isClose(glIntegral, intPoly(poly7, -1, 1)));

%% Done!

disp('Gauss-Lobatto quadrature rules PASSED');

