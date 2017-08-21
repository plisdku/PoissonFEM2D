% Test vandermonde matrix

import support.*

isClose_tol = @(a,b,tol) norm(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Test small one value-by-value

xx = [-1 0 0.5]; % Node positions, weird on purpose

v = vandermonde(length(xx), xx);

l = legendrePolynomialsNormalized(3);

for pp = 1:3
    
    polynomial = l(pp,:);
    
    p_of_x = polyval(polynomial, xx);
    
    assert(isClose(v(:,pp), p_of_x'));
end

%% Test use of smaller Vandermonde matrix

xx = [-3 -1 0 2];

v = vandermonde(length(xx), xx);
v1 = vandermonde(4, xx(1));
vEnd = vandermonde(4, xx(end));

assert(isClose(v(1,:), v1));
assert(isClose(v(end,:), vEnd));

%% Test gradient of Vandermonde

xx = [-3 -2 0 1];

dx = 1e-6;
xx2 = xx + dx;

V1 = vandermonde(length(xx), xx);
V2 = vandermonde(length(xx2), xx2);
DV_diff = (V2-V1)/dx;
DV = gradVandermonde(length(xx), xx);

assert(isClose_tol( DV, DV_diff, 1e-4 ) );

%% Done

disp('Vandermonde matrix tests PASSED')