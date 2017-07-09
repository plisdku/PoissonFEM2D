function [z, dz_dr, dz_ds] = jacobiNormalized2D(N, ii,jj, xs, ys, varargin)
% Return values of a Jacobi polynomial of order N at the given points.
% x and y should be in the big simplex in [-1,1]x[-1,1].

import support.*

jac_i = jacobiPolynomialsNormalized(N,0,0);
jac_j = jacobiPolynomialsNormalized(N, 2*ii+1, 0);

% Evaluate two Legendre polynomials (one is fancier than the other)
P_i = @(x) polyval(jac_i(ii+1,:), x);
P_j = @(x) polyval(jac_j(jj+1,:), x);

% Evaluate the derivatives of the Legendre polynomials
DP_i = @(x) polyval(diffPoly(jac_i(ii+1,:)), x);
DP_j = @(x) polyval(diffPoly(jac_j(jj+1,:)), x);

% this function does have a problem when s = 1... hmm!!!
a = @(r,s) 2*(1+r)./(1-s) - 1;
b = @(r,s) s;

Da_dr = @(r,s) 2./(1-s);
Da_ds = @(r,s) 2*(1+r)./(1-s).^2;
Db_dr = @(r,s) 0*r;
Db_ds = @(r,s) ones(size(s));

%da_dr = @(r,s) 

basisFunc = @(r,s) sqrt(2.0) .* P_i(a(r,s)) .* P_j(b(r,s)) .* (1-b(r,s)).^ii;

% its partial derivatives:
dBF_dr = @(r,s) ...
    sqrt(2.0) * DP_i(a(r,s)).*Da_dr(r,s).*P_j(b(r,s)).*(1-b(r,s)).^ii; % some terms are zero
dBF_ds = @(r,s) ...
    sqrt(2.0) * DP_i(a(r,s)).*Da_ds(r,s).*P_j(b(r,s)).*(1-b(r,s)).^ii ...
    + sqrt(2.0) * P_i(a(r,s)).*DP_j(b(r,s)).*Db_ds(r,s).*(1-b(r,s)).^ii ...
    + sqrt(2.0) * P_i(a(r,s)) .* P_j(b(r,s)) .* (-ii * (1-s).^(ii-1));


% Here's how we avoid the stupid s == 1 problem!  Well done Paul!
ys(ys == 1) = 1 - 1e-9;

z = basisFunc(xs, ys);
dz_dr = dBF_dr(xs,ys);
dz_ds = dBF_ds(xs,ys);
