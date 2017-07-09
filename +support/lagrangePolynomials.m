function coeffs = lagrangePolynomials(nodes)
% lagrangePolynomials(nodes)  Return all the coefficients of the Lagrange
% polynomials for the given nodes.
%
% lp = lagrangePolynomials([0 1])
%
% polyval(lp(1,:), 0)  % returns 1
% polyval(lp(1,:), 1)  % returns 0
% polyval(lp(2,:), 0)  % returns 0
% polyval(lp(2,:), 1)  % returns 1

coeffs = zeros(numel(nodes));

for n = 1:numel(nodes)
    zeds = nodes([1:n-1, n+1:end]);
    coeffs(n, :) = poly(zeds)/polyval(poly(zeds), nodes(n));
end

