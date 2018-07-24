function ii = intPoly(p, x0, x1)
% intPoly   Integrate a polynomial
%
% intPoly(p) returns the coefficients of an antiderivative of p
%
% intPoly(p, x0, x1) integrates the polynomial p from x0 to x1

if nargin == 3
    ii = sum(p .* (x1.^fliplr(1:numel(p))) ./ fliplr(1:numel(p))) - ...
        sum(p .* (x0.^fliplr(1:numel(p))) ./ fliplr(1:numel(p)));
else
    ii = [p ./ (numel(p):-1:1), 0];
end