function [xx, ww] = gaussLobatto(nPoints)
% gaussLobatto   Generate abscissae and weights for Gauss-Lobatto
% quadrature
%
% [xx ww] = gaussLobatto(nPoints)

%nPoints = 3;

import support.*

l = legendrePolynomials(nPoints);

%% The abscissae are the zeros of P'_(n-1)

abscissae = roots(diffPoly(l(end,:)));

rowVec = @(A) reshape(A, 1, []);

xx = [-1, rowVec(sort(abscissae)), 1];

%% The weights are some funky formula
% Thanks, Wikipedia.  I should study this myself, I realize.

ww = 0*xx;

for nn = 2:nPoints-1
    ww(nn) = 2./(nPoints*(nPoints-1)*paulyval(l(end,:), xx(nn))^2);
end

ww([1 end]) = 1 - 0.5*sum(ww);
