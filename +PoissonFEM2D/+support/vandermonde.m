function V = vandermonde(numNodes, xx)
% vandermonde   Create vandermonde matrix for evaluating normalized Legendre
% polynomials
%
% vandermonde(xx) uses length(xx) polynomials
% vandermonde(xx, numNodes) uses numNodes polynomials

import PoissonFEM2D.support.*

if nargin == 1
    %numNodes = numel(xx);
    xx = nodes1d(numNodes);
end


l = legendrePolynomialsNormalized(numNodes);

V = zeros(numel(xx), numNodes);

for cc = 1:numNodes
    V(:,cc) = polyval(l(cc,:), xx);
end

