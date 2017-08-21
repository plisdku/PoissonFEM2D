function V = vandermonde(numNodes, xx)
% vandermonde   Create vandermonde matrix for evaluating normalized Legendre
% polynomials
%
% vandermonde(xx) uses length(xx) polynomials
% vandermonde(xx, numNodes) uses numNodes polynomials

if nargin == 1
    %numNodes = numel(xx);
    xx = support.nodes1d(numNodes);
end

import support.*

l = legendrePolynomialsNormalized(numNodes);

V = zeros(numel(xx), numNodes);

for cc = 1:numNodes
    V(:,cc) = polyval(l(cc,:), xx);
end

