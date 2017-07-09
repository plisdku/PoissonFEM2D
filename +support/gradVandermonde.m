function DV = gradVandermonde(xx, numNodes)
% gradVandermonde   Create gradient of Vandermonde matrix
%
% gradVandermonde(xx)
% gradVandermonde(xx, numNodes)

if nargin == 1
    numNodes = numel(xx);
end

import support.*

l = legendrePolynomialsNormalized(numNodes);

DV = zeros(numel(xx), numNodes);

for cc = 1:numNodes
    Dl = diffPoly(l(cc,:));
    DV(:,cc) = polyval(Dl, xx);
end
