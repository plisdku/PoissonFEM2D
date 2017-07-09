function q = diffPoly(p)
% diffPoly  Return the derivative of a polynomial
%
% Example: the derivative of x^2 is 2x
%   diffPoly([1 0 0])     returns [2 0]

q = fliplr(1:numel(p)-1) .* p(1:end-1);
