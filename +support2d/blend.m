function y = blend(b1, b2, b3)
% y = blend(b1, b2, b3) in barycentric coords returns the value of a func
% which is 1 when b1 == 0 and 0 when b2 == 0 or b3 == 0.
% 
% Like H&W I put in the extra factors of 2 that seem sort of extraneous.
%

y = 4.*b2.*b3 ./ ( (2*b2 + b1).*(2*b3 + b1) );
y(isinf(y) | isnan(y)) = 0;