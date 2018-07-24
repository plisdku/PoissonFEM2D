function l = legendrePolynomials(N)
% legendrePolynomials    Return coefficients of first N Legendre
% polynomials
%
% legendrePolynomials(1) returns one polynomial, [1]
%
% legendrePolynomials(2) returns the polynomials {1, x} in a matrix:
%
%   [ 0 1 
%     1 0 ]
%

l = zeros(N,N);

l(1,end) = 1; % First poly is 1

if N == 1
    return
end

l(2,[end-1, end]) = [1 0];

for ll = 3:N
    
    n = ll - 1;
    
    l(ll,:) = (2*n-1)/n * conv(l(ll-1,2:end), [1 0]) - ...
        (n-1)/n * l(ll-2,:);
end



