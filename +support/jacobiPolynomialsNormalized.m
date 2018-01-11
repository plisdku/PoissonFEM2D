function polyCoeffs = jacobiPolynomialsNormalized(N, alpha, beta)
% jacobiPolynomials  Return coefficients of first N Jacobi polynomials
%
% jacobiPolynomials(1, 0, 0) returns one polynomial, the first (0,0) Jacobi
% polynomial, equal to L0, the first Legendre polynomial.
%
% jacobiPolynomials(N, 0, 0) returns the first N Legendre polynomials
%
% jacobiPolynomials(N, a, b) returns the first N Jacobi polynomials with
% weight parameters (a, b).
%

a_prefactor = @(n) 2./(2*n + alpha + beta);
a_numer = @(n) sqrt(n .* (n+alpha+beta).*(n+alpha).*(n+beta));
a_denom = @(n) sqrt((2*n + alpha + beta - 1).*(2*n + alpha + beta + 1));
a = @(n) a_prefactor(n) .* a_numer(n) ./ a_denom(n);

% coefficient "b"

b_numer = @(n) -(alpha^2 - beta^2);
b_denom = @(n) (2*n + alpha + beta).*(2*n + alpha + beta + 2);
b = @(n) b_numer(n)./b_denom(n);

% P0 and P1, the first two polynomials

P0 = sqrt( 2^(-alpha-beta-1) * gamma(alpha+beta+2)/gamma(alpha+1)/gamma(beta+1));

P1 = 0.5*P0 * sqrt((alpha+beta+3)/(alpha+1)/(beta+1)) * ...
    [alpha + beta + 2, alpha - beta];

% Now we recurse:

aa = a(1:N);
bb = b(1:N);

polyCoeffs = zeros(N);

polyCoeffs(1,end) = P0;
if N > 1
    polyCoeffs(2, end-1:end) = P1;
end
%%

for nn = 3:N
    polyCoeffs(nn,:) = (1/aa(nn-1)) * ...
        ( [polyCoeffs(nn-1,2:end), 0] - bb(nn-2)*polyCoeffs(nn-1,:) ...
        - aa(nn-2)*polyCoeffs(nn-2,:));
end
