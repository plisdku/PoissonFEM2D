function l = legendrePolynomialsNormalized(N)
% legendrePolynomialsNormalized(numPolynomials)
% Normalize polynomials so their squared integral is 1
%%
import support.*

l = legendrePolynomials(N);

assert(~any(imag(l(:))));

for nn = 1:N
    % It turns out that using conv on the polynomial with itself results in
    % inaccurate results verrrry quickly (for around N = 20).
    %poly2 = conv(l(nn,:), l(nn,:));
    %assert(~any(imag(poly2(:))));
    %i1(nn) = intPoly(poly2, -1, 1);
    %i2(nn) = quad(@(x) polyval(l(nn,:),x).^2, -1, 1);
    
    i2 = quad(@(x) polyval(l(nn,:),x).^2, -1, 1);
    l(nn,:) = l(nn,:) / sqrt(i2);
end

if any(imag(l(:)))
    error('Failure to normalize polynomials.');
end

