function Q = quadratureKernel(N, rr, varargin)
% Q = quadratureKernel(N, rr, ss)
% Q = quadratureKernel(N, rr, ss, jacobian)
import PoissonFEM2D.*

if nargin <= 3
    jacobian = 1;
else
    jacobian = varargin{1};
end

detJ = abs(jacobian);

V = support.vandermonde(N, rr);
V_inv = inv(V);

Q = V_inv' * V_inv * detJ;

