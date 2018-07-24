function Q = quadratureKernel(N, rr, ss, varargin)
% Q = quadratureKernel(N, rr, ss)
% Q = quadratureKernel(N, rr, ss, jacobian)
import PoissonFEM2D.*
if nargin <= 3
    jacobian = eye(2); 
else
    jacobian = varargin{1};
end

detJ = det(jacobian);

V = support2d.vandermonde(N, rr, ss);
V_inv = inv(V);

Q = V_inv' * V_inv * detJ;

