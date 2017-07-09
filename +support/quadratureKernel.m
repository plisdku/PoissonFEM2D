function Q = quadratureKernel(N, rr, varargin)
% Q = quadratureKernel(N, rr, ss)
% Q = quadratureKernel(N, rr, ss, jacobian)

if nargin <= 3
    jacobian = 1;
else
    jacobian = varargin{1};
end

detJ = abs(jacobian);

V = support.vandermonde(rr, N);
V_inv = inv(V);

Q = V_inv' * V_inv * detJ;

