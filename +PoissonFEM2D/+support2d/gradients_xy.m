function [Dx, Dy] = gradients_xy(N, rr, ss, varargin)
% [Dx, Dy] = partialDerivativeOperators_xy(N, rr, ss)
% [Dx, Dy] = partialDerivativeOperators_xy(N, rr, ss, invJac)
import PoissonFEM2D.*
if nargin < 4
    invJac = eye(2);
else
    invJac = varargin{1};
end

[Dr, Ds] = support2d.gradients(N, rr, ss);

Dx = Dr*invJac(1,1) + Ds*invJac(2,1);
Dy = Dr*invJac(1,2) + Ds*invJac(2,2);
