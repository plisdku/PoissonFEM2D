function [Dx, Dy] = gradientSensitivities_xy(N, rr, ss, invJacSensitivity)
% [dDx, dDy] = partialDerivativeOperators_xy(N, rr, ss, invJacobianSensitivity)

[Dr, Ds] = support2d.gradients(N, rr, ss);

Dx = Dr*invJacSensitivity(1,1) + Ds*invJacSensitivity(2,1);
Dy = Dr*invJacSensitivity(1,2) + Ds*invJacSensitivity(2,2);
