function [Dr, Ds] = gradients(N, rr, ss)
%[Dr, Ds] = gradients(N, rr, ss)

warning('DEPRECATED.  Use BasisNodes.');

import support2d.*

[V, dVdr, dVds] = vandermonde(N, rr, ss);

Dr = dVdr * inv(V);
Ds = dVds * inv(V);