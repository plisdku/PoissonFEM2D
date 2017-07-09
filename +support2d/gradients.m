function [Dr, Ds] = partialDerivativeOperators(N, rr, ss)

import support2d.*

[V, dVdr, dVds] = vandermonde(N, rr, ss);

Dr = dVdr * inv(V);
Ds = dVds * inv(V);