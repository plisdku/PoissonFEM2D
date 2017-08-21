function Dr = partialDerivativeOperator(N, rr)
% Dr = partialDerivativeOperator(N, rr)


V = support.vandermonde(N, rr);
dVdr = support.gradVandermonde(N, rr);

Dr = dVdr * inv(V);