function Dr = partialDerivativeOperator(N, rr)
% Dr = partialDerivativeOperator(N, rr)


V = support.vandermonde(rr, N);
dVdr = support.gradVandermonde(rr, N);

Dr = dVdr * inv(V);