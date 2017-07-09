%% Integration:

N = 20;
rr = support.nodes1d(N);
Q = support.quadratureKernel(N, rr);

%% Area of the single element:

uno = ones(N,1);

elementArea = uno'*Q*uno; % should be 2

fprintf('Element area %g should be 2\n', elementArea);

%% Sin squared

f = sin(rr'*pi);

integralVal = f'*Q*f;
integralExpected = 1;

fprintf('Sin squared integral %g should be %g\n', integralVal, integralExpected);

%% Make kernel for integrating Du * Dv for functions u and v

Dr = support.partialDerivativeOperator(N, rr);

K = Dr' * Q * Dr;

% Simplest: u = x, v = x

uu = rr';
vv = rr';

intVal = uu'*K*vv;

expected = 2;

fprintf('Integral is %g, expected %g\n', intVal, expected);

%% 