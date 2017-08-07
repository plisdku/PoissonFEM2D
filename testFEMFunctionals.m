%% Tests!

checkClose = @(a,b) assert(norm(a-b) < 1e-3*norm(a), 'Values are not close');
checkSmall = @(a) assert(norm(a) < 1e-6, 'Value is not small');

%% Setup

N = 4;
vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];

meshNodes = TriNodalMesh(N, faces, vertices);
fem = PoissonFEM2D(meshNodes);
numFaceNodes = meshNodes.basis.numNodes;

%% Test elementIntegralFunctional

% Integrate 1 and get the area of the triangle

f = ones(numFaceNodes,1);
jacobian = meshNodes.getLinearJacobian(1);
[F, ~, ~] = fem.elementIntegralFunctional(f, 0*f, jacobian);
checkClose(F, 0.5);

% Integrate x and get 1/6

xy = meshNodes.getFaceNodeCoordinates(1);
f = xy(:,1);
jacobian = meshNodes.getLinearJacobian(1);
[F, ~, ~] = fem.elementIntegralFunctional(f, 0*f, jacobian);
checkClose(F, 1.0/6);

% Perturb the jacobian and integrate 1; test dFdJ

f = ones(numFaceNodes, 1);
jacobian = meshNodes.getLinearJacobian(1);
delta = 1e-6;
jac2 = jacobian + delta*[1, -2; -3, 4];
[F, dFdJ, ~] = fem.elementIntegralFunctional(f, 0*f, jacobian);
[F2, ~, ~] = fem.elementIntegralFunctional(f, 0*f, jac2);

dF_meas = F2-F;
dF_calc = dot(dFdJ(:), jac2(:)-jacobian(:));
checkClose(dF_meas, dF_calc);

% Integrate a function of u, perturb u, test the sensitivity.
% The function is u(x,y) = y^2.

xy = meshNodes.getFaceNodeCoordinates(1);
u = xy(:,2); % u = y
f = u.^2;
dfdu = 2*u;
delta = 1e-6;
u2 = u + delta;
f2 = u2.^2;
[F, ~, dFdu] = fem.elementIntegralFunctional(f, dfdu, jacobian);
[F2, ~, ~] = fem.elementIntegralFunctional(f2, dfdu, jacobian);

dF_meas = F2-F;
dF_calc = dot(dFdu(:), u2(:)-u(:));
checkClose(dF_meas, dF_calc);


fprintf('Element integral and sensitivity tests PASSED\n');


%% Test surfaceIntegralFunctional

% Integrate 1 and get the triangle area.

u = ones(meshNodes.getNumNodes(), 1);
[F, ~, dFdu] = fem.surfaceIntegralFunctional(@multiplyByOne, [1], u);
checkClose(F, 0.5);

% Perturb u somehow and check the sensitivity.

u2 = u + 1e-3*(1:length(u))';
[F2, ~, ~] = fem.surfaceIntegralFunctional(@multiplyByOne, [1], u2);
dF_meas = F2-F;
dF_calc = dFdu*(u2-u);
checkClose(dF_meas, dF_calc);

fprintf('Mesh integral with sensitivity tests PASSED\n');


%% Test pointEvaluationFunctional

xy = meshNodes.getNodeCoordinates();
u = xy(:,1); % u = x

[F, ~, dFdu] = fem.pointEvaluationFunctional(@multiplyByOne, [0.2,0.2], u);
checkClose(F, 0.2);

u2 = u + 1e-3*(1:length(u))';
[F2, ~, ~] = fem.pointEvaluationFunctional(@multiplyByOne, [0.2,0.2], u2);

dF_meas = F2-F;
dF_calc = dFdu*(u2-u);
checkClose(dF_meas, dF_calc);

fprintf('Mesh point evaluation with sensitivity tests PASSED\n');

