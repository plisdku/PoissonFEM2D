%% Let's try to assemble a system matrix for a whole grid!!

N = 3;

% Here is a mesh
[vertices,faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

vertices(1,1) = 3;
meshNodes = MeshNodes(faces,vertices,N);

%% 

vv = VVMesh.fv2vv(faces,vertices);

figure(1); clf
VVMesh.plotVV(vv, vertices, 'b-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'o')

%% FEM starts here

fem = PoissonFEM2D(meshNodes);

%% Element matrix sensitivity

% face to test
ff = 1;

% Jacobian element to perturb
ii = 2;
jj = 2;

delta = 1e-6;
jac = meshNodes.getLinearJacobian(ff);
jac2 = jac;
jac2(ii,jj) = jac2(ii,jj) + delta;

% TEST GRADIENT SENSITIVITY

[Dx, Dy, dDx, dDy] = fem.elementGradientMatrix(jac);
[Dx2, Dy2] = fem.elementGradientMatrix(jac2);

dDx_meas = (Dx2-Dx)/delta;
dDy_meas = (Dy2-Dy)/delta;

dDx_calc = dDx{ii,jj};
dDy_calc = dDy{ii,jj};

% TEST QUADRATURE SENSITIVITY

detJ = det(jac);
detJ2 = det(jac2);
dDetJ_meas = (detJ2-detJ)/delta;
dDetJ_calc = detJ*transpose(inv(jac));

[Q, dQdJ] = fem.elementQuadratureMatrix(jac);
Q2 = fem.elementQuadratureMatrix(jac2);
dQ_meas = (Q2-Q)/delta;
dQ_calc = dQdJ{ii,jj};

% TEST POTENTIAL MATRIX SENSITIVITY

[A, dAdJ] = fem.elementPotentialMatrix(jac);
A2 = fem.elementPotentialMatrix(jac2);
dA_meas = (A2-A)/delta;
dA_calc = dAdJ{ii,jj};

% TEST CHARGE MATRIX SENSITIVITY

[B, dBdJ] = fem.elementChargeMatrix(jac);
B2 = fem.elementChargeMatrix(jac2);
dB_meas = (B2-B)/delta;
dB_calc = dBdJ{ii,jj};













%%

[M1, M2, DM1, DM2] = poissonSystemMatrices(meshNodes);

delta = 1e-6;
meshNodes.vertices = meshNodes.vertices + delta*meshNodes.Dvertices;

[M1b, M2b] = poissonSystemMatrices(meshNodes);
meshNodes.vertices = vertices;

DM1_meas = (M1b-M1)/delta;
DM2_meas = (M2b-M2)/delta;

%% Dirichlet boundary conditions

% Separate edges from centers
iEdgeNodes = meshNodes.getBoundaryNodes();
iCenterNodes = meshNodes.getInteriorNodes();

% The forward matrices
M1_center = M1(iCenterNodes, iCenterNodes);
M2_center = M2(iCenterNodes, iCenterNodes);
M1_edges = M1(iCenterNodes, iEdgeNodes);

M1b_center = M1b(iCenterNodes, iCenterNodes);
M2b_center = M2b(iCenterNodes, iCenterNodes);
M1b_edges = M1b(iCenterNodes, iEdgeNodes);

% And for the sensitivity matrices:
DM1_center = DM1(iCenterNodes, iCenterNodes);
DM2_center = DM2(iCenterNodes, iCenterNodes);
DM1_edges = DM1(iCenterNodes, iEdgeNodes);

%% Define a really dumb objective function: a'*u.

objectiveFuncWeights = zeros(numNodes,1);
objectiveFuncWeights(12) = 1; % just pick randomly... should be a center val

%% Solve the system

xy = meshNodes.getNodeCoordinates();
xy_edges = xy(iEdgeNodes,:);

u0 = zeros(numNodes,1);
u0_edges = u0(iEdgeNodes);
u0_edges = cos(3*pi*xy_edges(:,1));

f = zeros(numNodes,1);
f_center = f(iCenterNodes);

% Unperturbed forward system
u_center = M1_center \ (M2_center*f_center - M1_edges*u0_edges);
u = u0;
u(iEdgeNodes) = u0_edges;
u(iCenterNodes) = u_center;

% Objective function!
F = objectiveFuncWeights'*u;

% Perturbed forward system
ub_center = M1b_center \ (M2b_center*f_center - M1b_edges*u0_edges);
ub = u0;
ub(iEdgeNodes) = u0_edges;
ub(iCenterNodes) = ub_center;

Fb = objectiveFuncWeights'*ub;

%% Solve the primal (DD) system

Df_center = 0*f_center;
Du0_edges = 0*u0_edges;

RHS = DM2_center * f_center + M2_center * Df_center ...
    - DM1_center * u_center ...
    - DM1_edges * u0_edges ...
    - M1_edges*Du0_edges;

Du_center = M1_center \ RHS;

Du = 0*u0;
Du(iEdgeNodes) = Du0_edges;
Du(iCenterNodes) = Du_center;

% Objective function
DF = objectiveFuncWeights'*Du;

%% Adjoint system?  Separated edge and centers: barking up the wrong tree?

g_center = objectiveFuncWeights(iCenterNodes);
g_edges = objectiveFuncWeights(iEdgeNodes);

u_adj_center = M1_center' \ g_center;
u_adj_edge = g_edges - M1_edges' * u_adj_center;
u_adj = 0*u;
u_adj(iEdgeNodes) = u_adj_edge;
u_adj(iCenterNodes) = u_adj_center;
%%

% steal RHS from primal again.  Copy-pasted in case I overwrote it in
% the workspace while using cell mode or something.
RHS = DM2_center * f_center + M2_center * Df_center ...
    - DM1_center * u_center ...
    - DM1_edges * u0_edges ...
    - M1_edges*Du0_edges;

DF_dual = u_adj_center' * RHS + ...
    u_adj_edge' * Du0_edges;


%% Interpolate onto a triangular grid

xs = linspace(-1, 5, 100);
ys = linspace(-1, 1, 100);
[xx, yy] = ndgrid(xs, ys);

interpolant = scatteredInterpolant(xy, u, 'linear', 'none');
u_grid = interpolant(xx,yy);

Dinterpolant = scatteredInterpolant(xy, Du, 'linear', 'none');
Du_grid = Dinterpolant(xx,yy);

figure(2); clf
imagesc(xs, ys, u_grid');
hold on
VVMesh.plotVV(vv, vertices, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
title('Linear interpolation')
%%
figure(3); clf
imagesc(xs, ys, Du_grid');
set(gca, 'CLim', [-1,1]*max(abs(get(gca, 'CLim'))));
hold on
VVMesh.plotVV(vv, vertices, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
title('Sensitivity, linear interpolation')
