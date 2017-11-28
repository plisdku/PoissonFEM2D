%% Make a test mesh.

vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];

% Node orders
N_field = 4;
N_geom = 3;
N_quad = N_field;

lng = LinearNodalGeometry(faces, vertices, N_geom);
xyNodes = lng.getNodeCoordinates();
tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad);

%% Jacobian sensitivity

delta = 1e-8;

% Two test points
rr = [-0.25, 0.25];
ss = [-0.25, -0.25];

iGlobal = tnMesh.hGeomNodes.getFaceNodes(1); % Global indices of face nodes

J = tnMesh.getJacobianMatrix(1, rr, ss);
DJ = tnMesh.getJacobianSensitivity(rr, ss);

% Iterate over geometry nodes, perturb, test
for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        J2 = tnMesh.perturbed(iGlobal(mm), dirIdx, delta).getJacobianMatrix(1, rr, ss);
        DJ_meas = (J2-J)/delta;
        DJ_calc = DJ(:,:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DJ_calc, DJ_meas);
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end

%% Inverse Jacobian sensitivity

J = tnMesh.getJacobianMatrix(1, rr, ss);
K = tnMesh.getInverseJacobian(1, rr, ss);

for nn = 1:size(J,3)
    invJ = inv(J(:,:,nn));
    
    [same, relErr] = compareNorms(invJ, K(:,:,nn));
    if ~same
        fprintf('K rel err %0.4e ***\n', relErr);
    end
end

DK = tnMesh.getInverseJacobianSensitivity(1, rr, ss);

% Iterate over geometry nodes, perturb, test
for mm = 1:tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        K2 = tnMesh.perturbed(iGlobal(mm), dirIdx, delta).getInverseJacobian(1, rr, ss);
        J2 = tnMesh.perturbed(iGlobal(mm), dirIdx, delta).getJacobianMatrix(1, rr, ss);
        DK_meas = (K2-K)/delta;
        DK_calc = DK(:,:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DK_calc, DK_meas);
        if ~same
            fprintf('DK error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end


%% Edge Jacobian sensitivity

iGlobal = tnMesh.hGeomNodes.getEdgeNodes(1); % Global indices of face nodes

J = tnMesh.getEdgeJacobianMatrix(1, rr, -1);
DJ = tnMesh.getEdgeJacobianSensitivity(rr, -1);

% Iterate over geometry nodes, perturb, test
for mm = 1:tnMesh.hGeomNodes.N
    for dirIdx = 1:2
        J2 = tnMesh.perturbed(iGlobal(mm), dirIdx, delta).getEdgeJacobianMatrix(1, rr, -1);
        DJ_meas = (J2-J)/delta;
        DJ_calc = DJ(:,:,dirIdx,mm);
        
        [same, relErr, normDiff, normExact] = compareNorms(DJ_calc, DJ_meas);
        if ~same
            fprintf('DJ error %0.4e out of %0.4e ***\n', normDiff, normExact)
        end
    end
end

%% Face Jacobian determinant sensitivity







%% Let's try to assemble a system matrix for a whole grid!!

N = 4;

% Here is a mesh
[vertices,faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

meshNodes = TriNodalMesh(N,faces,vertices);

numNodes = meshNodes.hNodes.getNumNodes();

%% 

vv = VVMesh.fv2vv(faces,vertices);

figure(1); clf
VVMesh.plotVV(vv, vertices, 'b-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'o')

%% Jacobian sensitivity to vertex perturbations

ff = 1;
ii = 1; jj = 1;
delta = 1e-6;

% perturb this vertex
iVertInFace = 3;
tmp = meshNodes.hMesh.getFaceVertices(ff);
iVert = tmp(iVertInFace);
iXY = 2;

vertices2 = vertices;
vertices2(iVert, iXY) = vertices2(iVert, iXY) + delta;
meshNodes2 = TriNodalMesh(N, faces, vertices2);

jac = meshNodes.getLinearJacobian(ff);
jac2 = meshNodes2.getLinearJacobian(ff);
Djac = meshNodes.getLinearJacobianSensitivity(ff);

Djac_meas = (jac2-jac)/delta;
Djac_calc = Djac(:,:, iVertInFace, iXY);

fprintf('dJdv error: %g\n', norm(Djac_meas - Djac_calc));

%% FEM starts here

fem = PoissonFEM2D(meshNodes);

%% Element matrix sensitivity

% face to test
ff = 1;

% Jacobian element to perturb
ii = 2;
jj = 1;

fprintf('Testing element sensitivities to perturbing J(%d,%d)\n', ii, jj);

delta = 1e-6;
jac = meshNodes.getLinearJacobian(ff);
jac2 = jac;
jac2(ii,jj) = jac2(ii,jj) + delta;

% TEST GRADIENT SENSITIVITY

[Dx, Dy, dDx, dDy] = fem.elementGradientMatrix(jac);
[Dx2, Dy2] = fem.elementGradientMatrix(jac2);

dDx_meas = (Dx2-Dx)/delta;
dDy_meas = (Dy2-Dy)/delta;

dDx_calc = dDx(:,:,ii,jj);
dDy_calc = dDy(:,:,ii,jj);

fprintf('dDx error: %g\n', norm(dDx_meas - dDx_calc));
fprintf('dDy error: %g\n', norm(dDy_meas - dDy_calc));

% TEST QUADRATURE SENSITIVITY

detJ = det(jac);
detJ2 = det(jac2);
dDetJ_meas = (detJ2-detJ)/delta;
dDetJ_calc = detJ*transpose(inv(jac));

[Q, dQdJ] = fem.elementQuadratureMatrix(jac);
Q2 = fem.elementQuadratureMatrix(jac2);
dQ_meas = (Q2-Q)/delta;
dQ_calc = dQdJ(:,:,ii,jj);

fprintf('dQ error: %g\n', norm(dQ_meas - dQ_calc));

% TEST POTENTIAL MATRIX SENSITIVITY

[A, dAdJ] = fem.elementPotentialMatrix(jac);
A2 = fem.elementPotentialMatrix(jac2);
dA_meas = (A2-A)/delta;
dA_calc = dAdJ(:,:,ii,jj);

fprintf('dA error: %g\n', norm(dA_meas - dA_calc));

% TEST CHARGE MATRIX SENSITIVITY

[B, dBdJ] = fem.elementChargeMatrix(jac);
B2 = fem.elementChargeMatrix(jac2);
dB_meas = (B2-B)/delta;
dB_calc = dBdJ(:,:,ii,jj);

fprintf('dB error: %g\n', norm(dB_meas - dB_calc));

% TEST NEUMANN MATRIX SENSITIVITY
[C, dCdJ] = fem.elementNeumannMatrix(jac);
C2 = fem.elementNeumannMatrix(jac2);
dC_meas = (C2-C)/delta;
dC_calc = dCdJ(:,:,ii,jj);

fprintf('dC error: %g\n', norm(dC_meas - dC_calc));

%% Perturb vertices instead...

fprintf('Element matrices\n');

ff = 3;
iVertInFace = 3; % local vertex to perturb
tmp = meshNodes.hMesh.getFaceVertices(ff);
vv = tmp(iVertInFace);
xy = 2; % direction to perturb
delta = 1e-6;

vertices2 = vertices;
vertices2(vv, xy) = vertices2(vv, xy) + delta;
meshNodes2 = TriNodalMesh(N, faces, vertices2);

jacobian = meshNodes.getLinearJacobian(ff);
jacobian2 = meshNodes2.getLinearJacobian(ff);

dJdv = meshNodes.getLinearJacobianSensitivity(ff);

% Test dJdv
dJdv_meas = (jacobian2-jacobian)/delta;
dJdv_calc = dJdv(:,:,iVertInFace, xy);
fprintf('dJdv error: %g\n', norm(dJdv_meas - dJdv_calc));

% Sensitivity of system matrix to perturbation
% Two parts: A and B (potential and charge)

[A, dAdJ] = fem.elementPotentialMatrix(jacobian);
A2 = fem.elementPotentialMatrix(jacobian2);

dAdv = multiplyTensors.txt(dAdJ, 4, dJdv, 4, 3:4, 1:2);

dAdv_calc = dAdv(:,:,iVertInFace, xy);
dAdv_meas = (A2-A)/delta;
fprintf('dAdv error: %g\n', norm(dAdv_calc - dAdv_meas));

[B, dBdJ] = fem.elementChargeMatrix(jacobian);
B2 = fem.elementChargeMatrix(jacobian2);
dBdv = multiplyTensors.txt(dBdJ, 4, dJdv, 4, 3:4, 1:2);
dBdv_calc = dBdv(:,:,iVertInFace, xy);
dBdv_meas = (B2-B)/delta;
fprintf('dBdv error: %g\n', norm(dBdv_calc - dBdv_meas));

% Need to do this with a 1D Jacobian.
%[C, dCdJ] = fem.elementNeumannMatrix(jacobian);
%C2 = fem.elementNeumannMatrix(jacobian2);
%dCdv = multiplyTensors.txt(dCdJ, 3, dJdv, 3, 3, 1);
%dCdv_calc = dCdv(:,:,iVertInFace, xy);
%dCdv_meas = (C2-C)/delta;
%fprintf('dCdv error: %g\n', norm(dCdv_calc - dCdv_meas));

%% System matrices and sensitivities

fprintf('System matrices\n')

numVerts = meshNodes.hMesh.getNumVertices();

iVert = 1;
iXY = 1;
%for iVert = 1:numVerts

delta = 1e-3;
vertices2 = vertices;
vertices2(iVert, iXY) = vertices2(iVert, iXY) + delta;

meshNodes = TriNodalMesh(N, faces, vertices);
meshNodes2 = TriNodalMesh(N, faces, vertices2);

fem = PoissonFEM2D(meshNodes);
fem2 = PoissonFEM2D(meshNodes2);

[A, dAdv] = fem.systemMatrix();
[B, dBdv] = fem.rhsMatrix();
[C, dCdv] = fem.neumannMatrix();
A2 = fem2.systemMatrix();
B2 = fem2.rhsMatrix();
C2 = fem2.neumannMatrix();

dAdv_meas = (A2-A)/delta;
dBdv_meas = (B2-B)/delta;
dCdv_meas = (C2-C)/delta;

dAdv_calc = dAdv{iVert,iXY};
dBdv_calc = dBdv{iVert,iXY};
dCdv_calc = dCdv{iVert,iXY};

fprintf('dAdv error %g / %g\n', norm(full(dAdv_calc - dAdv_meas)), norm(full(dAdv_meas)));
fprintf('dBdv error %g / %g\n', norm(full(dBdv_calc - dBdv_meas)), norm(full(dBdv_meas)));
fprintf('dCdv error %g / %g\n', norm(full(dCdv_calc - dCdv_meas)), norm(full(dCdv_meas)));
%end


%% Dirichlet boundary conditions

% Separate edges from centers
iEdgeNodes = meshNodes.hNodes.getBoundaryNodes();
iCenterNodes = meshNodes.hNodes.getInteriorNodes();

% The forward matrices
M1_center = A(iCenterNodes, iCenterNodes);
M2_center = B(iCenterNodes, iCenterNodes);
M1_edges = A(iCenterNodes, iEdgeNodes);

M1b_center = A2(iCenterNodes, iCenterNodes);
M2b_center = B2(iCenterNodes, iCenterNodes);
M1b_edges = A2(iCenterNodes, iEdgeNodes);

% And for the sensitivity matrices:
DM1_center = dAdv_calc(iCenterNodes, iCenterNodes);
DM2_center = dBdv_calc(iCenterNodes, iCenterNodes);
DM1_edges = dAdv_calc(iCenterNodes, iEdgeNodes);

%% Define a really dumb objective function: a'*u.

iObjectiveFnNode = 6;
objectiveFuncWeights = zeros(numNodes,1);
objectiveFuncWeights(iObjectiveFnNode) = 1; % just pick randomly... should be a center val

doPlot = 1;
if doPlot
    vv = VVMesh.fv2vv(faces,vertices);
    
    nodeXY = meshNodes.getNodeCoordinates();
    
    figure(1); clf
    VVMesh.plotVV(vv, vertices, 'b-')
    ax = axis;
    hold on
    VVMesh.plotVV(vv, vertices + 10*(vertices2-vertices), 'b--')
    plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'o')
    plot(nodeXY(:,1), nodeXY(:,2), 'r.', 'MarkerSize', 20);
    plot(nodeXY(iObjectiveFnNode,1), nodeXY(iObjectiveFnNode,2), 'ko', 'MarkerSize', 10);
    axis(ax);
end

%% Solve the system

xy = meshNodes.getNodeCoordinates();
xy2 = meshNodes2.getNodeCoordinates();
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

fprintf('Measured gradient of F is %.15g\n', (Fb-F)/delta);

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

fprintf('Primal gradient of F is %.15g\n', DF);

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

fprintf('Dual gradient of F is %0.15g\n', DF_dual);

%% Interpolate onto a triangular grid

xs = linspace(-1, 5, 100);
ys = linspace(-1, 1, 100);
[xx, yy] = ndgrid(xs, ys);

interpolant = scatteredInterpolant(xy, u, 'linear', 'none');
u_grid = interpolant(xx,yy);
interpolant_b = scatteredInterpolant(xy2, ub, 'linear', 'none');
u_grid_b = interpolant_b(xx,yy);

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
axis xy
title('Linear interpolation')

figure(3); clf
imagesc(xs, ys, u_grid_b');
hold on
VVMesh.plotVV(vv, vertices2, 'w-')
hold on
plot(vertices2(:,1), vertices2(:,2), 'wo')
plot(xy2(:,1), xy2(:,2), 'w.', 'MarkerSize', 20)
colorbar
axis xy
title('Linear interpolation (2)')
%%
figure(4); clf
imagesc(xs, ys, Du_grid');
set(gca, 'CLim', [-1,1]*max(abs(get(gca, 'CLim'))));
hold on
VVMesh.plotVV(vv, vertices, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
title('Sensitivity, linear interpolation')
