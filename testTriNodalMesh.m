%% Test Jacobians

N_field = 3;
N_geom = 3;
[xyNodes, faces] = nodalWagonWheel(5, N_geom);

%% Test Jacobian on unit simplex
% Try changing N and the vertices, or moving a node.

N = 3;
vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];
[xyNodes, faces] = nodalMesh(faces, vertices, N);
%xyNodes(6,2) = xyNodes(6,2) + 0.1;

tnMesh = TriNodalMesh(faces, xyNodes, N, N, N);

%figure(1); clf; tnMesh.plotMesh();

% To test the Jacobian... pick a point I guess?
% Evaluate the coordinate transformation in one triangle
% at an initial point (r0,s0) and then at perturbed points (r0+delta,s0)
% and (r0,s0+delta).  Calculate the Jacobian from these, and compare to
% the analytical value.

delta = 1e-6;
r0 = -0.3;
s0 = -0.2;

xy0 = tnMesh.getFaceCoordinates(1, r0, s0);
xy1_r = tnMesh.getFaceCoordinates(1, r0+delta, s0);
xy1_s = tnMesh.getFaceCoordinates(1, r0, s0+delta);

dxy_dr_meas = (xy1_r - xy0)/delta;
dxy_ds_meas = (xy1_s - xy0)/delta;

[dxy_dr, dxy_ds] = tnMesh.getJacobian(1, r0, s0);

assert(norm(dxy_dr - dxy_dr_meas) < 1e-6);
assert(norm(dxy_ds - dxy_ds_meas) < 1e-6);

fprintf('Single-point Jacobian test PASSED\n');

%% Test 1D Jacobian

delta = 1e-6;

for iEdge = 1:3

    r0 = 0.1;
    xy0 = tnMesh.getEdgeCoordinates(iEdge, r0);
    xy1 = tnMesh.getEdgeCoordinates(iEdge, r0+delta);

    dxy_dr_meas = (xy1 - xy0)/delta;
    dxy_dr = tnMesh.getEdgeJacobian(iEdge, r0);

    assert(norm(dxy_dr_meas - dxy_dr) < 1e-6);
    
end

fprintf('Single-point 1D Jacobian test PASSED\n');

%% Interpolation operator sensitivity

[vertices,faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

%vertices = [0,0; 1,0; 0,1; 1,0.95];
%faces = [1,2,4];

N = 4;
tnMesh = TriNodalMesh(N, faces, vertices);

figure(1); clf
VVMesh.plotFV(faces, vertices, 'k');

%% Make the interpolation operator and its derivative

xs = linspace(-0.9,0.9, 10);
ys = xs;
[xx,yy] = ndgrid(xs,ys);

%xx = 0.2;
%yy = 0.25;

M = tnMesh.getInterpolationOperator(xx(:), yy(:));
DM = tnMesh.getInterpolationOperatorSensitivity(xx(:), yy(:));

%% Test perturbations

% Try a few vertices and directions
iVert = 2;
iXY = 2;

delta = 1e-6;
vertices2 = vertices;
vertices2(iVert, iXY) = vertices2(iVert, iXY) + delta;

tnMesh2 = TriNodalMesh(N, faces, vertices2);
M2 = tnMesh2.getInterpolationOperator(xx(:), yy(:));

DM_meas = (M2-M)/delta;
DM_calc = DM{iVert,iXY};

fprintf('Norm of DM_measured: %g\n', norm(full(DM_meas)));
fprintf('Norm of error: %e\n', norm(full(DM_meas - DM_calc)));