%% Test Jacobians

relErr = @(a,b) norm(a-b)/norm(a);

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

for orientation = [1,-1]
for iEdge = 1:3

    r0 = 0.1;
    xy0 = tnMesh.getEdgeCoordinates(iEdge, r0, orientation);
    xy1 = tnMesh.getEdgeCoordinates(iEdge, r0+delta, orientation);

    dxy_dr_meas = (xy1 - xy0)/delta;
    dxy_dr = tnMesh.getEdgeJacobian(iEdge, r0, orientation);

    assert(norm(dxy_dr_meas - dxy_dr) < 1e-6);
end  
end

fprintf('Single-point 1D Jacobian test PASSED\n');

%% Test differentiation

N_field = 14;
N_geom = 3;

vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];
[xyNodes, faces] = nodalMesh(faces, vertices, N_geom);
xyNodes(6,2) = xyNodes(6,2) + 0.1;
xyNodes(5,1) = xyNodes(5,1) - 0.1;

tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_field);
[Dx, Dy] = tnMesh.getGradientOperators();
xy = tnMesh.getNodeCoordinates(); % Field nodes

figure(1); clf; tnMesh.plotMesh(); hold on; plot(xy(:,1), xy(:,2), 'bo');

% Test 1: linear polynomial
%f = @(xy) xy(:,1) - xy(:,2);
%df_dx = @(xy) ones(size(xy(:,1)));
%df_dy = @(xy) -ones(size(xy(:,1)));

% Test 2: bilinear
%f = @(xy) xy(:,1).*xy(:,2);
%df_dx = @(xy) xy(:,2);
%df_dy = @(xy) xy(:,1);

% Test 3: quadratic
%f = @(xy) xy(:,1).^2 - xy(:,2).^2;
%df_dx = @(xy) 2*xy(:,1);
%df_dy = @(xy) -2*xy(:,2);

% Test 4: exponential
a = 0.5;
f = @(xy) exp(a*xy(:,1) - a*xy(:,2));
df_dx = @(xy) a*exp(a*xy(:,1) - a*xy(:,2));
df_dy = @(xy) -a*exp(a*xy(:,1) - a*xy(:,2));

zz = f(xy);

df_dx_meas = Dx*zz;
df_dx_calc = df_dx(xy);

df_dy_meas = Dy*zz;
df_dy_calc = df_dy(xy);

fprintf('Relative error of df/dx = %0.4e\n', relErr(df_dx_meas, df_dx_calc));
fprintf('Relative error of df/dy = %0.4e\n', relErr(df_dy_meas, df_dy_calc));

%assert(norm(df_dx_meas - df_dx_calc) < 1e-6);
%assert(norm(df_dy_meas - df_dy_calc) < 1e-6);

%% Interpolation operator


N_field = 4;
N_geom = 5;

vertices = [0,0; 1,0; 0,1];
faces = [1,2,3];

%[xyNodes, faces] = nodalMesh(faces, vertices, N_geom);
[xyNodes, faces] = nodalWagonWheel(5, N_geom);
xyNodes = xyNodes + 0.03*randn(size(xyNodes));

%xyNodes(6,2) = xyNodes(6,2) + 0.1;
%xyNodes(5,1) = xyNodes(5,1) - 0.1;
%xyNodes(11,1) = xyNodes(11,1) - 0.3;

tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_field);

figure(3); clf
tnMesh.plotMesh();

%%

% Test functions
%f = @(xy) xy(:,1).^2 - xy(:,2).^2;
f = @(xy) sin(2*xy(:,1)) + cos(3*xy(:,2));

tnMesh.plotMesh()

%xs = linspace(0.2, 0.4);
%ys = linspace(0.2, 0.4);
%xs = linspace(-0.2, 1.2);
%ys = linspace(-0.2, 1.2);

xs = linspace(-1.2, 1.2, 400);
ys = linspace(-1.2, 1.2, 400);

[xx,yy] = ndgrid(xs,ys);
xyDense = [xx(:)'; yy(:)'];

xyFields = tnMesh.getNodeCoordinates();

M = tnMesh.getInterpolationOperator(xx(:), yy(:));

zzNodes = f(xyFields);
zz = f(xyDense');
zzInterp = M*zzNodes;

%%

zzErr = zzInterp - zz;
zzErr(zzInterp == 0) = NaN;

%%

figure(1); clf
imagesc(xs, ys, reshape(zz, size(xx))', [-1.5, 1.5]);
axis xy image vis3d
hold on
tnMesh.plotMesh('Color', 'w')
plot(xyFields(:,1), xyFields(:,2), 'wo')
title('Function to interpolate')
colorbar

%%

figure(2); clf
imagesc(xs, ys, reshape(zzInterp, size(xx))', [-1.5, 1.5]);
%colormap cool
axis xy image vis3d
hold on
tnMesh.plotMesh('Color', 'w')
plot(xyFields(:,1), xyFields(:,2), 'wo')
title('Interpolated')
colorbar


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