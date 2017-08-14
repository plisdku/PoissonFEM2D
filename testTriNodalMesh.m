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