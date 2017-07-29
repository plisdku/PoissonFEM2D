%% Tests of functionals
%
% Test in one triangle

N = 4;
vertices = [0, 0; 2, 0; 0, 1];
faces = [1, 2, 3];
meshNodes = MeshNodes(faces, vertices, N);

x_eval = [0.4; 0.5];

xy_tri = vertices';
r_eval = support2d.xy2rs(xy_tri, x_eval);
phi_eval = support2d.vandermonde(N, r_eval(1), r_eval(2));

[rr,ss] = support2d.nodes2d(N);
V = support2d.vandermonde(N, rr, ss);

%%
xn = meshNodes.getNodeCoordinates();

figure(1); clf
plot(xn(:,1), xn(:,2), 'o')
hold on
VVMesh.plotFV(faces, vertices, 'b-');
plot(x_eval(1), x_eval(2), 'x')
axis equal
grid on
title('Point evaluation test')
%% Interpolation

% Fixed function.
f_func = @(xy) xy(1,:) - xy(2,:);
f = f_func(xn')';

V = support2d.vandermonde(N);

rs_eval = support2d.xy2rs(xy_tri, x_eval);

%% Evaluate on a grid

xs = linspace(0,2);
ys = linspace(0,1);
[xx,yy] = ndgrid(xs,ys);
xs_grid = [xx(:)'; yy(:)'];
rs_grid = support2d.xy2rs(xy_tri, [xx(:)'; yy(:)']);

f_exact = reshape(f_func(xs_grid), size(xx));

figure(2); clf
imagesc(xs, ys, f_exact'); axis xy image
colorbar

%% V on the grid

V_grid = support2d.vandermonde(N, rs_grid(1,:), rs_grid(2,:));

%% Modal coefficients of f

[rr,ss] = support2d.nodes2d(N);
xyn = support2d.rs2xy(xy_tri, [rr, ss]');

f_modal = V \ f_func(xyn)';

f_interp = reshape(V_grid * f_modal, size(xx));


%%

f_eval = support2d.interpolate(f, x_eval, xy_tri, V);