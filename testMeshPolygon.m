%% Meshing demos!!!!


lx = [0, 1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 0];
ly = [0, 0, 1, 1, 0, 0, 3, 3, 2, 2, 3, 3];

[domainV,domainF] = meshPolygon(lx, ly, 1);

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'k-');
patch('Faces', domainF, 'Vertices', domainV, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%%

top = MeshTopology(domainF, 4);

%% Make an FEM mesh

N = 8;
%meshNodes = MeshNodes(domainF, domainV, N);
meshNodes = TriNodalMesh(N, domainF, domainV);

xy = meshNodes.getNodeCoordinates();
%%
%iEdgeNodes = meshNodes.getBoundaryNodeCoordinates();
%iCenterNodes = meshNodes.getInteriorNodeCoordinates();
iEdgeNodes = meshNodes.topology.getBoundaryNodes();
iCenterNodes = meshNodes.topology.getInteriorNodes();

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'b-');
hold on
axis xy image
plot(xy(:,1), xy(:,2), 'b.');
plot(xy(iEdgeNodes,1), xy(iEdgeNodes,2), 'ro');
plot(xy(iCenterNodes,1), xy(iCenterNodes,2), 'go');

%% Set up matrices

fem = PoissonFEM2D(meshNodes);
numNodes = meshNodes.topology.getNumNodes();

[A, B] = fem.systemMatrix();

% Separate edges from centers
iEdgeNodes = meshNodes.topology.getBoundaryNodes();
iCenterNodes = meshNodes.topology.getInteriorNodes();

% The forward matrices
M1_center = A(iCenterNodes, iCenterNodes);
M2_center = B(iCenterNodes, iCenterNodes);
M1_edges = A(iCenterNodes, iEdgeNodes);

%% Set up boundary conditions

xy = meshNodes.getNodeCoordinates();
xy_edges = xy(iEdgeNodes,:);

u0 = zeros(numNodes,1);
u0_edges = u0(iEdgeNodes);
%u0_edges = cos(0.25*pi*xy_edges(:,1));
%u0_edges = xy_edges(:,1);

%% Set up free charge

xy_centers = xy(iCenterNodes,:);

f = zeros(numNodes,1);

% Blob of charge at 2.25, 1.25

x0 = 2.25;
y0 = 1.25;
sigma = 0.1;
fFunc = @(x,y) exp( (-(x-x0).^2 - (y-y0).^2)/(2*sigma^2));
f = fFunc(xy(:,1), xy(:,2));


f_center = f(iCenterNodes);

%% Solve!!

% Unperturbed forward system
u_center = M1_center \ (M2_center*f_center - M1_edges*u0_edges);
u = u0;
u(iEdgeNodes) = u0_edges;
u(iCenterNodes) = u_center;

%% Calculate the Ex and Ey fields

[Dx, Dy, count] = meshNodes.getGradientOperators();
Ex = Dx*u;
Ey = Dy*u;

magE = sqrt(abs(Ex).^2 + abs(Ey).^2);

%% Good interpolation

xs = linspace(-1, 4, 400);
ys = linspace(-1, 4, 400);
[xx, yy] = ndgrid(xs, ys);

II = meshNodes.getInterpolationOperator(xx(:), yy(:));

%%

u_grid = reshape(II*magE, size(xx));

%% Interpolate onto a triangular grid

xs = linspace(-1, 4, 400);
ys = linspace(-1, 4, 400);
[xx, yy] = ndgrid(xs, ys);

interpolant = scatteredInterpolant(xy, Ex, 'linear', 'none');
%interpolant = scatteredInterpolant(xy, Dx*xy(:,1), 'linear', 'none');
u_grid = interpolant(xx,yy);

iiIn = inpolygon(xx(:), yy(:), lx(:), ly(:));
u_grid(~iiIn) = 0;

%%
figure(2); clf
imagesc_centered(xs, ys, u_grid');
colormap orangecrush(0.7)
%colorbar
%
hold on
VVMesh.plotFV(domainF, domainV, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
axis xy image vis3d
%title('Linear interpolation')
title('Basis interpolation')
%title('Difference')