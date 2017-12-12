%% Test a whole adjoint thing

assertClose = @(a,b) assert(norm(a(:)-b(:))/(norm(a(:)) + norm(b(:))) < 1e-5 || norm(a(:) + b(:)) < 1e-50);

%% Set up the mesh

lx = [0, 1, 1, 0];
ly = [0, 0, 1, 1];

in_lx = 0.5 + 0.15*[-1, -1, 1, 1];
in_ly = 0.5 + 0.15*[-1, 1, 1, -1];

density = 4;
[domainV,domainF] = meshPolygon(lx, ly, density, in_lx, in_ly);

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'k-');
patch('Faces', domainF, 'Vertices', domainV, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%% Make an FEM object

% Node orders
N_field = 4;
N_geom = 2;
N_quad = N_field;

lng = LinearNodalGeometry(domainF, domainV, N_geom);
xyNodes = lng.getNodeCoordinates();

tnMesh = TriNodalMesh(domainF, xyNodes, N_field, N_geom, N_quad);
poi = PoissonFEM2D(tnMesh);

dirichletPredicate = @(x,y) norm(x-0.5) < 0.25 || norm(y-0.5) < 0.25;

femp = FEMProblem(poi, dirichletPredicate);


%% Define some shit

x0 = 0.55;
y0 = 0.85;
sigma = 0.05;
freeChargeFunc = @(x,y) exp( (-(x-x0).^2 - (y-y0).^2)/(2*sigma^2));

dirichletFunc = @(x,y) 0; %double(x>0);
neumannFunc = @(x,y) 0;

femp.setSources(freeChargeFunc, dirichletFunc, neumannFunc);


%% Objective function

iArbitraryFace = 13;
arbitraryInteriorNodes = tnMesh.hFieldNodes.getFaceInteriorNodes(iArbitraryFace);
iArbitraryNode = arbitraryInteriorNodes(1);

objFun = @(u) u(iArbitraryNode);
DobjFun = @(u) double( (1:length(u)) == iArbitraryNode );

%% Solve it

femp.solve(objFun);
fprintf('F = %0.4e\n', femp.F);

%femp.solveAdjoint(DobjFun);

%% Plot the field

xCoarse = linspace(0, 1, 40);
yCoarse = linspace(0, 1, 40);

u = femp.poi.tnMesh.rasterizeField(femp.u, xCoarse, yCoarse); % slow as hell
%%
figure(1); clf
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
%colormap orangecrush
hold on
femp.poi.tnMesh.plotMesh('color', 'w');

%% Your potential:

electromagneticPotential = femp.u;



