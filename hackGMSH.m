%% Write a GMSH file

%% Contour points (oriented correctly) and mesh size parameters

Lx = 6;
Ly = 5;

outerContour = [-Lx, 0; Lx, 0; Lx, Ly; -Lx, Ly];

rAperture1 = 0.5;
ap1 = [-Lx+1, rAperture1; -Lx+2, rAperture1; -Lx+2, Ly-1; -Lx+1, Ly-1];

rAperture2 = 1.5;
ap2 = [-Lx+3, rAperture2; -Lx+4, rAperture2; -Lx+4, Ly-1; -Lx+3, Ly-1];

rAperture3 = 1;
ap3 = [Lx-2, rAperture3; Lx-1, rAperture3; Lx-1, Ly-1; Lx-2, Ly-1];

contours = {outerContour, ap1}; %, ap2, ap3};
meshSizes = [1.0, 0.5, 0.5, 0.5]*4;

%%

writeGEO('fromMatlab.geo', contours, meshSizes);
!/usr/local/bin/gmsh -2 fromMatlab.geo
[mshFaceVertices, mshEdgeVertices, mshVerts, mshEdgeContour, mshEdgeLine] = readMSH('fromMatlab.msh');

%%

figure(101); clf
patch('Faces', mshFaceVertices, 'Vertices', mshVerts, 'FaceColor', 'r');
axis xy image

%%

% Node orders
N_field = 2;
N_geom = 2;
N_quad = N_field;

lng = LinearNodalGeometry(mshFaceVertices, mshVerts, N_geom);
xyNodes = lng.getNodeCoordinates();

tnMesh = TriNodalMesh(mshFaceVertices, xyNodes, N_field, N_geom, N_quad);
poi = PoissonFEM2D(tnMesh);

dirichletPredicate = @(x,y) norm(x-0.5) < 0.25 || norm(y-0.5) < 0.25;

femp = FEMProblem(poi, dirichletPredicate);


%% Questions

% 1. How do I assign different voltages to different contours?
%    Answer: use edgeVertices and edgeContour to find the desired vertices.
%    Find the edges in the TriNodalMesh which have those vertices.
%    Get the nodes corresponding to those edges.
%
% Big question is: which nodes are on which input contours?
% Second question is: how do geom nodes vary with input contour vertices?

iContour = 2;
iContourVertices = mshEdgeVertices(mshEdgeContour == iContour, :); % Mesh vertices on contour
iContourEdges = tnMesh.hMesh.getVertexEdgesExclusive(unique(iContourVertices(:))); % Mesh edges on contour
iContourNodes = tnMesh.hFieldNodes.getEdgeNodes(iContourEdges);

figure(1); clf
tnMesh.plotMesh();
hold on
plot(tnMesh.xyNodes(iEdgeVerts,1), tnMesh.xyNodes(iEdgeVerts, 2), 'ro')
plot(tnMesh.xyNodes(iContourNodes,1), tnMesh.xyNodes(iContourNodes,2), 'gx');

%% Here's how it's going to go, roughly

femp.setDirichlet(iContourNodes, @(x,y) 0.0);
femp.setNeumann(iOtherNodes, @(x,y) 1.0);

dF_dv = dF_dn * dn_dv;


%%

f = FEMInterface();
f.addContour(@(p) [-1, 1, 1, -1], @(p) [-1, -1, 1, 1], @(p) 1, 'dirichlet', @(p,x,y) 0.0);
f.addContour(@(p) 0.5*[-1, -1, 1, 1], @(p) 0.5*[-1, 1, 1, -1], @(p) 0.5, 'dirichlet', @(p,x,y) 1.0);

f.setFreeCharge(@(p,x,y) 0.0);

p = [];
[ff, vv] = f.solve(p, @(u) 1);

%%

xCoarse = linspace(-1, 1, 40);
yCoarse = linspace(-1, 1, 40);

u = f.femProblem.poi.tnMesh.rasterizeField(f.femProblem.u, xCoarse, yCoarse); % slow as hell

%%

figure(1); clf
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
%colormap orangecrush
hold on
f.femProblem.poi.tnMesh.plotMesh('color', 'w');

%%
s.outerBoundary(@(p) [1, 2, 2, 1], @(p) [1, 1, 2, 2], 'Dirichlet', @(p,x,y) 1);
s.innerBoundary(@(p) [0.1, 0.2, 0.2, 0.1], @(p) [0.1, 0.1, 0.2, 0.2], 'Neumann', @(p,x,y) 1);
s.freeCharge(@(p,x,y) 1);

%s.objectiveFunction(objfun, dobjfun); % is this a good way?

s.solve(p, objfun);
s.solveAdjoint(p, dobjfun);

% At this point, I have access to dF/dp, e.g.

do_something_with(s.dF_dp);
plot_something_like(s.dF_dvertices);
plot_something_else(s.dF_dnodes);













