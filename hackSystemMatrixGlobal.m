%% Let's try to assemble a system matrix for a whole grid!!

N = 4;

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

%%

[rr,ss] = support2d.nodes2d(N);
V = support2d.vandermonde(N, rr, ss);
Q = support2d.quadratureKernel(N, rr, ss);

%%
% How can I get the matrices I need?
% I need a jacobian for each triangle.  Ask MeshNodes for that.
% 

numNodes = meshNodes.getNumNodes();
systemMatrix = sparse(numNodes, numNodes);
rhsMatrix = sparse(numNodes, numNodes);

%%

numFaces = size(meshNodes.faces,1);
for ff = 1:numFaces
    jac = meshNodes.getLinearJacobian(ff);
    invJac = inv(jac);
    [Dx,Dy] = support2d.gradients_xy(N, rr, ss, invJac);
    
    Q_xy = Q*det(jac);
    
    M1 = (Dx'*Q_xy*Dx + Dy'*Q_xy*Dy);
    M2 = Q_xy;
    
    iGlobal = meshNodes.local2global(ff);
    
    systemMatrix(iGlobal,iGlobal) = systemMatrix(iGlobal, iGlobal) + M1;
    rhsMatrix(iGlobal, iGlobal) = rhsMatrix(iGlobal, iGlobal) + M2;
    
    doPlot = 0;
    if doPlot
        figure(1); clf
        VVMesh.plotVV(vv, vertices, 'b-');
        hold on
        axis xy image
        plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'o');

        xy = support2d.rs2xy(meshNodes.getFaceVertices(ff)', [rr'; ss']);
        for ii = 1:size(xy,2)
            text(xy(1,ii), xy(2,ii), num2str(iGlobal(ii)), 'FontSize', 24);
        end
        pause(0.1)
    end
end

%% Somehow I need to identify the edge numbers on the outer boundary
% Then I can grab their node indices and create boundary conditions.

iEdgeNodes = meshNodes.getBoundaryNodes();
iCenterNodes = meshNodes.getInteriorNodes();

%%
xy = meshNodes.getNodeCoordinates();

figure(1); clf
VVMesh.plotVV(vv, vertices, 'b-');
hold on
axis xy image
plot(xy(:,1), xy(:,2), 'b.');
plot(xy(iEdgeNodes,1), xy(iEdgeNodes,2), 'ro');
plot(xy(iCenterNodes,1), xy(iCenterNodes,2), 'go');

%% Dirichlet boundary conditions

M1_center = systemMatrix(iCenterNodes, iCenterNodes);
M2_center = rhsMatrix(iCenterNodes, iCenterNodes);

M1_edges = systemMatrix(iCenterNodes, iEdgeNodes);

%% Solve the system

xy_edges = xy(iEdgeNodes,:);

u0 = zeros(numNodes,1);
u0_edges = u0(iEdgeNodes);
%u0_edges(xy_edges(:,1) > 0.5) = 1;
u0_edges = cos(3*pi*xy_edges(:,1));

f = zeros(numNodes,1);
f_center = f(iCenterNodes);

RHS = M2_center * f_center;

u_center = M1_center \ (RHS - M1_edges * u0_edges);
u = u0;
u(iEdgeNodes) = u0_edges;
u(iCenterNodes) = u_center;

%% Interpolate onto a triangular grid

xs = linspace(-1, 5, 100);
ys = linspace(-1, 1, 100);
[xx, yy] = ndgrid(xs, ys);

interpolant = scatteredInterpolant(xy, u, 'linear', 'none');

u_grid = interpolant(xx,yy);


figure(2); clf
imagesc(xs, ys, u_grid');
hold on
VVMesh.plotVV(vv, vertices, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
title('Linear interpolation')

%% Do a proper interpolation

tr = triangulation(meshNodes.faces, meshNodes.vertices);

iFaces = tr.pointLocation(xx(:), yy(:));

for iFace = 1:meshNodes.getNumFaces()
    iFace1 = find(iFaces == iFace);
    xy1 = [xx(iFace1)'; yy(iFace1)'];

    xyTri1 = meshNodes.getFaceVertices(iFace)';

    rs1 = support2d.xy2rs(xyTri1, xy1);

    W = support2d.vandermonde(N, rs1(1,:), rs1(2,:));

    iFace1Nodes = meshNodes.local2global(iFace);
    u_tri1 = W*(V \ u(iFace1Nodes));

    u_grid(iFace1) = u_tri1;
end

figure(4); clf
imagesc(xs, ys, u_grid');
hold on
VVMesh.plotVV(vv, vertices, 'w-')
hold on
plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo')
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 20)
colorbar
title('FEM interpolation');

%%


for ii = 1:length(iFaces)
    iFace = iFaces(ii);
    
    if isnan(iFace)
        continue
    end
    
    xyTri = meshNodes.getFaceVertices(iFace);
    
    rs = support2d.xy2rs(xyTri', [xx(ii); yy(ii)]);
    
end

%%

figure(3); clf
surf(u_grid)


%%

[Dx,Dy] = support2d.partialDerivativeOperators_xy(N, rr, ss);
%[Dr,Ds] = support2d.partialDerivativeOperators(N, rr, ss);

numNodes = size(V,1);

%% First matrix

M1 = (Dx'*Q*Dx + Dy'*Q*Dy);

%% Second matrix

M2 = Q;

%% Edge indices and center indices and such

[iCorners, iEdgeCenters, iCenters] = support2d.classifyNodes(N);

iEdges = setdiff(1:numNodes, iCenters);

% Get the non-boundary part of the system matrix:

M1_center = M1(iCenters,iCenters);
M1_edges = M1(iCenters, iEdges);

%% Solve the system!

% Fixed values
u0 = zeros(numNodes, 1);
u0(iCorners(2)) = 1;
%u0(iEdgeCenters{3}) = 1;
%u0(iEdgeCenters{1}) = -1;
u0_edges = u0(iEdges);
%u0_edges(3:12) = 1;

% RHS: the electric charge density
f = zeros(numNodes,1);
f_center = f(iCenters);
%f_center(36) = 1;
RHS = M2(iCenters, iCenters)*f_center;

% Solve the fucking system

u_center = M1_center \ (RHS - M1_edges*u0_edges);
u = u0;
u(iEdges) = u0_edges;
u(iCenters) = u_center;

figure(1); clf
plot3(rr(iCenters), ss(iCenters), u(iCenters), 'o')