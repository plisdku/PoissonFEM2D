%% Meshing demos!!!!

lx = [0, 1, 1, 0];
ly = [0, 0, 1, 1];

in_lx = 0.5 + 0.05*[-1, -1, 1, 1];
in_ly = 0.5 + 0.05*[-1, 1, 1, -1];

density = 4;
[domainV,domainF] = meshPolygon(lx, ly, density, in_lx, in_ly);

%% The input to my FEM monster:
% domainV
% domainF
% choice of functional
% choice of boundary conditions
% choice of source terms

%%
figure(1); clf
VVMesh.plotFV(domainF, domainV, 'k-');
patch('Faces', domainF, 'Vertices', domainV, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%% Make an FEM mesh

N = 8;
%meshNodes = MeshNodes(domainF, domainV, N);
meshNodes = TriNodalMesh(N, domainF, domainV);

xy = meshNodes.getNodeCoordinates();
%%
%iEdgeNodes = meshNodes.getBoundaryNodeCoordinates();
%iCenterNodes = meshNodes.getInteriorNodeCoordinates();
iEdgeNodes = meshNodes.getBoundaryNodes();
iCenterNodes = meshNodes.getInteriorNodes();

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'b-');
hold on
axis xy image
%meshNodes.plotEdgeIndices();
%meshNodes.plotNodeIndices();
plot(xy(:,1), xy(:,2), 'b.');
plot(xy(iEdgeNodes,1), xy(iEdgeNodes,2), 'ro');
plot(xy(iCenterNodes,1), xy(iCenterNodes,2), 'go');

%% Set up matrices

fem = PoissonFEM2D(meshNodes);
numNodes = meshNodes.getNumNodes();

[A, dA_dv] = fem.systemMatrix();
[B, dB_dv] = fem.rhsMatrix();
[NM, dNM_dv] = fem.neumannMatrix();

% Separate edges from centers
%iEdgeNodes = meshNodes.getBoundaryNodes();
%iCenterNodes = meshNodes.getInteriorNodes();
%%
% Choose the Neumann edges somehow

% General idea: we have selected nodes that are Neumann and nodes that
% are Dirichlet.  We will make matrices... accordingly... :-O
%
% Key tasks: choosing edges, then finding nodes for those edges.

boundaryEdges = meshNodes.getBoundaryEdges();

predicate = @(v0,v1) norm(v0-0.5) < 0.25;

chooseMe = false(length(boundaryEdges), 1);
for ii = 1:length(boundaryEdges)
    iEdge = boundaryEdges(ii);
    verts = meshNodes.getVertexNodeCoordinates(meshNodes.getEdgeVertices(iEdge));
    
    chooseMe(ii) = predicate(verts(1,:), verts(2,:));
end

innerEdges = boundaryEdges(chooseMe);
outerEdges = setdiff(boundaryEdges, innerEdges);

innerNodes = meshNodes.getEdgeNodes(innerEdges);
outerNodes = meshNodes.getEdgeNodes(outerEdges);

%iNeumann = outerNodes;
%iDirichlet = [innerNodes, outerNodes];
iDirichlet = innerNodes;
iNeumann = setdiff(meshNodes.getBoundaryNodes(), iDirichlet);

%% Set up free charge

f = zeros(numNodes,1);
x0 = 0.25;
y0 = 0.5;
sigma = 0.05;
fFunc = @(x,y) exp( (-(x-x0).^2 - (y-y0).^2)/(2*sigma^2));
f = fFunc(xy(:,1), xy(:,2));

%% NEW WAY: Dirichlet and Neumann are both Robin boundary conditions.

if 0
    % Dirichlet: k is big, RHS = k*u0
    % Neumann: k is small, RHS = e_n

    % Put Dirichlet on the ... inner nodes?

    k_Dirichlet = 1e1*ones(numNodes, 1);
    k_Dirichlet(iNeumann) = 0;

    k_Neumann = 1e-1*ones(numNodes, 1);
    k_Neumann(iDirichlet) = 0;

    neumannMatrix = fem.neumannMatrix(k_Neumann);
    dirichletMatrix = fem.neumannMatrix(k_Dirichlet);

    lhsRobinMatrix = neumannMatrix + dirichletMatrix;

    rhsDirichletMatrix = dirichletMatrix(:,iDirichlet);
    u0_dirichlet = zeros(length(iDirichlet),1) + 4;
    u0_dirichlet(1:length(innerNodes)) = 0;
    
    u = (A - lhsRobinMatrix) \ (B*f - rhsDirichletMatrix*u0_dirichlet);
end

%% Old-style solve!!

iCenterNodes = setdiff(1:numNodes, iDirichlet);

A_center = A(iCenterNodes, iCenterNodes);
B_center = B(iCenterNodes, iCenterNodes);
NM_center = NM(iCenterNodes, iCenterNodes);
NM_neumann = NM(iCenterNodes, iNeumann);

A_dirichlet = A(iCenterNodes, iDirichlet);

f_center = f(iCenterNodes);
u0_dirichlet = zeros(length(iDirichlet),1); % + randn(length(iDirichlet),1);
%u0_dirichlet(1:length(innerNodes)) = 0;
en_neumann = 0*iNeumann;

if numel(iNeumann) > 0
    u_center = A_center \ (B_center*f_center - A_dirichlet*u0_dirichlet - NM_neumann*en_neumann);
else
    u_center = A_center \ (B_center*f_center - A_dirichlet*u0_dirichlet);
end
u = zeros(numNodes,1);
u(iCenterNodes) = u_center;
u(iDirichlet) = u0_dirichlet;

%% Evaluate a functional!!

% Point evaluation of u.
xyMeas = [0.2; 0.2];

[F, dFdv, dFdu] = fem.pointEvaluationFunctional(@multiplyByOne, xyMeas, u);

%%

%meshNodes.getInterpolationOperator();
dIdv = meshNodes.getInterpolationOperatorSensitivity([0.2], [0.3]);

%% Solve the adjoint system!

dFdu_center = dFdu(iCenterNodes);

v = A_center' \ dFdu_center';

% just for making a plot if you want to
v_big = 0*u;
v_big(iCenterNodes) = v;

%% Get some fucking sensitivities!!!!

% Sensitivity with respect to Dirichlet boundary values
dF_dud = -v'*A_dirichlet;
dF_dud_big = 0*u;
dF_dud_big(iDirichlet) = dF_dud;

% Sensitivity with respect to free charge
dF_df = v'*B_center;
dF_df_big = 0*u;
dF_df_big(iCenterNodes) = dF_df;

% Sensitivity with respect to normal E on Neumann boundary nodes
dF_den = -v'*NM_neumann;
dF_den_big = 0*u;
dF_den_big(iNeumann) = dF_den;

%% Sensitivity to geometry changes?

numVerts = size(domainV, 1);

dFdv_total = 0*domainV;

for vv = 1:numVerts
    
    wx = -dA_dv{vv,1}(iCenterNodes, iCenterNodes)*u_center ...
        - dNM_dv{vv,1}(iCenterNodes, iNeumann)*en_neumann ...
        - dA_dv{vv,1}(iCenterNodes, iDirichlet)*u0_dirichlet ...
        - dB_dv{vv,1}(iCenterNodes, iCenterNodes)*f_center;
    
    wy = -dA_dv{vv,2}(iCenterNodes, iCenterNodes)*u_center ...
        - dNM_dv{vv,2}(iCenterNodes, iNeumann)*en_neumann ...
        - dA_dv{vv,2}(iCenterNodes, iDirichlet)*u0_dirichlet ...
        - dB_dv{vv,2}(iCenterNodes, iCenterNodes)*f_center;
    
    dFdvx = v'*wx;
    dFdvy = v'*wy;
    
    dFdv_total(vv,1) = dFdvx;
    dFdv_total(vv,2) = dFdvy;
end

dFdv_total = dFdv_total + dFdv;

%%

figure(1); clf
VVMesh.plotFV(domainF, domainV, 'k-');
patch('Faces', domainF, 'Vertices', domainV, 'FaceColor', 'g', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
hold on
plot(domainV(:,1), domainV(:,2), 'o')
quiver(domainV(:,1), domainV(:,2), dFdv_total(:,1), dFdv_total(:,2), 'linewidth', 2);
%quiver(domainV(:,1), domainV(:,2), dFdv(:,1), dFdv(:,2), 'linewidth', 2);
axis xy image
title('Vertex sensitivities??')
plot(xyMeas(1), xyMeas(2), 'rx', 'MarkerSize', 10)


%% Calculate the Ex and Ey fields

[Dx, Dy, count] = meshNodes.getGradientOperators();
Ex = Dx*u;
Ey = Dy*u;

magE = sqrt(abs(Ex).^2 + abs(Ey).^2);

%% Good interpolation

%xs = linspace(-1, 4, 400);
%ys = linspace(-1, 4, 400);
xs = linspace(-0.1, 1.1, 400);
ys = linspace(-0.1, 1.1, 400);
[xx, yy] = ndgrid(xs, ys);

II = meshNodes.getInterpolationOperator(xx(:), yy(:));

%%

u_grid = reshape(II*u, size(xx));

%%

interpolant = scatteredInterpolant(xy, dF_df_big, 'linear', 'none');
u_grid = interpolant(xx,yy);
u_grid(isnan(u_grid)) = 0;

%%
figure(2); clf
imagesc_centered(xs, ys, u_grid'); %, [0, 0.1]);
colormap orangecrush(0.7)
%colorbar
%
hold on
%VVMesh.plotFV(domainF, domainV, 'w-', 'linewidth', 0.01)
hold on
%plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo');
plot(xy(:,1), xy(:,2), 'w.', 'MarkerSize', 2)
colorbar
axis xy image vis3d
title('Basis interpolation')

%%

plot(domainV(:,1), domainV(:,2), 'wo')
hold on
quiver(domainV(:,1), domainV(:,2), dFdv_total(:,1), dFdv_total(:,2), 'r');
