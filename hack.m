% To test the curvilinear Jacobian, I will want to be able to
% - plot a curvilinear element
% - carry out the mappings X(r,s) and Y(r,s)
% - verify the Jacobian numerically

%faces = [1 2 3; 3 2 4];
%vertices = [0,0; 1,0; 0,1; 1,1];

[vertices, faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

% Node orders
N_field = 4;
N_geom = 3;
N_quad = N_field;

lng = LinearNodalGeometry(faces, vertices, N_geom);
xyNodes = lng.getNodeCoordinates();
iBdy = lng.hNodes.getBoundaryNodes();

theta = atan2(xyNodes(:,2), xyNodes(:,1));

xyNodes(iBdy,1) = cos(theta(iBdy)); % .* (1 + 0.1*sin(8*theta(iBdy)));
xyNodes(iBdy,2) = sin(theta(iBdy)); % .* (1 + 0.1*sin(8*theta(iBdy)));

xyNodes(iBdy,1) = xyNodes(iBdy,1) - 0.5*sin(theta(iBdy));
xyNodes(iBdy,2) = xyNodes(iBdy,2) + 0.5*cos(theta(iBdy));

tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad);
poi = PoissonFEM2D(tnMesh);

figure(2); clf
tnMesh.plotMesh();
hold on
xs = linspace(-1, 1, 20);
ys = linspace(-1, 1, 20);
[xx,yy] = ndgrid(xs,ys);
plot(xx,yy,'k.');
axis image

%% Set up the problem

dirichletPredicate = @(x,y) 1; % everything is dirichlet
freeChargeFunc = @(x,y) x*y;
dirichletFunc = @(x,y) 0; %double(x>0);
neumannFunc = @(x,y) 0;

f = FEMProblem(poi);
[iDirichlet, iNeumann] = f.classifyBoundary(dirichletPredicate);
f.setDirichlet(iDirichlet, dirichletFunc);
f.setNeumann(iNeumann, neumannFunc);
f.setFreeCharge(freeChargeFunc);

%% Objective function

objFun = @(u) u(28);
DobjFun = @(u) double( (1:length(u)) == 28 );

%% Solve it

f.solve(objFun);
fprintf('F = %0.4e\n', f.F);

f.solveAdjoint(DobjFun);

%% Plot the field

xCoarse = linspace(-1.2, 1.2, 40);
yCoarse = linspace(-1.2, 1.2, 40);

figure(1); clf
u = f.poi.tnMesh.rasterizeField(f.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
hold on
f.poi.tnMesh.plotMesh('color', 'w');


%% Perturb Dirichlet

delta = 1e-8;

for ii = 1:length(f.iDirichlet)
    g = f.perturbedDirichlet(ii, delta);
    g.solve(objFun);
    
    dF_meas = (g.F - f.F)/delta;
    dF_calc = f.dF_dDirichlet(ii);
    
    fprintf('%0.4e vs %0.4e\n', dF_meas, dF_calc);
    
    xy = g.poi.tnMesh.getNodeCoordinates();
    tsi = scatteredInterpolant(xy(:,1), xy(:,2), g.u, 'linear', 'none');
    
    xs = linspace(-1.2, 1.2, 200);
    ys = linspace(-1.2, 1.2, 200);
    [xx,yy] = ndgrid(xs,ys);
    figure(1); clf
    u = tsi(xx,yy);
    imagesc_centered(xs, ys, u'); axis xy image
    colorbar
    hold on
    g.poi.tnMesh.plotMesh('color', 'w');
    %pause
end

%% Perturb free charge

xCoarse = linspace(-1.2, 1.2, 40);
yCoarse = linspace(-1.2, 1.2, 40);

xy = f.poi.tnMesh.getNodeCoordinates();

numNodes = length(f.u);
for ii = 1:numNodes
    g = f.perturbedFreeCharge(ii, delta);
    g.solve(objFun);
    
    dF_meas = (g.F - f.F)/delta;
    dF_calc = f.dF_dCharge(ii);
    
    fprintf('%0.4e vs %0.4e\n', dF_meas, dF_calc);
    
    %figure(1); clf
    %u = g.poi.tnMesh.rasterizeField(g.u, xCoarse, yCoarse);
    %imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    %hold on
    %g.poi.tnMesh.plotMesh('color', 'w');
    %plot(xy(ii,1), xy(ii,2), 'wo');
    %pause
end

%% Perturb Neumann

%% Perturb geometry nodes

delta = 1e-8;

xy = f.poi.tnMesh.xyNodes;

for mm = 1:f.poi.tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        g = f.perturbedMesh(mm, dirIdx, delta);
        g.solve(objFun);
        
        dF_meas = (g.F - f.F)/delta;
        dF_calc = f.dF_dxy(mm,dirIdx);
        
        dCharge_meas = (g.freeCharge - f.freeCharge)/delta;
        if dirIdx == 1
            dCharge_calc = f.dFreeCharge_dx(:,mm);
        else
            dCharge_calc = f.dFreeCharge_dy(:,mm);
        end
        
        disp([dCharge_meas, dCharge_calc]);
        
        fprintf('Measured %0.7e expected %0.7e\n', dF_meas, dF_calc);
        
        %figure(1); clf
        %u = g.poi.tnMesh.rasterizeField(g.u, xCoarse, yCoarse);
        %imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
        %hold on
        %g.poi.tnMesh.plotMesh('color', 'w');
        %plot(xy(mm,1), xy(mm,2), 'wo');
        %pause
    end
end

%% System matrices

xyNodes = poi.tnMesh.getNodeCoordinates();

A = poi.getSystemMatrix();
B = poi.getRhsMatrix();
NM = poi.getNeumannMatrix();

%% Partition into edge and center nodes

numNodes = poi.tnMesh.hFieldNodes.getNumNodes();

iBoundary = poi.tnMesh.hFieldNodes.getBoundaryNodes();

iDirichlet = iBoundary;
iNeumann = setdiff(iBoundary, iDirichlet);
iCenter = setdiff(1:numNodes, iDirichlet);

xyNodesCenter = xyNodes(iCenter,:);

%% Matrices to solve

A_center = A(iCenter, iCenter);
B_center = B(iCenter, iCenter);
A_dirichlet = A(iCenter, iDirichlet);
NM_center = NM(iCenter, iNeumann);

%% Inhomogeneous terms

freeCharge_center = zeros(length(iCenter),1);
freeCharge_center(2) = 1.0;
u0_dirichlet = zeros(length(iDirichlet),1);
en_neumann = zeros(length(iNeumann), 1);

%% Solve the forward problem

u_center = A_center \ (B_center*freeCharge_center + NM_center*en_neumann - A_dirichlet*u0_dirichlet);

u_all = zeros(numNodes,1);
u_all(iDirichlet) = u0_dirichlet;
u_all(iCenter) = u_center;

%% Evaluate the objective function

f_val = objFun(u_center);
Df_val = DobjFun(u_center);

%% Solve the adjoint system

v_center = A_center' \ Df_val';
v_all = 0*u_all;
v_all(iCenter) = v_center;

%% Matrix sensitivities

dA_dv = poi.getSystemMatrixSensitivity();
dB_dv = poi.getRhsMatrixSensitivity();
dNM_dv = poi.getNeumannMatrixSensitivity();

%% Sensitivity to free charge

dF_dFreeCharge_center = v_center' * B_center;
dF_dFreeCharge = 0*u_all; dF_dFreeCharge(iCenter) = dF_dFreeCharge_center;

%% Sensitivity to boundary electric field

dF_dEn_center = v_center' * NM_center;
dF_dEn = 0*u_all; dF_dEn(iNeumann) = dF_dEn_center;

%% Sensitivity to boundary voltage

dF_du0_center = -v_center' * A_dirichlet;
dF_du0 = 0*u_all; dF_du0(iDirichlet) = dF_du0_center;

%% Sensitivity to node perturbations

numNodes = tnMesh.hGeomNodes.getNumNodes();

dFdv_total = zeros(numNodes,2);

for vv = 1:numNodes
    
    wx = -dA_dv{1,vv}(iCenter, iCenter)*u_center ...
        - dA_dv{1,vv}(iCenter, iDirichlet)*u0_dirichlet ...
        + dNM_dv{1,vv}(iCenter, iNeumann)*en_neumann ...
        + dB_dv{1,vv}(iCenter, iCenter)*freeCharge_center; % ...
        %+ B_center*dFreeCharge_dv{vv,1}(iCenter);

    wy = -dA_dv{2,vv}(iCenter, iCenter)*u_center ...
        - dA_dv{2,vv}(iCenter, iDirichlet)*u0_dirichlet ...
        + dNM_dv{2,vv}(iCenter, iNeumann)*en_neumann ...
        + dB_dv{2,vv}(iCenter, iCenter)*freeCharge_center; %...
        %+ B_center*dFreeCharge_dv{vv,2}(iCenter);

    dFdvx = v_center'*wx;
    dFdvy = v_center'*wy;

    dFdv_total(vv,1) = dFdvx;
    dFdv_total(vv,2) = dFdvy;
end


%% See the answer!
tic
xs = linspace(-1.5, 1.5, 400);
ys = linspace(-1.5, 1.5, 400);
[xx,yy] = ndgrid(xs,ys);
%%
II = poi.tnMesh.getInterpolationOperator(xx(:),yy(:));
toc
%%

tic
II = poi.tnMesh.getRasterInterpolationOperator([-1.5, -1.5], [1.5, 1.5], [400, 400]);
toc
%%

[Dx, Dy] = poi.tnMesh.getGradientOperators();
ex = -Dx*u_all;
ey = -Dy*u_all;

%%

%z = reshape(II*u_all, size(xx));
z = reshape(II*u_all, size(xx));
z(z == 0) = nan;
figure(1); clf
imagesc_centered(xs, ys, z'); %, [-1, 1]);
colormap orangecrush
colorbar
axis xy image vis3d
hold on
%contour(xs, ys, z', 'w--');
%poi.tnMesh.plotMesh('Color', 'y');
plot(xyNodes(iDirichlet(1),1), xyNodes(iDirichlet(1),2), 'wo');
%plot(xyNodes(:,1), xyNodes(:,2), 'w.');
title('Poisson solution with Dirichlet boundary condition');
%quiver(xyNodes(:,1), xyNodes(:,2), ex, ey, 'w');
%sl = streamslice(xx', yy', reshape(II*ex, size(xx))', reshape(II*ey, size(xx))');
%set(sl, 'color', 'w');
%%

myFunc = @(x,y) sin(x) + cos(y);
f = poi.evaluateOnNodes(myFunc);

%%
f2 = reshape(II*f, size(xx));


%%

%hFineNodes = NodalTopology(tnMesh.hMesh, 10);
%[rrr,sss] = hFineNodes.basis.getNodes();
%xyFine = tnMesh.getFaceCoordinates(rrr, sss);

subFaces = tnMesh.hGeomNodes.getNodalMesh();
xyFine = xyNodes;
%subFaces = hFineNodes.getNodalMesh();
%lngFine = LinearNodalGeometry(subFaces, xyNodes(1:length(vertices),:), 10);
%xyFine = lngFine.getNodeCoordinates();

figure(1); clf
patch('Faces', subFaces, 'Vertices', xyFine, 'FaceAlpha', 0);
tnMesh.plotMesh()

%%

figure(1); clf
tnMesh.plotMesh('Color', 'b', 'LineStyle', '--');
%VVMesh.plotFV(faces, vertices, 'k--');
hold on
plot(xyNodes(:,1), xyNodes(:,2), 'bo');
title('Geometry nodes')

xy = tnMesh.getNodeCoordinates();
figure(2); clf
%VVMesh.plotFV(faces, vertices, 'k--');
tnMesh.plotMesh('Color', 'b', 'LineStyle', '--');
hold on
plot(xy(:,1), xy(:,2), 'bo');
title('Field nodes')

%%



%%

[dxy_dr, dxy_ds] = tnMesh.getQuadratureJacobian(1);
[dxy_dr, dxy_ds] = tnMesh.getFieldJacobian(1);

%%

tnMesh.hGeomNodes.basis1d.getInteriorNodes()


%%
% The geometry nodes need to be positioned at roughly the usual node
% positions and then displaced off of those by some amount.  The field
% and quadrature nodes' Jacobians depend on the geometry nodes and their
% coordinates depend on the geometry nodes and blahblah.
%
% My question is how do I actually position the geom nodes in the first
% place?  I have in TriNodalMesh the following useful functions:
% - getVertexNodeCoordinates
% - getEdgeInteriorNodeCoordinates
% - getEdgeNodeCoordinates() [COPYPASTA]
% - getFaceInteriorNodeCoordinates
% - getFaceNodeCoordinates
% - getNodeCoordinates
% which use (r,s) and the vertex coordinates to get their job done.
%
% In fact, I could begin with a vertex node topology and the vertex list
% and use that to fill out the other values.  Then the node coordinate
% functions would be
% - getVertexNodeCoordinates(geomNodes, outNodes)
% ...
% - getFieldNodeCoordinates()
% - getGeomNodeCoordinates()
% - getQuadNodeCoordinates()
%
% This is a little icky too, though.  I'm starting to have to make three
% versions of each major public function in TriNodalMesh because of my
% three nodal sets.  How can I refactor this stuff?
%
% NodalTopology is just an indexing tool.  These node coordinate
% calculators are only going to be used now for finding the geometry nodes
% and after that not used again.
%
% Perhaps I need a simple nodal mesh class and a curvilinear nodal mesh
% class.
%
% Arranging the nodes is sort of a meshing task.  I'd expect that for
% purposes of shape optimization I would somehow penalize the node
% positions so that they can't get too extreme.  Sometimes the nodes solve
% an elastostatics problem.  TriNodalMesh is no longer the appropriate
% place, maybe, for any of this node positioning stuff.
%
% Then TriNodalMesh must take the node positions in its constructor.
% 

% Now when I want to get node positions and sensitivities out...




