%% Test a whole adjoint thing

assertClose = @(a,b) assert(norm(a(:)-b(:))/(norm(a(:)) + norm(b(:))) < 1e-5 || norm(a(:) + b(:)) < 1e-50);

%% Set up the mesh

lx = [0, 1, 1, 0];
ly = [0, 0, 1, 1];

in_lx = 0.5 + 0.05*[-1, -1, 1, 1];
in_ly = 0.5 + 0.05*[-1, 1, 1, -1];

density = 2;
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


femp = FEMProblem(poi);
[iDirichlet, iNeumann] = femp.classifyBoundary(dirichletPredicate);

%% Define some shit

x0 = 0.55;
y0 = 0.85;
sigma = 0.05;
freeChargeFunc = @(x,y) exp( (-(x-x0).^2 - (y-y0).^2)/(2*sigma^2));

dirichletFunc = @(x,y) 0; %double(x>0);
neumannFunc = @(x,y) 0;

femp.setDirichlet(iDirichlet, dirichletFunc);
femp.setNeumann(iNeumann, neumannFunc);
femp.setFreeCharge(freeChargeFunc);

%% Objective function

iArbitraryFace = 13;
arbitraryInteriorNodes = tnMesh.hFieldNodes.getFaceInteriorNodes(iArbitraryFace);
iArbitraryNode = arbitraryInteriorNodes(1);
%%
for iArbitraryNode = 1:500
%iArbitraryNode = 1;

%objFun = @(u) sum(u);
%DobjFun = @(u) ones(size(u))';
objFun = @(u) u(iArbitraryNode);
DobjFun = @(u) double( (1:length(u)) == iArbitraryNode);

%% Solve it

femp.solve(objFun);
fprintf('F = %0.4e\n', femp.F);

femp.solveAdjoint(DobjFun);

%% Plot the field

xCoarse = linspace(-0.2, 1.2, 40);
yCoarse = linspace(-0.2, 1.2, 40);

figure(1); clf
u = femp.poi.tnMesh.rasterizeField(femp.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
hold on
femp.poi.tnMesh.plotMesh('color', 'w');


%% Perturb (Dirichlet)

delta = 1e-6;
femp2 = FEMProblem(poi);
femp2.setDirichlet(iDirichlet, dirichletFunc);
femp2.setNeumann(iNeumann, neumannFunc);
femp2.setFreeCharge(freeChargeFunc);
femp2.u0_dirichlet(1) = femp2.u0_dirichlet(1) + delta;

femp2.solve(objFun);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_dDirichlet(1);

fprintf('Measured %g, calculated %g\n', dF_meas, full(dF_calc));
assertClose(dF_meas, dF_calc);
end
%% Perturb (Neumann)

femp2 = FEMProblem(N, domainF, domainV, dirichletPredicate);
femp2.setSources(freeCharge, dirichletVal, neumannVal);
femp2.en_neumann(1) = femp2.en_neumann(1) + delta;
femp2.solve(measPt);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_den(1);

fprintf('Measured %g, calculated %g\n', dF_meas, full(dF_calc));
assertClose(dF_meas, dF_calc);

%% Perturb (free charge)

femp2 = FEMProblem(N, domainF, domainV, dirichletPredicate);
femp2.setSources(freeCharge, dirichletVal, neumannVal);
femp2.freeCharge(10) = femp2.freeCharge(1) + delta;
femp2.solve(measPt);

dF_meas = (femp2.F - femp.F)/delta;
dF_calc = femp.dF_df(10);

fprintf('Measured %g, calculated %g\n', dF_meas, full(dF_calc));
assertClose(dF_meas, dF_calc);

%% Perturb (vertices)

delta = 1e-6;

iXY = 1;
for iVert = 1:femp.fem.meshNodes.hMesh.getNumVertices()
    disp(iVert)

    femp2 = FEMProblem(N, domainF, domainV, dirichletPredicate);
    femp2.fem.meshNodes.vertices(iVert,iXY) = femp.fem.meshNodes.vertices(iVert,iXY) + delta;
    femp2.setSources(freeCharge, dirichletVal, neumannVal);
    femp2.solve(measPt);

    dF_meas = (femp2.F - femp.F)/delta;
    dF_calc = femp.dFdv_total(iVert,iXY);
    assertClose(dF_meas, dF_calc);
    
    nf = @(A) norm(full(A));
    fprintf('%g vs %g\n', nf(dF_meas), nf(dF_calc));
%    fprintf('Rel err %g/%f = %g\n', norm(full(dF_meas-dF_calc)), ...
%        norm(full(dF_calc)), ...
%        norm(full(dF_meas-dF_calc))/norm(full(dF_calc)));

    if 0
        figure(11); clf
        VVMesh.plotFV(domainF, domainV, 'k-');
        vertex = domainV(iVert,:);
        hold on
        plot(x0, y0, 'rx');
        plot(vertex(1), vertex(2), 'ro', 'MarkerSize', 10);
        title(sprintf('%g vs %g\n', nf(dF_meas), nf(dF_calc)));
        pause
    end
    
    dA_meas = (femp2.A - femp.A)/delta;
    dA_calc = femp.dA_dv{iVert,iXY};
    assertClose(dA_meas, dA_calc);
    
    dB_meas = (femp2.B - femp.B)/delta;
    dB_calc = femp.dB_dv{iVert,iXY};
    assertClose(dB_meas, dB_calc);
    
    dNM_meas = (femp2.NM - femp.NM)/delta;
    dNM_calc = femp.dNM_dv{iVert,iXY};
    assertClose(dNM_meas, dNM_calc);
    
    dFreeCharge_meas = (femp2.freeCharge - femp.freeCharge)/delta;
    dFreeCharge_calc = femp.dFreeCharge_dv{iVert, iXY};
    assertClose(dFreeCharge_meas, dFreeCharge_calc);
    
end


% TODO: check all the pieces

%% Picture?

xs = linspace(-0.1, 1.1, 400);
ys = linspace(-0.1, 1.1, 400);
[xx, yy] = ndgrid(xs, ys);

II = femp.fem.meshNodes.getInterpolationOperator(xx(:), yy(:));

%%

xy_nodes = femp.fem.meshNodes.getNodeCoordinates();
%%

[Dx, Dy] = femp.fem.meshNodes.getGradientOperators();
Ex = Dx*femp.u;
Ey = Dy*femp.u;

magE = sqrt(abs(Ex).^2 + abs(Ey).^2);

%%

dF_dud = 0*femp.u;
dF_dud(femp.iNodesDirichlet) = femp.dF_dud;

dF_den = 0*femp.u;
dF_den(femp.iNodesNeumann) = femp.dF_den;
%%
u_grid = II*femp.u;
%u_grid = II*femp.dF_df;
u_grid = reshape(u_grid, size(xx));

%%

figure(1); clf
imagesc_centered(xs, ys, u_grid'); %, [0, 0.1]);

%%
colormap orangecrush(0.7)
%colorbar
%
hold on
%VVMesh.plotFV(domainF, domainV, 'w-', 'linewidth', 0.001)
%plot(measPt(1), measPt(2), 'go', 'markersize', 10);
%plot(x0, y0, 'rx', 'markersize', 10);
hold on
%plot(meshNodes.vertices(:,1), meshNodes.vertices(:,2), 'wo');
%plot(xy_nodes(:,1), xy_nodes(:,2), 'w.', 'MarkerSize', 2)
colorbar
axis xy image vis3d
title('Basis interpolation')

%%
quiver(xy_nodes(:,1), xy_nodes(:,2), Ex, Ey, 'w', 'linewidth', 0.5);

%%

%plot(domainV(:,1), domainV(:,2), 'wo')
hold on
quiver(domainV(:,1), domainV(:,2), femp.dFdv_total(:,1), femp.dFdv_total(:,2), 'r', 'linewidth', 2);



