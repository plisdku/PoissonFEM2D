%% Test a whole adjoint thing

assertClose = @(a,b) assert(norm(a(:)-b(:))/(norm(a(:)) + norm(b(:))) < 1e-5 || norm(a(:) + b(:)) < 1e-50);

%% Set up the mesh

lx = [0, 1, 1, 0];
ly = [0, 0, 1, 1];

in_lx = 0.5 + 0.1*[-1, -1, 1, 1];
in_ly = 0.5 + 0.1*[-1, 1, 1, -1];

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

dirichletPredicate = @(x,y) norm(x-0.5) < 0.25 && norm(y-0.5) < 0.25;

femp = FEMProblem(poi);
[iDirichlet, iNeumann] = femp.classifyBoundary(dirichletPredicate);

freeChargeFunc = @(x,y) 0; %x*y;
dirichletFunc = @(x,y) 0; %double(x>0);
neumannFunc = @(x,y) 0;

femp.setDirichlet(iDirichlet, dirichletFunc);
femp.setNeumann(iNeumann, neumannFunc);
femp.setFreeCharge(freeChargeFunc);

%% Objective function

numFieldNodes = femp.poi.tnMesh.hFieldNodes.getNumNodes();
df = zeros(1, numFieldNodes);
df(femp.iCenter) = 1;

objFun = @(u) sum(u(femp.iCenter));  % TODO: make it work for sum(u), incl. Dirichlet nodes.  WILL NEED!!!
DobjFun = @(u) df;

objFun = @(u) sum(u);
DobjFun = @(u) ones(size(u))';

%% Solve it

femp.solve(objFun);
fprintf('F = %0.4e\n', femp.F);
femp.solveAdjoint(DobjFun);

%% Plot the field

xCoarse = linspace(-0.1, 1.1, 40);
yCoarse = linspace(-0.1, 1.1, 40);

figure(1); clf
u = femp.poi.tnMesh.rasterizeField(femp.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
hold on
femp.poi.tnMesh.plotMesh('color', 'w');



%% Perturb Dirichlet

delta = 1e-8;

for ii = 1:length(femp.iDirichlet)
    g = femp.perturbedDirichlet(ii, delta);
    g.solve(objFun);
    
    dF_meas = (g.F - femp.F)/delta;
    dF_calc = femp.dF_dDirichlet(ii);
    
    fprintf('%0.4e vs %0.4e\n', dF_meas, dF_calc);
    
    xy = g.poi.tnMesh.getNodeCoordinates();
    tsi = scatteredInterpolant(xy(:,1), xy(:,2), g.u, 'linear', 'none');
    
    xs = linspace(-0.1, 1.1, 200);
    ys = linspace(-0.1, 1.1, 200);
    [xx,yy] = ndgrid(xs,ys);
    figure(1); clf
    u = tsi(xx,yy);
    imagesc_centered(xs, ys, u'); axis xy image
    colorbar
    hold on
    g.poi.tnMesh.plotMesh('color', 'w');
    pause(0.01)
end


%% Perturb Neumann

delta = 1e-8;

for ii = 1:length(femp.iNeumann)
    g = femp.perturbedNeumann(ii, delta);
    g.solve(objFun);
    
    dF_meas = (g.F - femp.F)/delta;
    dF_calc = femp.dF_dNeumann(ii);
    
    fprintf('%0.4e vs %0.4e\n', dF_meas, dF_calc);
    
    xy = g.poi.tnMesh.getNodeCoordinates();
    tsi = scatteredInterpolant(xy(:,1), xy(:,2), g.u, 'linear', 'none');
    
    xs = linspace(-0.1, 1.1, 200);
    ys = linspace(-0.1, 1.1, 200);
    [xx,yy] = ndgrid(xs,ys);
    figure(1); clf
    u = tsi(xx,yy);
    imagesc_centered(xs, ys, u'); axis xy image
    colorbar
    hold on
    g.poi.tnMesh.plotMesh('color', 'w');
    pause
end


%% Perturb free charge

xCoarse = linspace(-0.1, 1.1, 40);
yCoarse = linspace(-0.1, 1.1, 40);

xy = femp.poi.tnMesh.getNodeCoordinates();

numNodes = length(f.u);
for ii = 1:numNodes
    g = femp.perturbedFreeCharge(ii, delta);
    g.solve(objFun);
    
    dF_meas = (g.F - femp.F)/delta;
    dF_calc = femp.dF_dCharge(ii);
    
    fprintf('%0.4e vs %0.4e\n', dF_meas, dF_calc);
    
    %figure(1); clf
    %u = g.poi.tnMesh.rasterizeField(g.u, xCoarse, yCoarse);
    %imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    %hold on
    %g.poi.tnMesh.plotMesh('color', 'w');
    %plot(xy(ii,1), xy(ii,2), 'wo');
    %pause
end


%% Perturb geometry nodes

delta = 1e-8;

xy = femp.poi.tnMesh.xyNodes;

for mm = 1:femp.poi.tnMesh.hGeomNodes.getNumNodes()
    for dirIdx = 1:2
        g = femp.perturbedMesh(mm, dirIdx, delta);
        g.solve(objFun);
        
        dF_meas = (g.F - femp.F)/delta;
        dF_calc = femp.dF_dxy(mm,dirIdx);
        
        dCharge_meas = (g.freeCharge - femp.freeCharge)/delta;
        if dirIdx == 1
            dCharge_calc = femp.dFreeCharge_dx(:,mm);
        else
            dCharge_calc = femp.dFreeCharge_dy(:,mm);
        end
        
        %disp([dCharge_meas, dCharge_calc]);
        
        fprintf('Measured %0.7e expected %0.7e\n', dF_meas, dF_calc);
        
%         figure(1); clf
%         u = g.poi.tnMesh.rasterizeField(g.u, xCoarse, yCoarse);
%         imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
%         hold on
%         g.poi.tnMesh.plotMesh('color', 'w');
%         plot(xy(mm,1), xy(mm,2), 'wo');
%         pause(0.01)
    end
end


