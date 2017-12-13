%% Test a whole adjoint thing


%% Contour points (oriented correctly) and mesh size parameters

Lx = 6;
Ly = 5;
rAperture1 = 0.5;
rAperture2 = 1.5;
rAperture3 = 1;

%%

f = FEMInterface();
f.addContour(@(p) [-Lx, Lx, Lx, -Lx], @(p) [0, 0, Ly, Ly], @(p) 1, 'neumann', @(p,x,y) 0.0);
f.addContour(@(p) [-Lx+1, -Lx+2, -Lx+2, -Lx+1], @(p) [rAperture1, rAperture1, Ly-1, Ly-1], @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);
f.addContour(@(p) [-Lx+3, -Lx+4, -Lx+4, -Lx+3], @(p) [rAperture2, rAperture2, Ly-1, Ly-1], @(p) 0.5, 'dirichlet', @(p,x,y) 1.0);
f.addContour(@(p) [Lx-2, Lx-1, Lx-1, Ly-1], @(p) [rAperture3, rAperture3, Ly-1, Ly-1], @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);

f.setFreeCharge(@(p,x,y) 0.0);

p = [];
[ff, vv] = f.solve(p, @(u) 1);

%%

xCoarse = linspace(-Lx, Lx, 20);
yCoarse = linspace(0, Ly, 20);

tic
u = f.femProblem.poi.tnMesh.rasterizeField(f.femProblem.u, xCoarse, yCoarse); % slow as hell
toc
%%

figure(1); clf
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
colormap orangecrush
hold on
f.femProblem.poi.tnMesh.plotMesh('color', 'w');

%% Your potential:

electromagneticPotential = f.femProblem.u;



