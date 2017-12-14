%% Contour points (oriented correctly) and mesh size parameters

Lx = 6;
Ly = 8;
rAperture1 = 0.5;
rAperture2 = 1.5;
rAperture3 = 1;
d = 4;

%%

N_field = 2;
N_geom = 2;
N_quad = N_field;

f = FEMInterface(N_field, N_geom, N_quad);
f.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], @(p) [0.5, 4, 0.5, 4, 4, 4], 'neumann', @(p,x,y) 0.0);
f.addContour(@(p) [-Lx+1, -Lx+2, -Lx+2, -Lx+1]+p(1:4), @(p) [rAperture1, rAperture1, Ly-d, Ly-d]+p(5:8), @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);
f.addContour(@(p) [-Lx+3, -Lx+4, -Lx+4, -Lx+3], @(p) [rAperture2, rAperture2, Ly-d, Ly-d], @(p) 0.5, 'dirichlet', @(p,x,y) 1.0);
f.addContour(@(p) [Lx-2, Lx-1, Lx-1, Lx-2], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);

f.setFreeCharge(@(p,x,y) 0.0);

p = [0,0,0,0,  0,0,0,0];
[femProblem, dnx_dp, dny_dp] = f.instantiateProblem(p);

%%

for nn = 1:size(dnx_dp, 2)
    figure(101); clf
    femProblem.poi.tnMesh.plotMesh();
    hold on
    
    xy = femProblem.poi.tnMesh.xyNodes;
    quiver(xy(:,1), xy(:,2), dnx_dp(:,nn), dny_dp(:,nn), 'linewidth', 2, 'color', 'r');
    pause
end


%%

tic
femProblem.solve(@(u) u(5));
toc

tic
femProblem.solveAdjoint(@(u) (1:length(u)) == 5);
toc

%%
figure(101); clf
femProblem.poi.tnMesh.plotMesh();
hold on
quiver(xy(:,1), xy(:,2), femProblem.dF_dxy(:,1), femProblem.dF_dxy(:,2), 'r-', 'linewidth', 2)
axis xy image

%%

xCoarse = linspace(-Lx, Lx, 100);
yCoarse = linspace(0, Ly, 100);

tic
u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse); % slow as hell
toc

%%

figure(1); clf
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
colormap orangecrush
hold on
quiver(xy(:,1), xy(:,2), femProblem.dF_dxy(:,1), femProblem.dF_dxy(:,2), 'w-', 'linewidth', 2)
quiver(xy(:,1), xy(:,2), femProblem.dF_dxy(:,1), femProblem.dF_dxy(:,2), 'g-', 'linewidth', 1)

femProblem.poi.tnMesh.plotMesh('color', 'w');



