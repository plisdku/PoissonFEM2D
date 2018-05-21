%% Contour points (oriented correctly) and mesh size parameters

Lx = 15e-3;
Ly = 20e-3;
rAperture1 = 3e-3;
rAperture2 = 8e-3;
rAperture3 = 2e-3;
d = 5e-3;
d2 = 15e-3;
L_lens = 3e-3;
%%

N_field = 2;
N_geom = 2;
N_quad = N_field;

f = FEMInterface(N_field, N_geom, N_quad);

f.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], @(p) [0.5, 4, 0.5, 4, 4, 4]*1e-3, 'neumann', @(p,x,y) 0.0);

f.addContour(@(p) [-Lx+5e-3, -Lx+6e-3, -Lx+6e-3, -Lx+5e-3]+p(1:4), @(p) [rAperture1, rAperture1, Ly-d, Ly-d]+p(5:8), @(p) 0.5e-3, 'dirichlet', @(p,x,y) -20);

f.addContour(@(p) [-Lx+10e-3, -Lx+12e-3, -Lx+12e-3, -Lx+10e-3], @(p) [rAperture2, rAperture2, Ly-d, Ly-d], @(p) 0.5e-3, 'dirichlet', @(p,x,y) 1.0);
f.addContour(@(p) [0, L_lens, L_lens, 0], @(p) [rAperture2, rAperture2, Ly-d2, Ly-d2], @(p) 0.5e-3, 'dirichlet', @(p,x,y) 5);
f.addContour(@(p) [5e-3, 5e-3+L_lens, 5e-3+L_lens, 5e-3], @(p) [rAperture2, rAperture2, Ly-d2, Ly-d2], @(p) 0.5e-3, 'dirichlet', @(p,x,y) 0.0);

f.addContour(@(p) [13e-3, 14e-3, 14e-3, 13e-3], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], @(p) 0.5e-3, 'dirichlet', @(p,x,y) 30);

f.setFreeCharge(@(p,x,y) 0.0);

p = [0,0,0,0,  0,0,0,0];
[femProblem, dnx_dp, dny_dp] = f.instantiateProblem(p);

xy = femProblem.poi.tnMesh.xyNodes;

%%

for nn = 1:size(dnx_dp, 2)
    figure(101); clf
    femProblem.poi.tnMesh.plotMesh();
    hold on
    
    quiver(xy(:,1), xy(:,2), dnx_dp(:,nn), dny_dp(:,nn), 'linewidth', 2, 'color', 'r');
    pause
end


%%

%[DoutI, outI] = femProblem.poi.tnMesh.getRasterInterpolationOperatorSensitivity([-1, 1], [1, 3], [10, 10]);
outI = femProblem.poi.tnMesh.getRasterInterpolationOperator([-8.5e-3, 0], [14.5e-3, 1.6e-3], [10, 10]);

%%
%lars = @(u_raster) sum(sum(u_raster));

objFun = @(u_nodal) sum(outI * u_nodal);
DobjFun = @(u_nodal) sum(outI, 1);
%%

tic
%femProblem.solve(@(u) u(5));
femProblem.solve(objFun);
toc
%%
tic
%femProblem.solveAdjoint(@(u) (1:length(u)) == 5);
femProblem.solveAdjoint(@(u) dGdV_xr_row*outI);
toc

%%
figure(1); clf


xCoarse = linspace(-Lx, Lx, 100);
yCoarse = linspace(0, Ly, 100);
%u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
colormap orangecrush
hold on
femProblem.poi.tnMesh.plotMesh();

id = femProblem.iDirichlet;


quiver(xy(id,1), xy(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
quiver(xy(id,1), xy(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
axis xy image
plot(xv(ix_x,1),xv(ix_y,1),'w')

