%% Contour points (oriented correctly) and mesh size parameters

L = 1;

%%

isAxisymmetric = 1;

N_field = 2;
N_geom = 2;
N_quad = N_field+1;

s = 0.05; % mesh scale

geom2d = ParameterizedGeometry2D();
geom2d.addContour(@(p) [-L, L, L, -L], @(p) [0, 0, L, L], s,...
    2, [2,4], ...
    1, [1,3]);

fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(2, @(p,x,y) 0.0);
fem.setDirichlet(1, @(p,x,y) y);
fem.setFreeCharge(@(p,x,y) 0.0);

%%
[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblemNew(p0);

measBox = [-1, 1, 0, 2]; %[-0.5, -0.5, 0.5, 0.5];
measNxy = [20, 20];

objFun = @(u_Cartesian) sum(u_Cartesian(:));
DobjFun = @(u_Cartesian) ones(size(u_Cartesian));

fprintf('Forward solution... ');
femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);

%%

figure(1); clf
xCoarse = linspace(-L, L, 200);
yCoarse = linspace(0, L, 200);
u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
colormap orangecrush
hold on
femProblem.poi.tnMesh.plotMesh();

plot([instance.geometry.vertices(instance.geometry.lines(:,1),1), instance.geometry.vertices(instance.geometry.lines(:,2),1)]', ...
    [instance.geometry.vertices(instance.geometry.lines(:,1),2), instance.geometry.vertices(instance.geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)

figure(2); clf
plot(yCoarse, u(1:5:end,:)', '-');
xlabel('r');
ylabel('V');