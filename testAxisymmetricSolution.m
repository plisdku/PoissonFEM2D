%% Contour points (oriented correctly) and mesh size parameters

L = 1;
r0 = 0.1;
r1 = 1.0;

%%

isAxisymmetric = 1;

N_field = 3;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 0.5/4; % mesh scale

geom2d = ParameterizedGeometry2D();
geom2d.addContour(@(p) [-L, L, L, -L], @(p) [r0, r0, r1, r1], s, 1, [1,3], 2, [2,4]);

fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(2, @(p,x,y) 0.0);
fem.setDirichlet(1, @(p,x,y) y);
fem.setFreeCharge(@(p,x,y) 0.0);

%%

[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblem(p0);

measBox = [-1, 1, 0, 2]; %[-0.5, -0.5, 0.5, 0.5];
measNxy = [20, 20];

objFun = @(u_Cartesian) sum(u_Cartesian(:));
DobjFun = @(u_Cartesian) ones(size(u_Cartesian));

fprintf('Forward solution... ');
femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);

%%

figure(1); clf
xCoarse = linspace(-L, L, 200);
yCoarse = linspace(r0, r1, 200);
u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
colormap orangecrush
hold on
femProblem.poi.tnMesh.plotMesh();

geometry = fem.instantiatedGeom.geometry;
plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
    [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)

%%

u_exact = (r1-r0)/(log(r1) - log(r0)) * log(yCoarse/r0) + 0.1;
u_fem = u([1,end/2,end],:);

u_error = bsxfun(@minus, u_fem, u_exact);

errorVal = norm(u_error(:));

figure(2); clf
plot(yCoarse, u_fem');
hold on
plot(yCoarse, u_exact, 'b--', 'linewidth', 2);
grid on
legend('u(1)', 'u(end/2)', 'u(end)', 'exact', 'location', 'best')
xlabel('r')
ylabel('V')
title(sprintf('N = %i, s = %0.2e, err = %0.2e', N_field, s, errorVal));