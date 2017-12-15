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
%f.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], @(p) [0.5, 4, 0.5, 4, 4, 4], 'neumann', @(p,x,y) 0.0);
%f.addContour(@(p) [-Lx+1, -Lx+2, -Lx+2, -Lx+1]+p(1:4), @(p) [rAperture1, rAperture1, Ly-d, Ly-d]+p(5:8), @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);
%f.addContour(@(p) [-Lx+3, -Lx+4, -Lx+4, -Lx+3], @(p) [rAperture2, rAperture2, Ly-d, Ly-d], @(p) 0.5, 'dirichlet', @(p,x,y) 1.0);
%f.addContour(@(p) [Lx-2, Lx-1, Lx-1, Lx-2], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], @(p) 0.5, 'dirichlet', @(p,x,y) 0.0);

f.addContour(@(p) [-3, 3, 3, -3], @(p) [-3, -3, 3, 3], @(p) 2, 'neumann', @(p,x,y) 0.0);
f.addContour(@(p) [-2, -1.5, -1.5, -2] + p(1), @(p) [-2, -2, 2, 2], @(p) 1, 'dirichlet', @(p, x, y) p(2));
f.addContour(@(p) [1.5, 2, 2, 1.5], @(p) [-2, -2, 2, 2], @(p) 1, 'dirichlet', @(p, x, y) -1.0);

f.setFreeCharge(@(p,x,y) 0.0);

p = [0,0];
[femProblem, geometry, dDirichlet_dp, dnx_dp, dny_dp] = f.instantiateProblem(p);

xy = femProblem.poi.tnMesh.xyNodes;

%%

plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
    [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'b', 'linewidth', 3)
hold on
femProblem.poi.tnMesh.plotMesh()

%% Sensitivity

p0 = [0, 1];
iParamToVary = 1;

deltas = linspace(0, 0.5, 10);

Fs = 0*deltas;
DFs = 0*deltas;

for nn = 1:length(deltas)
    p = p0;
    p(iParamToVary) = p(iParamToVary) + deltas(nn);
    
    fprintf('Instantiating...\n');
    [femProblem, geometry, dDirichlet_dp, dnx_dp, dny_dp] = f.instantiateProblem(p);
    xy = femProblem.poi.tnMesh.xyNodes;
    
    % Build and evaluate objective function and gradient
    outI = femProblem.poi.tnMesh.getRasterInterpolationOperator([-0.5, -0.5], [0.5, 0.5], [5, 5]);
    objFun = @(u_nodal) sum(outI * u_nodal);
    DobjFun = @(u_nodal) sum(outI, 1);
    
    fprintf('Forward solution... ');
    femProblem.solve(objFun);
    fprintf('Adjoint solution... ');
    femProblem.solveAdjoint(DobjFun);
    fprintf('complete.\n');
    
    % Sensitivity to parameters
    dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;
    
    Fs(nn) = femProblem.F;
    DFs(nn) = dFdp(iParamToVary);
    
    figure(1); clf
    xCoarse = linspace(-3, 3, 200);
    yCoarse = linspace(-3, 3, 200);
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    id = femProblem.iDirichlet;
    quiver(xy(id,1), xy(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
    quiver(xy(id,1), xy(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    axis xy image
    title(sprintf('Iteration %i', nn));
    
    if nn > 1
        figure(2); clf
        df_meas = gradient(Fs, deltas);
        plot(1:nn, DFs(1:nn));
        hold on
        plot(1:nn, df_meas(1:nn));
        legend('Adjoint', 'Calculated');
    end
    pause(0.01)
end

%%

