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


N_field = 3;
N_geom = 2;
N_quad = N_field;

f = FEMInterface(N_field, N_geom, N_quad);


f.addContour(@(p) [-3, 3, 3, -3], @(p) [-3, -3, 3, 3], @(p) 2, 'neumann', @(p,x,y) 0.0);
f.addContour(@(p) [-2, -1.5, -1.5, -2] + p(1), @(p) [-2, -2, 2, 2], @(p) 1, 'dirichlet', @(p, x, y) p(2));
f.addContour(@(p) [1.5, 2, 2, 1.5], @(p) [-2, -2, 2, 2], @(p) 1, 'dirichlet', @(p, x, y) -1.0);


% f.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], @(p) 0.5*10*[0.5, 4, 0.5, 4, 4, 4]*1e-3, 'neumann', @(p,x,y) 0.0);
% 
% f.addContour(@(p) [-Lx+5e-3, -Lx+6e-3, -Lx+6e-3, -Lx+5e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d], @(p) 0.5*10*0.5e-3, 'dirichlet', @(p,x,y) -20);
% 
% f.addContour(@(p) [-Lx+10e-3, -Lx+12e-3, -Lx+12e-3, -Lx+10e-3], @(p) [rAperture2, rAperture2, Ly-d, Ly-d], @(p) 5*0.5e-3, 'dirichlet', @(p,x,y) 1.0);
% f.addContour(@(p) [0, L_lens, L_lens, 0] + [p(1), 0, 0, 0], @(p) [rAperture2, rAperture2, Ly-d2, Ly-d2], @(p) 5*0.5e-3, 'dirichlet', @(p,x,y) 5);
% f.addContour(@(p) [5e-3, 5e-3+L_lens, 5e-3+L_lens, 5e-3], @(p) [rAperture2, rAperture2, Ly-d2, Ly-d2], @(p) 5*0.5e-3, 'dirichlet', @(p,x,y) 0.0);
% 
% f.addContour(@(p) [13e-3, 14e-3, 14e-3, 13e-3], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], @(p) 5*0.5e-3, 'dirichlet', @(p,x,y) 30);

f.setFreeCharge(@(p,x,y) 0.0);

p = [0,0,0,0,0,0,0,0];
[femProblem, geometry, dnx_dp, dny_dp] = f.instantiateProblem(p);

xy = femProblem.poi.tnMesh.xyNodes;

%%
figure()
plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
    [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'b', 'linewidth', 3)
hold on
femProblem.poi.tnMesh.plotMesh()




%% Sensitivity

p0 = [0];
iParamToVary = 1;

deltas = linspace(0, 0.5e-3, 10) + 0.5e-3;

Fs = 0*deltas;
DFs = 0*deltas;

for nn = 1:length(deltas)
    p = p0;
    p(iParamToVary) = p(iParamToVary) + deltas(nn);
    
    fprintf('Instantiating...\n');
    [femProblem, geometry, dDirichlet_dp, dnx_dp, dny_dp] = f.instantiateProblem(p);
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    % Build and evaluate objective function and gradient
    %outI = femProblem.poi.tnMesh.getRasterInterpolationOperator([-8.5e-3, 0], [14.5e-3, 1.6e-3], [10, 10]);
    %outI = femProblem.poi.tnMesh.getRasterInterpolationOperator([-6e-3 -3.5e-3],[3e-3 5e-3],[10, 10]);
    %objFun = @(u_nodal) sum(outI * u_nodal);
    %DobjFun = @(u_nodal) sum(outI, 1);
    
    measBox = [-6e-3, 3e-3, -3.5e-3, 5e-3];
    measNxy = [5, 5];
    
    objFun = @(u_Cartesian) sum(u_Cartesian(:));
    DobjFun = @(u_Cartesian) ones(size(u_Cartesian));
    
    
    fprintf('Forward solution... ');
    %femProblem.solve(objFun);
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    F = objFun(femProblem.uCartesian);
    fprintf('Adjoint solution... ');
    femProblem.solveAdjointCartesian(DobjFun(femProblem.uCartesian));
    %femProblem.solveAdjoint(DobjFun);
    fprintf('complete.\n');
    
    
    % Sensitivity to parameters
    dFdp = femProblem.dF_dxy(:,1)' * dnx_dp +femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;
    
    %Fs(nn) = femProblem.F;
    Fs(nn) = F;
    DFs(nn) = dFdp(iParamToVary);
    
    figure(1); clf
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    %id = femProblem.iDirichlet;
    id = ':';
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    axis xy image
    title(sprintf('Iteration %i', nn));
    
    if nn > 1
        figure(2); clf
        df_meas = gradient(Fs, deltas);
        plot(1:nn, DFs(1:nn));
        yl = ylim;
        hold on
        plot(1:nn, df_meas(1:nn));
        legend('Adjoint', 'Calculated');
        ylim(yl)
    end
    pause(0.01)
end

%%
