%% Contour points (oriented correctly) and mesh size parameters
close all
Lx = 6e-3;
Ly = 8e-3;
rAperture1 = 0.5e-3;
rAperture2 = 1.5e-3;
rAperture3 = 1e-3;
d = 4e-3;

isAxisymmetric = 1;

%%

N_field = 3;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this

p0 = [0,0,0,0,0,0,0,0]';
s = 1e-3; % mesh scale

geom2d = ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], s*[0.5, 4, 0.5, 4, 4, 4], 1, 1:6);
geom2d.addContour(@(p) [-Lx+1e-3, -Lx+2e-3, -Lx+2e-3, -Lx+1e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d]+p(5:8)', s*0.5, 2, 1:4);
geom2d.addContour(@(p) [-Lx+3e-3, -Lx+4e-3, -Lx+4e-3, -Lx+3e-3]+p(1:4)', @(p) [rAperture2, rAperture2, Ly-d, Ly-d], s*0.5, 3, 1:4);
geom2d.addContour(@(p) [Lx-2e-3, Lx-1e-3, Lx-1e-3, Lx-2e-3], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], s*0.5, 2, 1:4);
geom2d.addContour(@(p) [Lx-5e-3, Lx-4e-3, Lx-4e-3, Lx-5e-3], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], s*0.5, 2, 1:4);

fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0.0);
fem.setDirichlet(3, @(p,x,y) 1.0);
fem.setFreeCharge(@(p,x,y) 0.0);

%%

[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblem(p0);

%% Sensitivity

p0 = [0,0,0,0,0,0,0,0]';
iParamToVary = 3;

deltas = linspace(0.0, 0.1, 80)*1e-3;
%deltas = linspace(0.0, 2.5, 10);
%deltas = deltas(2:end);

Fs = 0*deltas;
DFs = 0*deltas;
ps = 0*deltas;

for nn = 1:length(deltas)
    p = p0;
    p(iParamToVary) = p(iParamToVary) + deltas(nn);
    ps(nn) = p(iParamToVary);
    
    fprintf('Instantiating...\n');
    
    DO_ADJUST_MESH = 1;  % Lars pay attention to this!!!!
    
    if DO_ADJUST_MESH
        [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);
    else
        [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblem(p);
    end
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    measBox = [-5.5, 0, 5.5, 0.45]*1e-3; %[-0.5, -0.5, 0.5, 0.5];
    measNxy = [150, 50];
    
    objFun = @(u_Cartesian) sum(u_Cartesian(:));
    DobjFun = @(u_Cartesian) ones(size(u_Cartesian));
    
    fprintf('Forward solution... ');
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    %F = objFun(femProblem.uCartesian);
    [F, DF, DFDEx, DFDEr, E_x] = ElectronSetupdemo2D(femProblem.uCartesian, xCoarse, yCoarse, measBox, measNxy); 

    fprintf('Adjoint solution... ');
    %femProblem.solveAdjoint(DobjFun, DoutI);
    %femProblem.solveAdjointCartesian(DobjFun(femProblem.uCartesian));
    femProblem.solveAdjointCartesian(DF);
    fprintf('complete.\n');
    
    % Sensitivity to parameters
    dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;
    
    Fs(nn) = F; %femProblem.F;
    DFs(nn) = dFdp(iParamToVary);
    
    figure(1); clf
    %xCoarse = linspace(-3, 3, 200);
    %yCoarse = linspace(-3, 3, 200);
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    %id = femProblem.iDirichlet;
    id = ':';
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
    axis xy image
    title(sprintf('Iteration %i', nn));
    
    figure(2); clf
    subplot(2,1,1);
    plot(ps(1:nn), DFs(1:nn), 'o-');
    %yl = ylim;
    hold on
    if nn >= 2
        df_meas = gradient(Fs(1:nn), deltas(1:nn));
        plot(ps(1:nn), df_meas, 'o-');
        %ylim(yl);
    end
    xlabel('p')
    ylabel('gradient')
    legend('Adjoint', 'Calculated');

    subplot(2,1,2);
    plot(ps(1:nn), Fs(1:nn), 'o-');
    xlabel('p')
    ylabel('F')

    pause(0.01)
end



