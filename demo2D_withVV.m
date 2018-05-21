%% Contour points (oriented correctly) and mesh size parameters

Lx = 6;
Ly = 8;
rAperture1 = 0.5;
rAperture2 = 1.5;
rAperture3 = 1;
d = 4;

%%

N_field = 3;
N_geom = 2;
N_quad = N_field;

p0 = [0,0,0,0,0,0,0,0]';
s = 2; % mesh scale

geom2d = ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], @(p) s*[0.5, 4, 0.5, 4, 4, 4], 'neumann', @(p,x,y) 0.0);
geom2d.addContour(@(p) [-Lx+1, -Lx+2, -Lx+2, -Lx+1], @(p) [rAperture1, rAperture1, Ly-d, Ly-d]+p(5:8)', @(p) s*0.5, 'dirichlet', @(p,x,y) 0.0);
geom2d.addContour(@(p) [-Lx+3, -Lx+4, -Lx+4, -Lx+3]+p(1:4)', @(p) [rAperture2, rAperture2, Ly-d, Ly-d], @(p) s*0.5, 'dirichlet', @(p,x,y) 1000.0*1e3);
geom2d.addContour(@(p) [Lx-2, Lx-1, Lx-1, Lx-2], @(p) [rAperture3, rAperture3, Ly-d, Ly-d], @(p) s*0.5, 'dirichlet', @(p,x,y) 0.0);

instance = InstantiatedGeometry2D(geom2d, N_field, N_geom, N_quad);
instance.instantiateMesh(p0);

%%

fem = FEMInterface(N_field, N_geom, N_quad);
fem.setFreeCharge(@(p,x,y) 0.0);
[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblemNew(p0, instance, geom2d);

%%


%% Sensitivity

p0 = [0,0,0,0,0,0,0,0]';
iParamToVary = 3;

%deltas = linspace(0.0, 0.2, 15);
deltas = linspace(0.1, 0.15, 10);
%deltas = linspace(0.0, 2.5, 10);
deltas = deltas(2:end);

Fs = 0*deltas;
DFs = 0*deltas;
ps = 0*deltas;

for nn = 1:length(deltas)
    
    if nn > 1
        u_old = femProblem.uCartesian; 
        F_old = F; 
        E_x_old = E_x; 
        %E_r_old = E_r;
    end
    p = p0;
    p(iParamToVary) = p(iParamToVary) + deltas(nn);
    ps(nn) = p(iParamToVary);
    
    fprintf('Instantiating...\n');
    
    DO_ADJUST_MESH = 1;  % Lars pay attention to this!!!!
    
    if DO_ADJUST_MESH
        instance.adjustMesh(p);
    else
        instance.instantiateMesh(p);
    end
    [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblemNew(p, instance, geom2d);
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    measBox = [-5.5, 0, 5.5, 0.45]; %[x_min y_min (1) x_max y_max][-0.5, -0.5, 0.5, 0.5];
    measNxy = [100, 30];
    
    objFun = @(u_Cartesian) sum(u_Cartesian(:));
    DobjFun = @(u_Cartesian) ones(size(u_Cartesian));
    
    
    fprintf('Forward solution... ');
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    %UU2 = femProblem.uCartesian;
    %save('UU2_save','UU2')
    %F = objFun(femProblem.uCartesian);
    [F, DF, DFDEx, DFDEr, E_x] = ElectronSetupdemo2D(femProblem.uCartesian, xCoarse, yCoarse, measBox, measNxy); 

    fprintf('Adjoint solution... ');
    %femProblem.solveAdjoint(DobjFun, DoutI);
    %femProblem.solveAdjointCartesian(DobjFun(femProblem.uCartesian));
    femProblem.solveAdjointCartesian(DF);%DobjFun(femProblem.uCartesian));
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
    
    plot([instance.vertices(instance.lines(:,1),1), instance.vertices(instance.lines(:,2),1)]', ...
        [instance.vertices(instance.lines(:,1),2), instance.vertices(instance.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
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
    
    %x_grid = linspace(measBox(1),measBox(3),measNxy(1));
    %y_grid = linspace(measBox(2),measBox(4),measNxy(2));
    %r_grid = [-y_grid(end:-1:2) y_grid];

    %d_x = x_grid(2) - x_grid(1);
    %d_r = y_grid(2) - y_grid(1);
    %E_x = -centeredDiff(femProblem.uCartesian,1) / d_x; 
    %E_r = -centeredDiff(femProblem.uCartesian,2) / d_r; 

    
    if nn > 1
        %delta_Ex = E_x - E_x_old;
        %delta_F_adj_Ex = sum(sum(DFDEx .* delta_Ex))
        
        %delta_Er = E_r - E_r_old; 
        %delta_F_adj_Er = sum(sum(DFDEr .* delta_Er))

        delta_u = femProblem.uCartesian - u_old; 
        delta_F_adj = sum(sum(DF .* delta_u))
        
        delta_F = F - F_old
        delta_delta_F = delta_F_adj./ delta_F  
        %delta_delta_F_Ex = delta_F_adj_Ex ./delta_F
        %delta_delta_F_Er = delta_F_adj_Er ./delta_F
    end
    
    pause(0.01)
end



