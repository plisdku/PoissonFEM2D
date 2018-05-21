%% Contour points (oriented correctly) and mesh size parameters
%close all
Lx = 60e-3;
Lx_outer = 20e-3;
Ly = 45e-3;
rAperture1 = 3e-3;
rAperture2 = 1.5e-3;
rAperture3 = 1e-3;
d = 34e-3;
Lx_outer_2 = 120e-3;


isAxisymmetric = 1;


%%

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this
geom2d = ParameterizedGeometry2D();
%p0 = [0,0,0,0,0,0,0,0]';
s = 2.1e-3; % mesh scale
ratio = 0.3;

l2 = 5e-3;
l1 = 2e-3;
l3 = 2e-3;
Vb = 29e3;
D1 = 3e-3;
L1 = 5e-3;
e_d = 7.5e-3;
Wd = 43e-3;
Wd2 = 164e-3;


        geom2d.addContour(@(p) [-Lx-Lx_outer-Lx_outer_2, -Lx+e_d+Wd, 0,-Lx+e_d+Wd+l1+L1+l2+L1+l3, Lx+Lx_outer, Lx+Lx_outer, 0, -Lx-Lx_outer-Lx_outer_2], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[5, ratio, ratio, ratio, 5, 8, 8, 8], 1, 1:8);
        geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd]+ p([1 2 2 1])', @(p) [D1, D1, Ly-d, Ly-d] + p([3 3 4 4])', s*ratio*[1, 1, 2, 2], 2, 1:4);
        geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1]+ p([5 6 6 5])', @(p) [D1, D1, Ly-d, Ly-d] + p([7 7 8 8])', s*ratio*[1, 1, 2, 2], 3, 1:4);
        geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1] + p([9 10 10 9])', @(p) [D1, D1, Ly-d, Ly-d] + p([11 11 12 12])', s*ratio*[1, 1, 2, 2], 4, 1:4);



fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb );
fem.setDirichlet(4, @(p,x,y) 0 );
%fem.setDirichlet(5, @(p,x,y) 0.0 + p(4));

fem.setFreeCharge(@(p,x,y) 0.0);
%%


%% Sensitivity

p0 = zeros(12,1);
iParamToVary = 3;
[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.instantiateProblem(p0);
figure(10764); clf

        %imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
        ax = axis;
        colormap orangecrush
        hold on
        femProblem.poi.tnMesh.plotMesh();

        geometry = fem.instantiatedGeom.geometry;
        plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
            [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)

deltas = linspace(0.0, 0.3, 5)*1e-3;
%deltas = linspace(0.0, 2.5, 10);
%deltas = deltas(2:end);

%deltas = 10e-3;

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
    
    figure(101010); clf
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    title('temp');
    %pause
    
    measBox = [-Lx-Lx_outer_2-Lx_outer*0.5, 0, 65e-3, 0.7e-3]; %[-0.5, -0.5, 0.5, 0.5];
    measNxy = [250, 300];
    
    objFun = @(u_Cartesian) sum(u_Cartesian(:));
    DobjFun = @(u_Cartesian) ones(size(u_Cartesian));
    
    fprintf('Forward solution... ');
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    xCoarse = linspace(measBox(1), measBox(3), measNxy(1));
    yCoarse = linspace(measBox(2), measBox(4), measNxy(2));
    %F = objFun(femProblem.uCartesian);
    [particles, hit_objective] = SetupParticles_gradient();
    [VV] = ElectronSetup_obj(femProblem.uCartesian, measBox, measNxy, particles, hit_objective); 
figure(1); clf

    F = VV.Fval;
    DF = VV.dFdV;
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    
    
    
    %id = femProblem.iDirichlet;
   % id = ':';
   % quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
   % quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
    axis xy image
    title(sprintf('Iteration %i', nn));
    %Nt = 1500;
    %[ix_x, ix_y, ix_z, ~] =...
    %        get_Index3D(VV.ParticleArray(1).Nt);
    for ii = 1:length(VV.ParticleArray)
            plot(VV.ParticleArray(ii).xx,VV.ParticleArray(ii).yy,'g', 'LineWidth', 1)
            plot(VV.ParticleArray(ii).xx,VV.ParticleArray(ii).yy,'rx')
    end
    
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
    
    figure(19); clf
    %xCoarse = linspace(-3, 3, 200);
    %yCoarse = linspace(-3, 3, 200);
    xCoarse = linspace(-Lx-Lx_outer-Lx_outer_2, Lx+Lx_outer, 200);
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
   % Nt = 400;
    [ix_x, ix_y, ix_z, ~] =...
            get_Index3D(Nt);
    for ii = 1:size(xv_all,2)
            plot(xv_all(ix_x,ii),xv_all(ix_y,ii),'g', 'LineWidth', 1)
            %plot(xv(ix_x_all(1),ii),xv_all(ix_y(1),ii),'rx')
    end
    
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



