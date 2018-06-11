function [F, dFdp] = return_FEM_functionobj(fem, p, Lx, Ly)
     
    [~] = fem.instantiateProblem(p);
   
    [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);
    figure(107)
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    measBox = [-55e-3, 0, 55e-3, 1e-3];
    measNxy = [350, 18000];
    
    fprintf('Forward solution... ');
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    fprintf('Solved Forward ')
    
    u_cartesian = femProblem.poi.tnMesh.rasterizeField(femProblem.u,...
        linspace(measBox(1),measBox(3),measNxy(1)), ...
        linspace(measBox(2),measBox(4),measNxy(2)));
    %[F, DF, xv, x_p, y_p, Nt] = ElectronSetupdemo2D_new(femProblem.uCartesian, xCoarse, yCoarse, measBox, measNxy); 
        [particles, hit_objective] = SetupParticles_spottest();
        [VV] = ElectronSetup_obj(u_cartesian, measBox, measNxy, particles, hit_objective); 
    
     fprintf('Adjoint solution... ');

     femProblem.solveAdjointCartesian(VV.dFdV);
     fprintf('complete.\n');
     dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;
     F = VV.Fval;
     
     figure(100); clf
     

    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    ax = axis;
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
   
    for i = 1:length(VV.ParticleArray)
    plot(VV.ParticleArray(i).xx,VV.ParticleArray(i).yy,'w', 'LineWidth', 1)
    %plot(VV.ParticleArray(i).xv(ix_x(1)),VV.ParticleArray(i).xv(ix_y(1)),'rx')
    end
    %plot(x_p,y_p,'bx')
%     xlim([-Lx Lx])
%     ylim([0 Ly])


    %id = femProblem.iDirichlet;
    id = ':';
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
    %axis xy image
    axis(ax)

end 