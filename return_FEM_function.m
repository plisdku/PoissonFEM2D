function [F, dFdp] = return_FEM_function(fem, p, Lx, Ly)
     
    [~] = fem.instantiateProblem(p);
   
    %disp('Instantiated Problem')
    [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);
    figure(107)
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    measBox = [-190e-3, 0, 55e-3, 0.7e-3]; %[-0.5, -0.5, 0.5, 0.5];
    measNxy = [250, 600];
    
    fprintf('Forward solution... ');
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);
    fprintf('Solved Forward ')
    [F, DF, xv, x_p, y_p, Nt] = ElectronSetupdemo2D_new(femProblem.uCartesian, xCoarse, yCoarse, measBox, measNxy); 

     fprintf('Adjoint solution... ');

     femProblem.solveAdjointCartesian(DF);
     fprintf('complete.\n');
     dFdp = femProblem.dF_dxy(:,1)' * dnx_dp + femProblem.dF_dxy(:,2)' * dny_dp ...
        + femProblem.dF_dDirichlet * dDirichlet_dp;

     figure(100); clf
     
%      xCoarse = linspace(-Lx, Lx, 200);
%     yCoarse = linspace(0, Ly, 200);
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    ax = axis;
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
    [ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
    for i = 1:size(xv,2)
    plot(xv(ix_x,i),xv(ix_y,i),'w', 'LineWidth', 1)
    plot(xv(ix_x(1),i),xv(ix_y(1),i),'rx')
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