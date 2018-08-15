function FullResultsPlot(nFigure, femProblem, VV, geometry, idx_arrows,...
    measBox, xyGeomNodes, Lx, Ly)

    figure(nFigure); clf
    
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u,...
        xCoarse, yCoarse);
    
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    ax = axis;
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();

    plot([geometry.vertices(geometry.lines(:,1),1), ...
        geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2),...
        geometry.vertices(geometry.lines(:,2),2)]', ...
        'color', [0.8 0.8 0.8], 'linewidth', 2)


    for i = 1:length(VV.ParticleArray)
        plot(VV.ParticleArray(i).xx,VV.ParticleArray(i).yy,'w',...
            'LineWidth', 1)
    end

    id = idx_arrows;
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), ...
        femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2),...
        'w-', 'linewidth', 2)
    quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), ...
        femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2),...
        'g-', 'linewidth', 1)
    
    plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
    axis(ax)




end 