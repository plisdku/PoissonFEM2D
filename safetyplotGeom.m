%% Plot Initial Geometry with a p vector out of zeros

function safetyplotGeom(geom2d, p)

    if isrow(p) == 1 
        p = p';
    end
     geometry = geom2d.evaluateGeometry(p);
    figure(77)
       plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
    [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o-', 'color', [0.8 0.8 0.8], 'linewidth', 2)
end 