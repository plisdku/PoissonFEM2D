function [v, f] = meshPolygon(lx, ly, density)

lx = reshape(lx, 1, []);
ly = reshape(ly, 1, []);

numEdges = length(lx);

%% subdivide the edges as needed

% wrap
lx2 = [lx, lx(1)];
ly2 = [ly, ly(1)];
lxx = [];
lyy = [];

for ee = 1:numEdges
    edgeLength = sqrt( (lx2(ee+1)-lx2(ee))^2 + (ly2(ee+1)-ly2(ee))^2 );
    x = linspace(lx2(ee), lx2(ee+1), density*edgeLength + 1);
    y = linspace(ly2(ee), ly2(ee+1), density*edgeLength + 1);
    
    lxx = [lxx, x(1:end-1)];
    lyy = [lyy, y(1:end-1)];
end

%% the innards

xBounds = [min(lxx), max(lxx)];
yBounds = [min(lyy), max(lyy)];

Lx = diff(xBounds);
Ly = diff(yBounds);

Nx = density*Lx;
Ny = density*Ly;

xs = xBounds(1) + Lx*rand(1, Nx*Ny);
ys = yBounds(1) + Ly*rand(1, Nx*Ny);

i_inside = inpolygon(xs, ys, lx, ly);

xs = xs(i_inside);
ys = ys(i_inside);

%% Delaunay triangulation

numEdgePts = length(lxx);
constraints = [ 1:numEdgePts; [2:numEdgePts,1] ]';

xxx = [lxx, xs]';
yyy = [lyy, ys]';

numIters = 100;

for iter = 1:numIters
    tri = delaunayTriangulation( [ xxx, yyy ], constraints );

    xxx = tri.Points(:,1);
    yyy = tri.Points(:,2);
    numPoints = length(xxx);

    numTriEdges = size(tri.edges, 1);
    edges = tri.edges;
    adjacency = sparse(edges(:,1), edges(:,2), ones(numTriEdges,1), numPoints, numPoints);
    adjacency = adjacency + adjacency';
    
    % averaging step

    distx = diff(xxx(edges), [], 2);
    disty = diff(yyy(edges), [], 2);
    dist_total2 = distx.^2 + disty.^2;
    dist_total2(dist_total2 == 0) = 1;

    A = sparse(edges(:,1), edges(:,2), dist_total2, numPoints, numPoints);
    A = A + A';

    meanX = A*tri.Points(:,1) ./ sum(A, 2);
    meanY = A*tri.Points(:,2) ./ sum(A, 2);

    newX = xxx + 0.5*(meanX - xxx);
    newY = yyy + 0.5*(meanY - yyy);
    
    
    oldx = xxx;
    oldy = yyy;
    
    xxx(numEdgePts+1:end) = newX(numEdgePts+1:end);
    yyy(numEdgePts+1:end) = newY(numEdgePts+1:end);
    
    deltaX = oldx - xxx;
    deltaY = oldy - yyy;
    delta2 = deltaX.^2 + deltaY.^2;
    
    convergence = 1e-2;
    %fprintf('Test %g against %g\n', max(delta2), convergence/(density^2));
    if max(delta2) < convergence/(density^2)
        %fprintf('Iter %d\n', iter);
        break
    end    
end

v = tri.Points;
f = tri.ConnectivityList;
f = f(tri.isInterior(), :);