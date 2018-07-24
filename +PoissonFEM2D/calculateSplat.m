function idxSplatTime = calculateSplat(xs, ys, vertices, contours)
%%

% contours = geometry.contourVertexIndices;
% vertices = geometry.vertices;
% 
% numContours = numel(contours);
% 
% [xs,ys] = fakeTrajectories(-6, 0, 6, 0, 100, 10);
% ys = abs(ys);

%%

%xs = xTraj;
%ys = yTraj;
%vertices = geometry.vertices;
%contours = geometry.contourVertexIndices;

%%

numT = size(xs,1);
numParticles = size(xs, 2);
numContours = length(contours);

hitMask = zeros(numT, numParticles);

idxSplatTime = zeros(numParticles,1);

for nn = 1:numContours
    xyContour = vertices(contours{nn},:);
    
    if nn == 1
        hit = ~inpolygon(xs, ys, xyContour(:,1), xyContour(:,2));
    else
        hit = inpolygon(xs, ys, xyContour(:,1), xyContour(:,2));
    end
    
    hitMask = hitMask | hit;
end

for pp = 1:numParticles
    idxHit = find(hitMask(:,pp), 1, 'first');

    if idxHit
        idxSplatTime(pp) = idxHit;
    end
end

%%

%idxSplatTime = calculateSplat(xTraj, yTraj, geometry.vertices, geometry.contourVertexIndices);
%[idxSplatTraj,~,idxSplatTime] = find(idxSplatTime);