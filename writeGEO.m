function writeGEO(fname, contours, meshSizes)

fh = fopen(fname, 'w');

%% Write the points and lines

pointIdx = 1;
lineIdx = 1;

lineLoops = [];

contourLines = {};

for cc = 1:length(contours)
    points = contours{cc};
    numPoints = size(points, 1);
    
    meshSize = meshSizes(cc);
    
    iPointsInContour = pointIdx + (1:numPoints) - 1;
    iLinesInContour = lineIdx + (1:numPoints) - 1;
    
    contourLines{cc} = iLinesInContour;
    
    for pp = 1:numPoints
        fprintf(fh, 'Point(%i) = {%i, %i, 0.0, %i};\n', pointIdx, points(pp,1), points(pp,2), meshSize);
        pointIdx = pointIdx + 1;
    end
    
    for ll = 1:numPoints
        p0 = iPointsInContour(ll);
        p1 = iPointsInContour(1 + mod(ll, numPoints));
        fprintf(fh, 'Line(%i) = {%i, %i};\n', lineIdx, p0, p1);
        lineIdx = lineIdx + 1;
    end
    
    fprintf(fh, 'Line Loop(%i) = {', lineIdx);
    fprintf(fh, '%i, ', iLinesInContour(1:end-1));
    fprintf(fh, '%i};\n', iLinesInContour(end));
    lineLoops(end+1) = lineIdx;
    lineIdx = lineIdx + 1;
    
end

%% Surface to mesh

fprintf(fh, 'Plane Surface(1) = {');
fprintf(fh, '%i, ', lineLoops(1:end-1));
fprintf(fh, '%i};\n', lineLoops(end));

%% Physical entities (labels)

fprintf(fh, 'Physical Surface("Everwhere") = {1};\n');

for cc = 1:length(contourLines)
    fprintf(fh, 'Physical Line(%i) = {', cc);
    fprintf(fh, '%i, ', contourLines{cc}(1:end-1));
    fprintf(fh, '%i};\n', contourLines{cc}(end));
end
