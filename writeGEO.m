function writeGEO(fname, contours, meshSizes)

fh = fopen(fname, 'w');

if ~iscell(meshSizes)
    originalMeshSizes = meshSizes;
    
    meshSizes = cell(length(contours), 1);
    for cc = 1:length(contours)
        meshSizes{cc} = repmat(originalMeshSizes(cc), size(contours{cc},1), 1);
    end
end

%% Write the points and lines

pointIdx = 1;
lineIdx = 1;
lineLoopIdx = 100000;

lineLoops = [];

contourLines = {};

for cc = 1:length(contours)
    points = contours{cc};
    numPoints = size(points, 1);

    
    iPointsInContour = pointIdx + (1:numPoints) - 1;
    iLinesInContour = lineIdx + (1:numPoints) - 1;
    
    contourLines{cc} = iLinesInContour;
    
    for pp = 1:numPoints
        fprintf(fh, 'Point(%i) = {%0.9g, %0.9g, 0.0, %0.9g};\n', pointIdx, points(pp,1), points(pp,2), meshSizes{cc}(pp));
        pointIdx = pointIdx + 1;
        if pointIdx > 1856
            
            fprintf('Point(%i) = {%0.9g, %0.9g, 0.0, %0.9g};\n', pointIdx, points(pp,1), points(pp,2), meshSizes{cc}(pp));

        end
        if pointIdx == 1857
            disp('at error')
            %break
        end
    end
    
    for ll = 1:numPoints
        p0 = iPointsInContour(ll);
        p1 = iPointsInContour(1 + mod(ll, numPoints));
        fprintf(fh, 'Line(%i) = {%i, %i};\n', lineIdx, p0, p1);
        lineIdx = lineIdx + 1;
    end
    
    fprintf(fh, 'Line Loop(%i) = {', lineLoopIdx);
    fprintf(fh, '%i, ', iLinesInContour(1:end-1));
    fprintf(fh, '%i};\n', iLinesInContour(end));
    lineLoops(end+1) = lineLoopIdx;
    lineLoopIdx = lineLoopIdx + 1;
    
end

%% Surface to mesh

fprintf(fh, 'Plane Surface(1) = {');
for nn = 1:length(lineLoops)-1
    fprintf(fh, '%i, ', lineLoops(nn));
end
fprintf(fh, '%i};\n', lineLoops(end));

%% Physical entities (labels)

fprintf(fh, 'Physical Surface("Everwhere") = {1};\n');

for cc = 1:length(contourLines)
    fprintf(fh, 'Physical Line(%i) = {', cc);
    fprintf(fh, '%i, ', contourLines{cc}(1:end-1));
    fprintf(fh, '%i};\n', contourLines{cc}(end));
end

try
fclose(fh);
catch
    warning('Closing failed')
end
end
