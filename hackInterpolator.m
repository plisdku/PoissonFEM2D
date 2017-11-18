%% Point in element grid?

[vertices, faces] = VVMesh.wagonWheel(5);
vertices = vertices(:,1:2);

% Node orders
N_field = 4;
N_geom = 5;
N_quad = N_field;

lng = LinearNodalGeometry(faces, vertices, N_geom);
xyNodes = lng.getNodeCoordinates();
xyNodes(:,1) = xyNodes(:,1) - 0.15*sin(xyNodes(:,2)*2*pi);
xyNodes(:,2) = xyNodes(:,2) + 0.1*sin(xyNodes(:,2)*2*pi) + 0.1*cos(xyNodes(:,1)*2*pi);

tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad);

figure(1); clf
tnMesh.plotMesh();

%% Test points

xs = linspace(min(xyNodes(:,1))-0.1, max(xyNodes(:,1))+0.1, 200);
ys = linspace(min(xyNodes(:,2))-0.1, max(xyNodes(:,2))+0.1, 200);
[xx,yy] = ndgrid(xs,ys);
%%

enclosingFace = zeros(size(xx));

figure(1); clf
tnMesh.plotMesh();
hold on

numFaces = tnMesh.hMesh.getNumFaces();

for ff = 4 %1:numFaces
    xy = tnMesh.getFaceBoundary(ff);
    extents = max(xy) - min(xy);
    
    xy0 = min(xy) - 0.1*extents;
    xy1 = max(xy) + 0.1*extents;
    
    ii = find(xs > xy0(1) & xs < xy1(1));
    jj = find(ys > xy0(2) & ys < xy1(2));
    [iii,jjj] = ndgrid(ii,jj);
    
    [rs, bad] = tnMesh.inverseCoordinateTransform(ff, xx(ii,jj), yy(ii,jj));
    
    good = ~bad;
    fprintf('%i good\n', nnz(good));
    
    enclosingFace(sub2ind(size(xx), iii(good), jjj(good))) = ff;
    
    xxx = xx(ii,jj); yyy = yy(ii,jj);
    
    if exist('q', 'var')
        %delete(p);
        delete(q);
        delete(r);
    end
    
    q = plot(xy(:,1), xy(:,2), 'ro');
    %p = plot(xx(ii,jj), yy(ii,jj), 'k.');
    r = plot(xxx(good), yyy(good), 'g.');
    pause(0.1)
end

