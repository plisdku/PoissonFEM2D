%% Plot a curvilinear mesh

[vertices, faces] = VVMesh.wagonWheel(3);
vertices = vertices(:, 1:2);

N_displ = 4;
N_field = N_displ;

fieldMesh = TriNodalMesh(N_field, faces, vertices);
displacementMesh = TriNodalMesh(N_displ, faces, vertices);
xyNodes_original = displacementMesh.getNodeCoordinates();
numEdges = displacementMesh.getNumEdges();
numFaces = displacementMesh.getNumFaces();


%% Set the displacements

ux = zeros(displacementMesh.getNumNodes(), 1);
uy = zeros(size(ux));

%ux = 0.05*randn(size(ux));
%uy = 0.05*randn(size(uy));

% Let's try to make a circle.

iBoundaryNodes = displacementMesh.getBoundaryNodes();

for nn = iBoundaryNodes'
    
    xyNode = xyNodes_original(nn,:);
    
    theta = atan2(xyNode(2), xyNode(1));
    
    ux(nn) = cos(theta) - xyNode(1);
    uy(nn) = sin(theta) - xyNode(2);
    
end


%% Plot it

%%

rEdge = linspace(-1, 1, 50); % r-coordinate along a single edge
V_low = displacementMesh.basis1d.V;
V_high = support.vandermonde(N_displ, rEdge);
V_nodes = support.vandermonde(N_displ, fieldMesh.basis1d.r);

interpToLine= V_high*inv(V_low);
interpToFieldNodes = V_nodes*inv(V_low);

figure(1); clf
hold on
for ee = 1:numEdges
    iNodes = displacementMesh.getEdgeNodes(ee);
    
    xyNodes = xyNodes_original(iNodes,:) + [ux(iNodes), uy(iNodes)];
    
    % Move along the edge curvy-like
    xys = interpToLine * xyNodes;
    
    %xyFieldNodes = interpToFieldNodes * xyNodes;
    
    plot(xys(:,1), xys(:,2), 'b-');
    %plot(xyFieldNodes(:,1), xyFieldNodes(:,2), 'bo');
end

%%


V_low = displacementMesh.basis.V;
V_nodes = support2d.vandermonde(N_displ, fieldMesh.basis.r, fieldMesh.basis.s);
interpToFieldNodes = V_nodes*inv(V_low);

for ff = 1:numFaces
    
    iNodes = displacementMesh.getFaceNodes(ff);
    
    xyNodes = xyNodes_original(iNodes,:) + [ux(iNodes), uy(iNodes)];
    
    xys = interpToFieldNodes * xyNodes;
    plot(xys(:,1), xys(:,2), 'bo');
    
end
axis xy image
grid on
