%% Element matrix tests

faces = [1 2 3];
vertices = [0,0; 1,0; 0,1];

% Node orders
N_field = 4;
N_geom = 2;
N_quad = 4;

lng = LinearNodalGeometry(faces, vertices, N_geom);
xyNodes = lng.getNodeCoordinates();

tnMesh = TriNodalMesh(faces, xyNodes, N_field, N_geom, N_quad);
poi = PoissonFEM2D(tnMesh);

xy = poi.tnMesh.getNodeCoordinates();

%% Test element potential matrix
% -(Dx'*Q*Dx + Dy'*Q*Dy)

% f = x, g = x
% f'*A*g = -integral(1*1) = -0.5
f = xy(:,1);
g = xy(:,1);
expected = -0.5;

%A = poi.getElementPotentialMatrix(1);
A = poi.getSystemMatrix();

result = f'*A*g;

fprintf('Got %2.6f, expected %2.6f\n', result, expected);

%%
% f = y, g = y
% also should get -0.5
f = xy(:,2);
g = xy(:,2);
expected = -0.5;

result = f'*A*g;
fprintf('Got %2.6f, expected %2.6f\n', result, expected);

%%
% f = x^2 + y
% g = x
% f'*A*g = -integral(2x) = -2*(1/6) = -1/3
f = xy(:,1).^2 + xy(:,2);
g = xy(:,1);
expected = -1/3;

result = f'*A*g;
fprintf('Got %2.6f, expected %2.6f\n', result, expected);

%% Center nodes?

numNodes = poi.tnMesh.hFieldNodes.getNumNodes();

iBoundary = poi.tnMesh.hFieldNodes.getBoundaryNodes();
iCenter = setdiff(1:numNodes, iBoundary);

figure(1); clf
plot(xy(iBoundary,1), xy(iBoundary,2), 'bx');
hold on
plot(xy(iCenter,1), xy(iCenter,2), 'ro');
legend('Boundary', 'Center');




