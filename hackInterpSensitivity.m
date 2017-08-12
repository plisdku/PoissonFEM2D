%% Hack together my interpolation matrix sensitivity thinggggg

N = 3;

triVerts = [0,0; 1,0; 0,1];
DxyTri = 0*triVerts;
DxyTri(1,1) = 1;
delta = 1e-6;
triVerts2 = triVerts + delta*DxyTri;
face = [1,2,3];

basis = BasisNodes(N);
rsNodes = basis.getNodes();

%% Interpolation
xs = linspace(-1, 2);
ys = linspace(-1, 2);
[xx,yy] = ndgrid(xs,ys);
xy_query = [xx(:)'; yy(:)'];

%%
rNodes = rsNodes(:,1);
interped = basis.interpolate(rNodes, triVerts', xy_query(1,:), xy_query(2,:));
interped = reshape(interped, size(xx));

interped2 = basis.interpolate(rNodes, triVerts2', xy_query(1,:), xy_query(2,:));
interped2 = reshape(interped2, size(xx));

Dinterped_meas = (interped2-interped)/delta;

%% Plot the interpolated thing

figure(1); clf
imagesc(xs, ys, interped');
axis xy image
hold on
VVMesh.plotFV(face, triVerts, 'w-', 'linewidth', 2);

%% Now let's try to get a derivative

%[V, dVdr, dVds] = support2d.vandermonde(N, rsNodes(:,1), rsNodes(:,2));

[T, x0] = support2d.rs2xy_affineParameters(triVerts');

rs_query = support2d.xy2rs(triVerts', xy_query);
rs_query2 = support2d.xy2rs(triVerts2', xy_query);
Drs_query_meas = (rs_query2-rs_query)/delta;
Drs = support2d.xy2rsSensitivity(xyTri', DxyTri', xy_query);

[V, dVdr, dVds] = support2d.vandermonde(N, rs_query(1,:), rs_query(2,:));
V2 = support2d.vandermonde(N, rs_query2(1,:), rs_query2(2,:));
DV_meas = (V2-V)/delta;
DV = bsxfun(@times, dVdr, Drs(1,:)') + bsxfun(@times, dVds, Drs(2,:)');

interpMatrix = V*basis.invV;
interpMatrix2 = V2*basis.invV;
DinterpMatrix_calc = (interpMatrix2-interpMatrix)/delta;
DinterpMatrix = DV*basis.invV;

%% Let's try it another way

DM = basis.interpolationMatrixSensivitity_nodal_xy(triVerts', DxyTri', xy_query(1,:), xy_query(2,:));
Dinterped_basis = DM*rNodes;
