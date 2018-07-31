import PoissonFEM2D.*
%Lx = 60e-3;
Lx_outer = 60e-3;
%Ly = 45e-3;
%rAperture1 = 3e-3;
%rAperture2 = 1.5e-3;
%rAperture3 = 1e-3;
%d = 34e-3;

isAxisymmetric = 1;

%%

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this

%p0 = [0,0,0,0,0,0,0,0]';
s = 1.3e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.3;
geom2d = PoissonFEM2D.ParameterizedGeometry2D();

%l2 = 5e-3;
%l1 = 2e-3;
%l3 = 2.5e-3;
Vb = 29e3;
D1 = 3e-3;
%L1 = 5e-3;
%L2 = 2e-3;
%e_d = 7.5e-3;
%Wd = 43e-3;

%Lx  = 80e-3;
%D_outer = Ly - d;

Ly = 45e-3;

D12 = 5e-3;
D23 = D12;
L1 = 2e-3;
L2 = 5e-3;
L3 = L1;

r_i = 3e-3;
r_o = 11e-3;

Ncontour_vec = [30 30 30 30]/3;

E2_V1 = [-L2/2, r_i];
E2_V2 = [L2/2, r_i];
E2_V3 = [L2/2, r_o];
E2_V4 = [-L2/2, r_o];
E2 = [E2_V1; E2_V2; E2_V3; E2_V4];
[E2_x, E2_y, s2_vec] = BoxVector(E2, Ncontour_vec, @(N) s_vec_lin(N,1,2));
l2 = length(E2_x);
l22 = 2*l2;

E1_V1 = [-L2/2, r_i] - [D12+L1 0];
E1_V2 = [-L2/2, r_i] - [D12 0];
E1_V3 = [-L2/2, r_o] - [D12 0];
E1_V4 = [-L2/2, r_o] - [D12+L1 0];
E1 = [E1_V1; E1_V2; E1_V3; E1_V4];
[E1_x, E1_y, s1_vec] = BoxVector(E1, Ncontour_vec, @(N) s_vec_lin(N,1,2));
l1 = length(E1_x);
l12 = 2*l1;
%
E3_V1 = [L2/2, r_i] + [D23 0];
E3_V2 = [L2/2, r_i] + [D23+L3 0];
E3_V3 = [L2/2, r_o] + [D23+L3 0];
E3_V4 = [L2/2, r_o] + [D23 0];
E3 = [E3_V1; E3_V2; E3_V3; E3_V4];
[E3_x, E3_y, s3_vec] = BoxVector(E3, Ncontour_vec, @(N) s_vec_lin(N,1,2));
l3 = length(E3_x);
l32 = 2*l3;
%
plotElectrodePoints(E1_x, E1_y, E2_x, E2_y, E3_x, E3_y)

[E1x_func, E1y_func, end_p1] = createGeomFunction(E1_x, E1_y, 0);
[E2x_func, E2y_func, end_p2] = createGeomFunction(E2_x, E2_y, end_p1);
[E3x_func, E3y_func, end_p3] = createGeomFunction(E3_x, E3_y, end_p2);


geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);

geom2d.addContour(E1x_func, E1y_func, s*ratio*s1_vec, 2, 1:l1);
geom2d.addContour(E2x_func, E2y_func, s*ratio*s2_vec, 3, 1:l2);
geom2d.addContour(E3x_func, E3y_func, s*ratio*s3_vec, 4, 1:l3);

safetyplotGeom(geom2d, zeros(1,end_p3))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb );
fem.setDirichlet(4, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p3)';


fn_handle = @(fem) @(p) return_FEM_functionobjruntime(fem, p, Lx_outer, Ly); 

minX = ones(1,end_p3)'*(-5e-3);
maxX = ones(1,end_p3)'*(10e-3);
t1 = tic;

[x, fval, iter, xHist, fHist, DfHist] = extremize12(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', @ plotfunc, 'MaxIter', 25);
toc(t1)
save('ShowcaseData' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist')

