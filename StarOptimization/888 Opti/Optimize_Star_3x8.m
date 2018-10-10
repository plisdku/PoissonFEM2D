import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.28e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.3;
geom2d = PoissonFEM2D.ParameterizedGeometry2D();


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;
L1 = 2e-3;
L2 = 5e-3;
L3 = L1;

r_i = 3e-3;
r_o = 11e-3;

Ncontour_vec = [3 3 3 3];%[100 100 100 100];%[30 30 30 30];

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

[E1x_func, E1y_func, end_p11] = createGeomFunction(E1_x, E1_y, 0);
[E2x_func, E2y_func, end_p12] = createGeomFunction(E2_x, E2_y, end_p11);
[E3x_func, E3y_func, end_p13] = createGeomFunction(E3_x, E3_y, end_p12);


r_alpha1 = sqrt( (0.5*L1)^2 + ((r_o-r_i)*0.5)^2);
r_alpha2 = sqrt( (0.5*L2)^2 + ((r_o-r_i)*0.5)^2);

r1_start = [r_alpha1 ((r_o-r_i)*0.5) r_alpha1 L1/2 r_alpha1 ((r_o-r_i)*0.5) r_alpha1 L1/2];
r2_start = [r_alpha2 ((r_o-r_i)*0.5) r_alpha2 L2/2 r_alpha2 ((r_o-r_i)*0.5) r_alpha2 L2/2];
r3_start = r1_start;

center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
center2_start = [0 ((r_o-r_i)*0.5)+r_i];
center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];

alpha = atan(((r_o-r_i)*0.5)/(0.5*L1));
alpha2 = atan(((r_o-r_i)*0.5)/(0.5*L2));
angles1_start = [pi+alpha 3/2*pi -alpha 0 alpha 1/2*pi pi-alpha pi];
angles2_start = [pi+alpha2 3/2*pi -alpha2 0 alpha2 1/2*pi pi-alpha2 pi];
angles3_start = angles1_start;

x_vec = r1_start .* cos(angles1_start) + center1_start(1);
y_vec = r1_start .* sin(angles1_start) + center1_start(2);



%%

[x1_star, y1_star, l1, end_p1] = xy_star(length(angles1_start), 0, r1_start,...
    center1_start, angles1_start);
[x2_star, y2_star, l2, end_p2] = xy_star(length(angles1_start), end_p1, r2_start,...
    center2_start', angles2_start);
[x3_star, y3_star, l3, end_p3] = xy_star(length(angles1_start), end_p2, r3_start,...
    center3_start, angles3_start);


geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);


geom2d.addContour(x1_star, y1_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 2, 1:l1);

geom2d.addContour(x2_star, y2_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 3, 1:l2);

geom2d.addContour(x3_star, y3_star, s*ratio*[1 1 1 1.5 2 2 2 1.5], 4, 1:l3);

safetyplotGeom(geom2d, zeros(1,end_p3))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb );
fem.setDirichlet(4, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p3)';


fn_handle = @(fem) @(p) return_FEM_function_star_timer(fem, p, Lx_outer, Ly); 

minX = [ones(1,8)'*(-r1_start(1)); -5e-3; -r_i; -pi/8*ones(1,8)';ones(1,8)'*(-r2_start(1)); -5e-3; -r_i; -pi/8*ones(1,8)';ones(1,8)'*(-r3_start(1)); -5e-3; -r_i; -pi/8*ones(1,8)'];
%minX = ones(1,end_p3)'*(-5e-3);
%maxX = ones(1,end_p3)'*(10e-3);
maxX = [ones(1,8)'*(r1_start(1)); 5e-3; 10*r_i; pi/8*ones(1,8)';ones(1,8)'*(r2_start(1)); 5e-3; 10*r_i; pi/8*ones(1,8)';ones(1,8)'*(r3_start(1)); 5e-3; 10*r_i; pi/8*ones(1,8)'];
t1 = tic;

plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [1 1 1], 1);
[x, fval, iter, xHist, fHist, DfHist] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 15);
toc(t1)
save('StarOpti3x8' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist')
% 
