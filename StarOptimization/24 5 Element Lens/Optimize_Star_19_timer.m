

import PoissonFEM2D.*
%%

Lx_outer = 60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;
geom2d = PoissonFEM2D.ParameterizedGeometry2D();


Vb = 29e3;
D1 = 3e-3;

Ly = 45e-3;

D12 = 5e-3;
D23 = D12;
D14 = D12;
D35 = D12;
D46 = D12;
D57 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L1;
L4 = L1;
L5 = L1;
L6 = L1;
L7 = L1;

r_i = 3e-3;
r_o = 11e-3;

[angles1_start, r1_start, s1_vec] = getInitialVectorsStarBox(L1,r_i,r_o,5);
[angles2_start, r2_start, s2_vec] = getInitialVectorsStarBox(L2,r_i,r_o,5);
[angles3_start, r3_start, s3_vec] = getInitialVectorsStarBox(L3,r_i,r_o,5);
[angles4_start, r4_start, s4_vec] = getInitialVectorsStarBox(L4,r_i,r_o,5);
[angles5_start, r5_start, s5_vec] = getInitialVectorsStarBox(L5,r_i,r_o,5);
[angles6_start, r6_start, s6_vec] = getInitialVectorsStarBox(L6,r_i,r_o,5);
[angles7_start, r7_start, s7_vec] = getInitialVectorsStarBox(L7,r_i,r_o,5);

center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
center2_start = [0 ((r_o-r_i)*0.5)+r_i];
center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];
center4_start = [-L1-D12-L2/2-D14-L4/2  ((r_o-r_i)*0.5)+r_i];
center5_start = [L2/2+L3+D23+L5/2+D35 ((r_o-r_i)*0.5)+r_i];
center6_start = [-L1-D12-L2/2-D14-L4-D46-L6/2  ((r_o-r_i)*0.5)+r_i];
center7_start = [L2/2+L3+D23+L5+D35+D57+L7/2 ((r_o-r_i)*0.5)+r_i];



[x1_star, y1_star, l1, end_p1] = xy_star(length(angles1_start), 0, r1_start,...
    center1_start, angles1_start);
[x2_star, y2_star, l2, end_p2] = xy_star(length(angles1_start), end_p1, r2_start,...
    center2_start', angles2_start);
[x3_star, y3_star, l3, end_p3] = xy_star(length(angles1_start), end_p2, r3_start,...
    center3_start, angles3_start);
[x4_star, y4_star, l4, end_p4] = xy_star(length(angles1_start), end_p3, r4_start,...
    center4_start, angles4_start);
[x5_star, y5_star, l5, end_p5] = xy_star(length(angles1_start), end_p4, r5_start,...
    center5_start, angles5_start);
[x6_star, y6_star, l6, end_p6] = xy_star(length(angles1_start), end_p5, r6_start,...
    center6_start, angles6_start);
[x7_star, y7_star, l7, end_p7] = xy_star(length(angles1_start), end_p6, r7_start,...
    center7_start, angles7_start);


%%

geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);

geom2d.addContour(x1_star, y1_star, s*ratio*s1_vec, 2, 1:l1);

geom2d.addContour(x2_star, y2_star, s*ratio*s2_vec, 3, 1:l2);

geom2d.addContour(x3_star, y3_star, s*ratio*s3_vec, 4, 1:l3);

geom2d.addContour(x4_star, y4_star, s*ratio*s4_vec, 5, 1:l4);

geom2d.addContour(x5_star, y5_star, s*ratio*s5_vec, 6, 1:l5);

geom2d.addContour(x6_star, y6_star, s*ratio*s6_vec, 6, 1:l6);

geom2d.addContour(x7_star, y7_star, s*ratio*s7_vec, 6, 1:l7);


safetyplotGeom(geom2d, zeros(1,end_p7))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
fem.setDirichlet(4, @(p,x,y) 0 );
fem.setDirichlet(5, @(p,x,y) 0 );
fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p7)';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_star_timer(fem, p, Lx_outer, Ly); 
[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p7, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 25);
toc(t1)
save('StarOpti_fiveElement_shapeOnlychromtighter' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')




