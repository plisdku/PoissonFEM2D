import PoissonFEM2D.*
%%

Lx_outer = 60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.28e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.3;
geom2d = PoissonFEM2D.ParameterizedGeometry2D();


Vb = 29e3;
D1 = 3e-3;

Ly = 45e-3;

D12 = 5e-3;
D23 = D12;
L1 = 2e-3;
L2 = 5e-3;
L3 = L1;

r_i = 3e-3;
r_o = 11e-3;

[angles1_start, r1_start, s1_vec] = getInitialVectorsStarBox(L1,r_i,r_o,5);
[angles2_start, r2_start, s2_vec] = getInitialVectorsStarBox(L2,r_i,r_o,5);
[angles3_start, r3_start, s3_vec] = getInitialVectorsStarBox(L3,r_i,r_o,5);

center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
center2_start = [0 ((r_o-r_i)*0.5)+r_i];
center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];

[x1_star, y1_star, l1, end_p1] = xy_star(length(angles1_start), 0, r1_start,...
    center1_start, angles1_start);
[x2_star, y2_star, l2, end_p2] = xy_star(length(angles1_start), end_p1, r2_start,...
    center2_start', angles2_start);
[x3_star, y3_star, l3, end_p3] = xy_star(length(angles1_start), end_p2, r3_start,...
    center3_start, angles3_start);

geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);

geom2d.addContour(x1_star, y1_star, s*ratio*s1_vec, 2, 1:l1);

geom2d.addContour(x2_star, y2_star, s*ratio*s2_vec, 3, 1:l2);

geom2d.addContour(x3_star, y3_star, s*ratio*s3_vec, 4, 1:l3);

safetyplotGeom(geom2d, zeros(1,end_p3))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb + p(end_p3 + 1));
fem.setDirichlet(4, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p3 + 1)';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_star_chrom(fem, p, Lx_outer, Ly); 
[minX, maxX] = getBoundsStar_voltage([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;


plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 1);

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 25);
toc(t1)
save('StarOpti_' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')
% 
