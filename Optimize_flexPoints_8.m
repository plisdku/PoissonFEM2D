tic
Lx = 60e-3;
Lx_outer = 20e-3;
Ly = 45e-3;
rAperture1 = 3e-3;
rAperture2 = 1.5e-3;
rAperture3 = 1e-3;
d = 34e-3;

isAxisymmetric = 1;

%%

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this

%p0 = [0,0,0,0,0,0,0,0]';
s = 1.3e-3; % mesh scale
ratio = 0.3;
geom2d = ParameterizedGeometry2D();

l2 = 5e-3;
l1 = 2e-3;
l3 = 2.5e-3;
Vb = 29e3;
D1 = 3e-3;
L1 = 5e-3;
L2 = 2e-3;
e_d = 7.5e-3;
Wd = 43e-3;

%Lx  = 80e-3;
D_outer = Ly - d;

Ly = 45e-3;

D12 = 5e-3;
D23 = D12;
L1 = 2e-3;
L2 = 5e-3;
L3 = L1;

r_i = 3e-3;
r_o = 11e-3;

E2_V1 = [-L2/2, r_i];
E2_V2 = [L2/2, r_i];
E2_V3 = [L2/2, r_o];
E2_V4 = [-L2/2, r_o];
E2 = [E2_V1; E2_V2; E2_V3; E2_V4];
[E2_x, E2_y] = BoxVector(E2, [8 8 8 8]);
l2 = length(E2_x);
l22 = 2*l2;

E1_V1 = [-L2/2, r_i] - [D12+L1 0];
E1_V2 = [-L2/2, r_i] - [D12 0];
E1_V3 = [-L2/2, r_o] - [D12 0];
E1_V4 = [-L2/2, r_o] - [D12+L1 0];
E1 = [E1_V1; E1_V2; E1_V3; E1_V4];
[E1_x, E1_y] = BoxVector(E1, [8 8 8 8]);
l1 = length(E1_x);
l12 = 2*l1;

E3_V1 = [L2/2, r_i] + [D23 0];
E3_V2 = [L2/2, r_i] + [D23+L3 0];
E3_V3 = [L2/2, r_o] + [D23+L3 0];
E3_V4 = [L2/2, r_o] + [D23 0];
E3 = [E3_V1; E3_V2; E3_V3; E3_V4];
[E3_x, E3_y] = BoxVector(E3, [8 8 8 8]);
l3 = length(E3_x);
l32 = 2*l3;

plotElectrodePoints(E1_x, E1_y, E2_x, E2_y, E3_x, E3_y)

[E1x_func, E1y_func, end_p1] = createGeomFunction(E1_x, E1_y, 0);
[E2x_func, E2y_func, end_p2] = createGeomFunction(E2_x, E2_y, end_p1);
[E3x_func, E3y_func, end_p3] = createGeomFunction(E3_x, E3_y, end_p2);


      %  geom2d.addContour(@(p) [-Lx-Lx_outer, -Lx+e_d+Wd, 0,-Lx+e_d+Wd+l1+L1+l2+L1+l3, Lx+Lx_outer+25e-3, Lx+Lx_outer+30e-3, 0, -Lx-Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
      %  geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd]- 2e-3+ p([1 2 3 4])', @(p) [D1, D1, Ly-d, Ly-d]+ p([5 6 7 8])', s*ratio*[1, 1, 2, 2], 2, 1:4);
      %  geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+l2+l3, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1]+ p([9 10 11 12])', @(p) [D1 + rAperture1, D1+ rAperture1, Ly-d, Ly-d] + p([13 14 15 16])', s*ratio*[1, 1, 2, 2], 3, 1:4);
      %  geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1+L2, -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2, -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2, -Lx+e_d+Wd+l1+L1+l2+L1+L2]+ p([17 18 19 20])', @(p) [D1, D1, Ly-d, Ly-d]+ p([21 22 23 24])', s*ratio*[1, 1, 2, 2], 4, 1:4);

        geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
        
        geom2d.addContour(E1x_func, E1y_func, s*ratio*[1, 1, 1, 1, 1, 1, 1, 1, 1.14, 1.28, 1.42, 1.57, 1.71, 1.86 , 2, 2, 2, 2, 2, 2, 2, 2, 1.86, 1.71, 1.57, 1.42, 1.28, 1.14], 2, 1:l1);
        geom2d.addContour(E2x_func, E2y_func, s*ratio*[1, 1, 1, 1, 1, 1, 1, 1, 1.14, 1.28, 1.42, 1.57, 1.71, 1.86 , 2, 2, 2, 2, 2, 2, 2, 2, 1.86, 1.71, 1.57, 1.42, 1.28, 1.14], 3, 1:l2);
        geom2d.addContour(E3x_func, E3y_func, s*ratio*[1, 1, 1, 1, 1, 1, 1, 1, 1.14, 1.28, 1.42, 1.57, 1.71, 1.86 , 2, 2, 2, 2, 2, 2, 2, 2, 1.86, 1.71, 1.57, 1.42, 1.28, 1.14], 4, 1:l3);
        
       safetyplotGeom(geom2d, zeros(1,end_p3))
%         geom2d.addContour(@(p) [-Lx-Lx_outer, -Lx+e_d+Wd, 0,-Lx+e_d+Wd+l1+L1+l2+L1+l3, Lx+Lx_outer+25e-3, Lx+Lx_outer+25e-3, 0, -Lx-Lx_outer],...
%            @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
%         geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+1/3*l1,-Lx+e_d+Wd+2/3*l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+2/3*l1, -Lx+e_d+Wd+1/3*l1, -Lx+e_d+Wd, -Lx+e_d+Wd, -Lx+e_d+Wd]- 2e-3 +p([1 2 3 4 5 6 7 8 9 10 11 12])',...
%             @(p) [D1, D1, D1, D1, D1 + 1/3*(D_outer - D1), D1 + 2/3*(D_outer - D1), D_outer, D_outer, D_outer, D_outer, D1 + 2/3*(D_outer - D1), D1 + 1/3*(D_outer - D1) ] + p([13 14 15 16 17 18 19 20 21 22 23 24])', s*ratio*[1, 1, 1, 1, 1.3, 1.6, 2, 2, 2, 2, 1.6, 1.3], 2, 1:12);
%         geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+1/3*l2, -Lx+e_d+Wd+l1+L1+2/3*l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+2/3*l2, -Lx+e_d+Wd+l1+L1+1/3*l2, -Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1]+ p([25 26 27 28 29 30 31 32 33 34 35 36])',...
%             @(p) [D1, D1, D1, D1, D1 + 1/3*(D_outer - D1), D1 + 2/3*(D_outer - D1), D_outer, D_outer, D_outer, D_outer, D1 + 2/3*(D_outer - D1), D1 + 1/3*(D_outer - D1) ]+ p([37 38 39 40 41 42 43 44 45 46 47 48])', s*ratio*[1, 1, 1, 1, 1.3, 1.6, 2, 2, 2, 2, 1.6, 1.3], 3, 1:12);
%         geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1+L2, -Lx+e_d+Wd+l1+L1+l2+L1+L2+1/3*l3, -Lx+e_d+Wd+l1+L1+l2+L1+L2+2/3*l3, -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2, -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2, -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2,...
%             -Lx+e_d+Wd+l1+L1+l2+L1+l3+L2,  -Lx+e_d+Wd+l1+L1+l2+L1+L2+2/3*l3, -Lx+e_d+Wd+l1+L1+l2+L1+L2+1/3*l3, -Lx+e_d+Wd+l1+L1+l2+L1+L2, -Lx+e_d+Wd+l1+L1+l2+L1+L2, -Lx+e_d+Wd+l1+L1+l2+L1+L2]+ p([49 50 51 52 53 54 55 56 57 58 59 60])',...
%             @(p) [D1, D1, D1, D1, D1 + 1/3*(D_outer - D1), D1 + 2/3*(D_outer - D1), D_outer, D_outer, D_outer, D_outer, D1 + 2/3*(D_outer - D1), D1 + 1/3*(D_outer - D1) ]+ p([61 62 63 64 65 66 67 68 69 70 71 72])', s*ratio*[1, 1, 1, 1, 1.3, 1.6, 2, 2, 2, 2, 1.6, 1.3], 4, 1:12);

fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb );
fem.setDirichlet(4, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p3)';


fn_handle = @(fem) @(p) return_FEM_functionobj(fem, p, Lx, Ly); 

minX = ones(1,end_p3)'*(-5e-3);
maxX = ones(1,end_p3)'*(10e-3);



[x, fval, iter, xHist, fHist, DfHist] = extremize12(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', @ plotfunc);

save('ShowcaseData' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist')

