import PoissonFEM2D.*
%%
clear all

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

[E1x_func, E1y_func, end_p1] = createGeomFunction(E1_x, E1_y, 0);
[E2x_func, E2y_func, end_p2] = createGeomFunction(E2_x, E2_y, end_p1);
[E3x_func, E3y_func, end_p3] = createGeomFunction(E3_x, E3_y, end_p2);


%%
% 
% p1_add = (rand(end_p3,1) - 0.5)*0.75e-3;
% p2_add = (rand(end_p3,1) - 0.5)*0.75e-3;
% p3_add = (rand(end_p3,1) - 0.5)*0.75e-3;
% 
% p0 = zeros(1, end_p3)';
% p1 = p0 + p1_add;
% p2 = p1 + p2_add;
% p3 = p2 + p3_add;
% 
% 
% %%
% [Ex0, Ey0, ExSpline0, EySpline0] = createBoxSpline(Ncontour_vec, ...
%     E2x_func, E2y_func, p0, 1000);
% 
% [Ex1, Ey1, ExSpline1, EySpline1] = createBoxSpline(Ncontour_vec, ...
%     E2x_func, E2y_func, p1, 1000);
% 
% [Ex2, Ey2, ExSpline2, EySpline2] = createBoxSpline(Ncontour_vec, ...
%     E2x_func, E2y_func, p2, 1000);
% 
% [Ex3, Ey3, ExSpline3, EySpline3] = createBoxSpline(Ncontour_vec, ...
%     E2x_func, E2y_func, p3, 1000);
% 
% figure(567); clf
% hold on
% plot(ExSpline0,EySpline0,'rx')
% plot(ExSpline1,EySpline1,'bx')
% plot(ExSpline2,EySpline2,'gx')
% plot(ExSpline3,EySpline3,'kx')
% 
% figure(568); clf
% hold on 
% plot(Ex0, Ey0, 'ro')
% plot(ExSpline0, EySpline0, 'r-')
% plot(Ex1, Ey1, 'bo')
% plot(ExSpline1, EySpline1, 'b-')
% plot(Ex2, Ey2, 'go')
% plot(ExSpline2, EySpline2, 'g-')
% plot(Ex3, Ey3, 'ko')
% plot(ExSpline3, EySpline3, 'k-')

%%
Nspline = 10;

[Ex1Spline, l1, s1_vec] =  xVectorBoxSpline(Ncontour_vec, E1x_func,...
    E1y_func, Nspline);
[Ey1Spline, l1, s1_vec] =  yVectorBoxSpline(Ncontour_vec, E1x_func,...
    E1y_func, Nspline);
[Ex2Spline, l2, s2_vec] =  xVectorBoxSpline(Ncontour_vec, E2x_func,...
    E2y_func, Nspline);
[Ey2Spline, l2, s2_vec] =  yVectorBoxSpline(Ncontour_vec, E2x_func,...
    E2y_func, Nspline);
[Ex3Spline, l3, s3_vec] =  xVectorBoxSpline(Ncontour_vec, E3x_func,...
    E3y_func, Nspline);
[Ey3Spline, l3, s3_vec] =  yVectorBoxSpline(Ncontour_vec, E3x_func,...
    E3y_func, Nspline);


geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);

geom2d.addContour(Ex1Spline, Ey1Spline, s*ratio*s1_vec, 2, 1:l1);

geom2d.addContour(Ex2Spline, Ey2Spline, s*ratio*s2_vec, 3, 1:l2);

geom2d.addContour(Ex3Spline, Ey3Spline, s*ratio*s3_vec, 4, 1:l3);

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
% 
