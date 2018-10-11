ScriptsPath = '/home/users/larstn/DesignOptimisation/ChargedParticles/Scripts';
VVPath = '/home/users/larstn/DesignOptimisation/ChargedParticles/Velocity-Verlet';
FEMPath = '/home/users/larstn/DesignOptimisation/ChargedParticles/PoissonFEM2D';
%rmpath(genpath('/Volumes/GoogleDrive/My Drive/Research/Design Optimization /Charged Particle Optics/Velocity-Verlet_Objective'), genpath('/Volumes/GoogleDrive/My Drive/Research/Design Optimization /Charged Particle Optics/Scripts'), genpath('/Volumes/GoogleDrive/My Drive/Research/Design Optimization /Charged Particle Optics/FEM'));
addpath(genpath(ScriptsPath), genpath(VVPath), genpath(FEMPath));
%%

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

% %%
% 
% 
% N_points = 8;
% 
% alpha_vec = pi/4 + linspace(0,2*pi,N_points+1);
% alpha_vec = alpha_vec(1:(end-1));
% 
% r_vec = 1*ones(1,N_points);
% 
% x_vec = r_vec .* cos(alpha_vec);
% y_vec = r_vec .* sin(alpha_vec);
% 
% x_center = 2;
% y_center = -2;
% 
% x_vec = x_vec + x_center;
% y_vec = y_vec + y_center;
% figure()
% plot(x_vec, y_vec, 'rx')
% 
% 
% [xs, ys] = getStarContourVertex(4,ones(1,4), [-3 5]);
% 
% figure()
% plot(xs,ys)

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

% % %%
% Nspline = 2;
% 
% [Ex1Spline, l11, s11_vec] =  xVectorBoxSpline(Ncontour_vec, E1x_func,...
% E1y_func, Nspline);
% [Ey1Spline, l11, s11_vec] =  yVectorBoxSpline(Ncontour_vec, E1x_func,...
% E1y_func, Nspline);
% [Ex2Spline, l21, s21_vec] =  xVectorBoxSpline(Ncontour_vec, E2x_func,...
% E2y_func, Nspline);
% [Ey2Spline, l12, s21_vec] =  yVectorBoxSpline(Ncontour_vec, E2x_func,...
% E2y_func, Nspline);
% [Ex3Spline, l13, s13_vec] =  xVectorBoxSpline(Ncontour_vec, E3x_func,...
% E3y_func, Nspline);
% [Ey3Spline, l13, s31_vec] =  yVectorBoxSpline(Ncontour_vec, E3x_func,...
% E3y_func, Nspline);


r1_start = sqrt( (0.5*L1)^2 + ((r_o-r_i)*0.5)^2)*ones(1,4);
r2_start = sqrt( (0.5*L2)^2 + ((r_o-r_i)*0.5)^2)*ones(1,4);
r3_start = r1_start;

center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
center2_start = [0 ((r_o-r_i)*0.5)+r_i];
center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];

alpha = atan(((r_o-r_i)*0.5)/(0.5*L1));
alpha2 = atan(((r_o-r_i)*0.5)/(0.5*L2));
angles1_start = [pi+alpha -alpha alpha pi-alpha];
angles2_start = [ pi+alpha2 -alpha2 alpha2 pi-alpha2];
angles3_start = angles1_start;

x_vec = r1_start .* cos(angles1_start) + center1_start(1);
y_vec = r1_start .* sin(angles1_start) + center1_start(2);

% figure(44)
% plot(x_vec, y_vec, 'rx')
% 
% x_vec = r2_start .* cos(angles2_start) + center2_start(1);
% y_vec = r2_start .* sin(angles2_start) + center2_start(2);
% 
% figure(45)
% plot(x_vec, y_vec, 'rx')
% 
% x_vec = r3_start .* cos(angles3_start) + center3_start(1);
% y_vec = r3_start .* sin(angles3_start) + center3_start(2);
% 
% figure(46)
% plot(x_vec, y_vec, 'rx')

%%

[x1_star, y1_star, l1, end_p1] = xy_star(4, 0, r1_start,...
    center1_start, angles1_start);
[x2_star, y2_star, l2, end_p2] = xy_star(4, end_p1, r2_start,...
    center2_start', angles2_start);
[x3_star, y3_star, l3, end_p3] = xy_star(4, end_p2, r3_start,...
    center3_start, angles3_start);


geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[15, ratio, ratio, ratio, 15, 25, 25, 25], 1, 1:8);
%geom2d.addContour(Ex1Spline, Ey1Spline, s*ratio*s11_vec, 2, 1:l11);
%geom2d.addContour(Ex2Spline, Ey2Spline, s*ratio*s21_vec, 3, 1:l21);
%geom2d.addContour(Ex3Spline, Ey3Spline, s*ratio*s31_vec, 4, 1:l13);
% p0 = zeros(1,end_p3)';
% p01 = zeros(1,end_p13)';
% figure(65)
% hold on
% plot(x1_star(p0),y1_star(p0),'xr')
% plot(x2_star(p0),y2_star(p0),'ko')
%plot(Ex1Spline(p01),Ey1Spline(p01),'br')

geom2d.addContour(x1_star, y1_star, s*ratio*[1 1 2 2], 2, 1:l1);

geom2d.addContour(x2_star, y2_star, s*ratio*[1 1 2 2], 3, 1:l2);

geom2d.addContour(x3_star, y3_star, s*ratio*[1 1 2 2], 4, 1:l3);

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

minX = [ones(1,4)'*(-r1_start(1)); -5e-3; -r_i; -pi/4*ones(1,4)';ones(1,4)'*(-r2_start(1)); -5e-3; -r_i; -pi/4*ones(1,4)';ones(1,4)'*(-r3_start(1)); -5e-3; -r_i; -pi/4*ones(1,4)'];
%minX = ones(1,end_p3)'*(-5e-3);
%maxX = ones(1,end_p3)'*(10e-3);
maxX = [ones(1,4)'*(r1_start(1)); 5e-3; 10*r_i; pi/4*ones(1,4)';ones(1,4)'*(r2_start(1)); 5e-3; 10*r_i; pi/4*ones(1,4)';ones(1,4)'*(r3_start(1)); 5e-3; 10*r_i; pi/4*ones(1,4)'];
t1 = tic;

plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [0 0 0], 1);

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 2);
toc(t1)
save('StarOpti444_bigBox' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alpha')
% 
