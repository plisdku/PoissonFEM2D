%% Contour points (oriented correctly) and mesh size parameters
Lx = 60e-3;
Lx_outer = 20e-3;
Ly = 45e-3;
rAperture1 = 3e-3;
rAperture2 = 1.5e-3;
rAperture3 = 1e-3;
d = 34e-3;
Lx_outer2 = 120e-3;


isAxisymmetric = 1;

%%

N_field = 6;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this


s = 2.1e-3; % mesh scale
ratio = 0.3;
geom2d = ParameterizedGeometry2D();

l2 = 5e-3;
l1 = 2e-3;
l3 = 2e-3;
Vb = 29e3;
D1 = 3e-3;
L1 = 5e-3;
e_d = 7.5e-3;
Wd = 43e-3;
Wd2 = 164e-3;

        geom2d.addContour(@(p) [-Lx-Lx_outer-Lx_outer2, -Lx+e_d+Wd, 0,-Lx+e_d+Wd+l1+L1+l2+L1+l3, Lx+Lx_outer, Lx+Lx_outer, 0, -Lx-Lx_outer-Lx_outer2], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[5, ratio, ratio, ratio, 5, 8, 8, 8], 1, 1:8);
        geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd]+ p([1 2 2 1])', @(p) [D1, D1, Ly-d, Ly-d] + p([3 3 4 4])', s*ratio*[1, 1, 2, 2], 2, 1:4);
        geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1]+ p([5 6 6 5])', @(p) [D1, D1, Ly-d, Ly-d] + p([7 7 8 8])', s*ratio*[1, 1, 2, 2], 3, 1:4);
        geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1] + p([9 10 10 9])', @(p) [D1, D1, Ly-d, Ly-d] + p([11 11 12 12])', s*ratio*[1, 1, 2, 2], 4, 1:4);


fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
% fem.setDirichlet(2, @(p,x,y) 0.0 + x_old(1) + p(17));
% fem.setDirichlet(3, @(p,x,y) 1.0 + x_old(2) + p(18));
% fem.setDirichlet(4, @(p,x,y) 1.0 + x_old(3) + p(19));
% fem.setDirichlet(5, @(p,x,y) 0.0 + x_old(4) + p(20));
x_voltages = zeros(1,3);
fem.setDirichlet(2, @(p,x,y) 0 + x_voltages(1) );
fem.setDirichlet(3, @(p,x,y) Vb + x_voltages(2) );
fem.setDirichlet(4, @(p,x,y) 0 + x_voltages(3));
%fem.setDirichlet(5, @(p,x,y) 0.0 + p(4));

fem.setFreeCharge(@(p,x,y) 0.0);

%%


x0 = zeros(1,12)';


fn_handle = @(fem) @(p) return_FEM_functionobj(fem, p, Lx, Ly); 
minX = ones(1,12)'*(-5e-3);
% minX(1:16) = minX(1:16)*(-5e-3);
maxX = ones(1,12)'*(10e-3);
% maxX(1:16) = maxX(1:16)*(10e-3);
% 
% minX(17:20) = minX(17:20)*(-25);
% maxX(17:20) = maxX(17:20)*25;
%maxX = ones(1,4)'*25;
%minX = ones(1,4)'*-25;


[x2, fval2, iter2] = extremize12(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', @ plotfunc);

