%% Contour points (oriented correctly) and mesh size parameters
close all
Lx = 60e-3;
Ly = 15e-3;
rAperture1 = 3e-3;
rAperture2 = 1.5e-3;
rAperture3 = 1e-3;
d = 4e-3;

isAxisymmetric = 1;

%%

N_field = 3;
N_geom = 2;
N_quad = N_field + isAxisymmetric; % there is a reason for this

%p0 = [0,0,0,0,0,0,0,0]';
s = 1e-3; % mesh scale

geom2d = ParameterizedGeometry2D();
% geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], s*[0.5, 4, 0.5, 4, 4, 4], 1, 1:6);
% geom2d.addContour(@(p) [-Lx+1e-3, -Lx+2e-3, -Lx+2e-3, -Lx+1e-3] + p([1 2 2 1])', @(p) [rAperture2, rAperture2, Ly-d, Ly-d] + p([3 3 4 4])', s*0.5, 2, 1:4);
% geom2d.addContour(@(p) [-Lx+3e-3, -Lx+4e-3, -Lx+4e-3, -Lx+3e-3] + p([5 6 6 5])', @(p) [rAperture2, rAperture2, Ly-d, Ly-d] + p([7 7 8 8])', s*0.5, 3, 1:4);
% geom2d.addContour(@(p) [Lx-2e-3, Lx-1e-3, Lx-1e-3, Lx-2e-3] + p([9 10 10 9])', @(p) [rAperture3, rAperture3, Ly-d, Ly-d] + p([11 11 12 12])', s*0.5, 4, 1:4);
% geom2d.addContour(@(p) [Lx-5e-3, Lx-4e-3, Lx-4e-3, Lx-5e-3] + p([13 14 14 13])', @(p) [rAperture3, rAperture3, Ly-d, Ly-d] + p([15 15 16 16])', s*0.5, 5, 1:4);
% 
% geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], s*[0.5, 4, 0.5, 4, 4, 4], 1, 1:6);
% geom2d.addContour(@(p) [-Lx+2e-3, -Lx+3e-3, -Lx+3e-3, -Lx+2e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d], s*0.5, 2, 1:4);
% geom2d.addContour(@(p) [-Lx+5e-3, -Lx+6e-3, -Lx+6e-3, -Lx+5e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d], s*0.5, 3, 1:4);
% geom2d.addContour(@(p) [-Lx+10e-3, -Lx+11e-3, -Lx+11e-3, -Lx+10e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d], s*0.5, 4, 1:4);
% geom2d.addContour(@(p) [-Lx+13e-3, -Lx+14e-3, -Lx+14e-3, -Lx+13e-3], @(p) [rAperture1, rAperture1, Ly-d, Ly-d], s*0.5, 5, 1:4);
l2 = 5e-3;
l1 = 2e-3;
l3 = 2e-3;
Vb = 29e3;
D1 = 3e-3;
L1 = 5e-3;
e_d = 7.5e-3;
Wd = 43e-3;

geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], s*[0.5, 4, 0.5, 4, 4, 4], 1, 1:6);
geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 2, 1:4);
geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 3, 1:4);
geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 4, 1:4);


fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
% fem.setDirichlet(2, @(p,x,y) 0.0 + x_old(1) + p(17));
% fem.setDirichlet(3, @(p,x,y) 1.0 + x_old(2) + p(18));
% fem.setDirichlet(4, @(p,x,y) 1.0 + x_old(3) + p(19));
% fem.setDirichlet(5, @(p,x,y) 0.0 + x_old(4) + p(20));
fem.setDirichlet(2, @(p,x,y) 0 + p(1));
fem.setDirichlet(3, @(p,x,y) Vb + p(2));
fem.setDirichlet(4, @(p,x,y) 0 + p(3));
%fem.setDirichlet(5, @(p,x,y) 0.0 + p(4));

fem.setFreeCharge(@(p,x,y) 0.0);

%%

p = [0,0,0]';


x0 = zeros(1,4)';
[~] = fem.instantiateProblem(p);

[femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);

xyGeomNodes = femProblem.poi.tnMesh.xyNodes;

measBox = [-55e-3, 0, 55e-3, 3e-3]; %[-0.5, -0.5, 0.5, 0.5];

fprintf('Forward solution... ');
xCoarse = linspace(-Lx, Lx, 200);
yCoarse = linspace(0, Ly, 200);
femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);

%u_cell{i}  = femProblem.uCartesian;  




Nt_vec = [100 200 400 800];
Nx_vec = [50 100 200 400];
Ny_vec = [20 40 80 160];

% Nt_vec = 200;
% Nx_vec = 200;
% Ny_vec = 80;

max_j = length(Nt_vec);
max_i = length(Nx_vec);

u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);


xv_cell = cell(max_i,max_j);

for i = 1:max_i

    measNxy = [Nx_vec(i), Ny_vec(i)];

    u_cartesian = femProblem.poi.tnMesh.rasterizeField(femProblem.u,...
        linspace(measBox(1),measBox(3),measNxy(1)), ...
        linspace(measBox(2),measBox(4),measNxy(2)));
    
    elementary_charge   = 1.60217662e-19;
    for j = 1:max_j
        
        Nt = Nt_vec(j);



        Nx          = measNxy(1);
        Ny          = measNxy(2);
        Nz          = 2;

        uu = zeros(Nx,Ny,Nz);
        uu(:,:,1) = u_cartesian;
        uu(:,:,2) = u_cartesian;


        E_center = 30*1e3*elementary_charge;
        ion_mass = 5.1477e-26;
        v_center = sqrt(2*E_center / ion_mass);
        delta_E = 10*elementary_charge; 
        v_spread = sqrt(2*delta_E / ion_mass);

        %angle = linspace(-2.5e-3,2.5e-3,3);


        %y_radius = linspace(0,0.1e-6,3);
        %velocities = linspace(v_center - v_spread, v_center + v_spread, 3);


        %x_pos = -52.5e-3;

%         [xx, yy, angles, vv] = ndgrid(x_pos, y_radius, angle, velocities);
%         angles = angles(:);
%         vv = vv(:);
%        xv0 = [xx(:)'; yy(:)'; zeros(length(xx(:)),1)'; (vv.*cos(angles))'; (vv.*sin(angles))'; zeros(length(xx(:)),1)'];
        
        xv0 = [x_pos; 0.5e-6; 0; v_center; 0; 0];
        
        x_grid = linspace(measBox(1),measBox(3),measNxy(1));
        y_grid = linspace(measBox(2),measBox(4),measNxy(2));
        z_grid = [-2, 0];

        x_p = 52.5e-3;
        y_p  = 0;
        z_p  = 0;
        vx_p = 1e-3;
        vy_p = 0;
        vz_p = 0;

        obj_weights = [1, 1, 0, 0, 0, 0];

        U = [uu(:,end:-1:2,:) uu(:,:,:)];

        r_grid = [-y_grid(end:-1:2) y_grid];

        d_x = x_grid(2) - x_grid(1);
        d_r = y_grid(2) - y_grid(1);
        d_z = z_grid(2) - z_grid(1);


        E_x         = -centeredDiff(U, 1) / d_x;
        E_r         = -centeredDiff(U, 2) / d_r;
        E_z         = -centeredDiff(U, 3) / d_z;


        end_time = (120e-3 ./ (v_center));
        end_time = end_time - 0.2*end_time;
        electron_mass       = 1.6605e-27;

        n_charges = -1;
        n_masses = ion_mass/electron_mass;

        accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
            n_charges, n_masses);

        ts = linspace(0, end_time, Nt);


        [xv, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x, E_r, E_z));

        xv_cell{i,j} = xv;
        
        figure(100); clf

        imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
        ax = axis;
        colormap orangecrush
        hold on
        femProblem.poi.tnMesh.plotMesh();

        geometry = fem.instantiatedGeom.geometry;
        plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
            [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)

        [ix_x, ix_y, ix_z, ~] =...
            get_Index3D(Nt);
        for ii = 1:size(xv,2)
            plot(xv(ix_x,ii),xv(ix_y,ii),'w', 'LineWidth', 1)
            plot(xv(ix_x(1),ii),xv(ix_y(1),ii),'rx')
        end
        plot(x_p,y_p,'bx')



        plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
        axis(ax)
        
        figure(110)
        plot(xv(ix_x,1),xv(ix_y,1),'k','LineWidth',1)
        xlabel('x')
        ylabel('y')       
    end
end
%%
xv_fin = xv_cell{max_i,max_j};
Nt_fin = size(xv_fin,1)./6;
[ix_x_fin, ix_y_fin, ~] = get_Index3D(Nt_fin);
x_fin = xv_fin(ix_x_fin(end),1);
y_fin = xv_fin(ix_y_fin(end),1);
%%
delta_pos = zeros(max_i,max_j);
for i=1:max_i
    for j=1:max_j
        xv_c = xv_cell{i,j};
        Nt_c = size(xv_c,1)./6;
        [ix_x_c, ix_y_c, ~] = get_Index3D(Nt_c);
        x_c = xv_c(ix_x_c(end),1);
        y_c = xv_c(ix_y_c(end),1);
        
        delta_pos(i,j) = sqrt( (x_fin-x_c).^2 + (y_fin - y_c).^2 );        
    end
end

%%
Nx_new = linspace(Nx_vec(1),...
    Nx_vec(end),100);
Nt_new = linspace(Nt_vec(1),Nt_vec(end),100);

[NxNx NtNt] = ndgrid(Nx_new, Nt_new);

delta_pos_lin = interpn(Nx_vec,Nt_vec,delta_pos,NxNx,NtNt);
figure(400)
clf
hold on
imagesc(Nx_new,Nt_new, delta_pos_lin')
axis image xy
[C,h] = contour(Nx_new, Nt_new, delta_pos_lin', 'w')
clabel(C,h, 'Color', 'w')
xlabel('Nx')
ylabel('Nt')
colorbar
colormap hot
title('Nx, Nt - VVConvergance')






