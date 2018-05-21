%% Input parameters for this test case
%clc 
%%
% Setting natural constants
% in case that they are needed
function [VV, Nt] = ElectronSetupdemo2DobjVV(u, measBox, measNxy)


elementary_charge   = 1.60217662e-19;
%electron_mass       = 9.1093856e-31;
%atomic_mass = 1.6605e-27;

Nx          = measNxy(1);%size(xCoarse,2);
Ny          = measNxy(2);%size(yCoarse,2);
Nz          = 2;

%Nt = 100;


%u2_reshape = reshape(u2,[10, 10]);
%uu = zeros(Nx,Ny,Nz);
%uu(:,:,1) = u;%u2_reshape;
%uu(:,:,2) = u;%u2_reshape;
%uu(:,:,3) = u;
% E_center = 30*1e3*elementary_charge;
% ion_mass = 5.1477e-26;
% v_center = sqrt(2*E_center / ion_mass);
% delta_E = 10*elementary_charge; 
% v_spread = sqrt(2*delta_E / ion_mass);
% 
% angle = linspace(-2.5e-3,2.5e-3,3);
% % cos_angles = cos(angles);
% % sin_angles = sin(angles);
% % v_x_center = v_center * ;
% % v_y_center = v_center * sin(2.5);
% 
% y_radius = linspace(0,0.1e-6,3);
% velocities = linspace(v_center - v_spread, v_center + v_spread, 3);
% % 
% % v_y_vec = [v_y_center 0 -v_y_center];
% % v_x_vec = [v_x_center v_center v_x_center];
% x_pos = -52.5e-3;
% 
% [xx, yy, angles, vv] = ndgrid(x_pos, y_radius, angle, velocities);
% angles = angles(:);
% vv = vv(:);
% xv0 = [xx(:)'; yy(:)'; zeros(length(xx(:)),1)'; (vv.*cos(angles))'; (vv.*sin(angles))'; zeros(length(xx(:)),1)'];


  Nr = 100;
        
  E_center = 30*1e3*elementary_charge;
ion_mass = 5.1477e-26;
v_center = sqrt(2*E_center / ion_mass);
delta_E = 10*elementary_charge; 
v_spread = sqrt(2*delta_E / ion_mass);

angle = linspace(-2.5e-3,2.5e-3,Nr);
% cos_angles = cos(angles);
% sin_angles = sin(angles);
% v_x_center = v_center * ;
% v_y_center = v_center * sin(2.5);

y_radius = 0;
velocities = v_center;%linspace(v_center - v_spread, v_center + v_spread, 3);
% 
% v_y_vec = [v_y_center 0 -v_y_center];
% v_x_vec = [v_x_center v_center v_x_center];
%x_pos = -52.5e-3;
x_pos = -173.5e-3;

[xx, yy, angles, vv] = ndgrid(x_pos, y_radius, angle, velocities);
angles = angles(:);
vv = vv(:);
xv0 = [xx(:)'; yy(:)'; zeros(length(xx(:)),1)'; (vv.*cos(angles))'; (vv.*sin(angles))'; zeros(length(xx(:)),1)'];






%xv0 = [x_pos; 0.1e-6; 0; v_center; 0;0];

%x_grid = linspace(-1, 1, Nx)*(elementary_charge/atomic_mass);
%y_grid = linspace(0, 1, round(0.5*Ny))*(elementary_charge/atomic_mass);
%z_grid = linspace(-1, 1, Nz)*(elementary_charge/atomic_mass);
x_grid = linspace(measBox(1),measBox(3),measNxy(1));%xCoarse;
y_grid = linspace(measBox(2),measBox(4),measNxy(2));%yCoarse;
%x_grid = linspace(-8.5e-3,14.5e-3,10);
%y_grid = linspace(0,1.6e-3,10);
z_grid = [-2, 0];

%xv0 = [-8e-3 -8e-3; 1e-3 -1e-3; -1e-3 -1e-3;10000 10000;0 0; 0 0];
%xv0 = [-10 -10 -10;0.2 0.15 0.1; -0.5 -0.5 -0.5; ;  0 0 0;  0 0 0]*1e-3;
x_p = 52.5e-3;
y_p  = 0;
z_p  = 0;
vx_p = 1e-3;
vy_p = 0;
vz_p = 0;

obj_weights = [1, 1, 0, 0, 0, 0];

%n_charges = 1;
%n_masses = 1;
%nParticle = 2;

%% 
% Creating test environment with given Input



%V = repmat(linspace(5,-5,11)',1,6,2);
%V = V + repmat(linspace(-2.5,2.5,6),11,1,2);
%%
%V = [V(:,end:-1:2,:) V];
% figure()
% imagesc(u')
% axis xy image
% figure()
% imagesc(uu(:,:,2)')
% axis xy image
% figure()
% imagesc(uu(:,end:-1:2,2)')
% axis xy image

%%
%U = [uu(:,end:-1:2,:) uu(:,:,:)];
%V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
%r_grid = [-y_grid(end:-1:2) y_grid];

%d_x = x_grid(2) - x_grid(1);
%d_r = y_grid(2) - y_grid(1);
%d_z = z_grid(2) - z_grid(1);

%V(:,:,2) = V(:,:,1);

%E_x         = -centeredDiff(U, 1) / d_x;
%E_r         = -centeredDiff(U, 2) / d_r;
%E_z         = -centeredDiff(U, 3) / d_z;

% E_x = ones(Nx,Ny,Nz);
% E_r = -ones(Nx,Ny,Nz);
% E_z = zeros(Nx,Ny,Nz);
% 
% E_2r = [E_r(:,end:-1:1,:) E_r(:,:,:)];
% r2_grid = [r_grid(end:-1:1) r_grid];

% %%
%     figure(20)
%     subplot(3,1,1)
%     imagesc(E_x(:,:,2)')
%     colorbar
%     subplot(3,1,2)
%     imagesc(E_r(:,:,2)')
%     colorbar
%     subplot(3,1,3)
%     imagesc(E_z(:,:,2)')
%     colorbar
%     figure(21)
%     imagesc(U(:,:,2)')
%     
%%    

%hit_objective = @(x_v) hitObjective3D_wrap(...
%            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);


%%
Nt = 1500;

end_time = (120e-3 ./ (v_center - v_spread))*2.5;
electron_mass       = 1.6605e-27;

%nParticle = 2;
%dGdEx_meas = 0*E_x(:,:,1);
n_charges = -1;
n_masses = ion_mass/electron_mass; %electron_mass / ion_mass;
%n_charges = 1./elementary_charge;
%n_masses  = 1./atomic_mass;
% accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
%     n_charges, n_masses);
%Nt = 4*Nt;
ts = linspace(0, end_time, Nt);


% [xv, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x, E_r, E_z));
% [xv2, accel] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_r, E_z));
% 
% %%
% [ix_x, ix_y, ix_z, ~] =...
%         get_Index3D(Nt);
% figure(401)
% clf
% hold on
% imagesc(x_grid,r_grid,U(:,:,2)')
% axis xy image
% xlabel('x')
% ylabel('r')
% 
% %%
% 
% %plot3(xv(ix_x,1),xv(ix_y,1)*cos(linspace(0,2*pi,11)),xv(ix_y,1)*sin(linspace(0,2*pi,11)),'k', 'LineWidth', 3)
% 
% %%
% plot(xv(ix_x,1),xv(ix_y,1),'k', 'LineWidth', 3)
% plot(xv2(ix_x,1),xv2(ix_y,1),'k', 'LineWidth', 3)
% plot(xv(1,1),xv(ix_y(1),1),'rx')
% plot(xv2(1,1),xv2(ix_y(1),1),'rx')
% plot(x_p,y_p,'wx')
% plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
% plot(measBox([1,3,3,1,1]), -measBox([2,2,4,4,2]), 'w--');


 %%

%nParticle = size(xv0,2);


hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
% hit_objective = @(x_v) hitObjective3D_wrap(...
%             x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
%V_Comsol = zeros(Nx, Ny, Nz, 3);
% V_Comsol = cell(3,1);
% V_Comsol{1} = E_x;
% V_Comsol{2} = E_r; 
% V_Comsol{3} = E_z;
fprintf('Setup for Trajectory and Dual Trajectory Calculations done')
clear particles
particles(size(xv0,2)) = Particle(xv0(4:6,size(xv0,2)),'Velocity',n_masses,n_charges,xv0(1:3,size(xv0,2)), 0, end_time, Nt);
for zz = 1:size(xv0,2)
     particles(zz) = Particle(xv0(4:6,zz),'Velocity',n_masses,n_charges,xv0(1:3,zz), 0, end_time, Nt);
end

hit_objective = @(x_v) hitObjective3D_wrap(...
    x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
xyz = [measBox(1), measBox(3); measBox(2),measBox(4); -2,0];
Nxy = [measNxy(1), measNxy(2), 2];
VV = cylindricalVelocity_Verlet(xyz, Nxy, particles, hit_objective, u);
VV = VV.calculateF_dF;
%[dGdEx0, dGdEy0, dGdEz0, G0, xv00, DG, xv_dual] = VV_get_dual_E_final(n_charges, n_masses, E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, nParticle, hit_objective);
%[dGdEx_sum, dGdEy_sum, dGdEz_sum, dGdV, dGdV_xr, G_sum, xv_all, DG, xv_dual, Nt] ...
    %= VV_get_dual_E_cylindrical(n_charges, n_masses, V_Comsol, ...
   % x_grid, r_grid, z_grid, xv0, nParticle, hit_objective, ts);
fprintf('Trajectory and Dual Trajectory Calculations done')

%%

% figure(4040)
% clf
% imagesc(x_grid, r_grid, dGdEx_sum(:,:)');
% axis xy image
% hold on
% title('dual E_x')
%plot(x_p,y_p,'rx')
%plot(xv(1:100,1),xv(101:200,1),'k')
%plot(xv3(1:100,1),xv3(101:200,1),'b')
%axis xy image
% colorbar
% title('dual E_x')
% figure(500)
% imagesc((dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_meas')
% % % dGdEx_sum./dGdEx_meas
% % % (dGdEx_sum(:,:,1)'-dGdEx_meas')./dGdEx_sum(:,:,1)'
%%
% figure(501)
% imagesc(x_grid,r_grid,dGdEy_sum(:,:)');
% axis xy image
% title('Dual E_y')
% 
% %%
% 
% figure(500)
% imagesc(x_grid,r_grid,dGdV(:,:)')
% axis xy image
% 
% xlabel('x')
% ylabel('y')
% title('dual V')

% %%
% figure(502)
% imagesc(x_grid,y_grid,dGdV_xr')
% axis xy image
% 
% xlabel('x')
% ylabel('r')
% title('dual V in x,r')
%% Reshape 
%dGdV_xr_oneslice = dGdV_xr(:,:,1);
%dGdV_xr_row = reshape(dGdV_xr_oneslice,1,[]);


% F = G_sum;
% DF = dGdV_xr;
% DFDEx = reflect_back(dGdEx_sum(:,:));
% DFDEy = reflect_back(dGdEy_sum(:,:));
% E_x_out = 0.5*reflect_back(E_x(:,:,2));

end


%%
% 
% G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
% %G
% delta = 1e-3;
% 
% 
% for xx = 1:(length(x_grid)-5)%1:length(x_grid)
%     for yy = 35:45%1:length(y_grid)
%         fprintf('%i, %i\n', xx, yy);
% 
%         Ex2 = E_x;
%         Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
%         %max(max(Ex2-E_x))
%         [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(Ex2, E_r, E_z));
%         [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(Ex2, E_r, E_z));
% 
%         G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
%         G2
%         G-G2
%         max(xv-xv3)
%         max(xv2-xv4)
%         dGdEx_meas(xx,yy) = dGdEx_meas(xx,yy) + (G2-G)/(delta);
% 
% 
% 
% 
%     end
% end
% 
% 
% %%
% G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
% %%
% delta = 1e-9;
% dGdV_meas = 0*U(:,:,1);
% 
% 
% for xx = 5:12%1:length(x_grid)
%     for yy = 10:(length(r_grid)-10)%1:length(y_grid)
%           fprintf('%i, %i\n', xx, yy);
%             U2 = U;
%             U2(xx,yy,1) = U2(xx,yy,1) + delta;
%             
%             E_x2         = -centeredDiff(U2, 1);
%             E_r2         = -centeredDiff(U2, 2);
%             E_z2         = -centeredDiff(U2, 3);
%             
%             [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
%             [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));
%             
%             G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
%             G2
%             G-G2
%             max(xv-xv3)
%             max(xv2-xv4)
%             dGdV_meas(xx,yy) = dGdV_meas(xx,yy) + (G2-G)/(delta);
%             
%             
%     end
% end
% 
% %%
% figure(607)
% imagesc(dGdV_meas')
% axis xy image
% 
% figure(608)
% imagesc(dGdV(:,:,1)')
% axis xy image
% 
% figure(609)
% imagesc(dGdV_meas'./dGdV(:,:,1)')
% axis xy image
% %%
% % 
% % dGdEx_meas2 = 0*E_x(:,:,1);
% % 
% % accelFunc = accelerationFunction( x_grid, y_grid, z_grid, ...
% %     n_charges, n_masses);
% % 
% % 
% % [xv2, accel] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_y, E_z));
% % 
% % 
% % G = sqrt(hit_objective(xv).^2);
% % 
% % delta = 1e-6;
% % 
% % 
% % for xx = 4:9
% %     for yy = 4:9
% %         fprintf('%i, %i\n', xx, yy);
% % 
% %         Ex2 = E_x;
% %         Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
% %         [xv2, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
% %         G2 = sqrt(hit_objective(xv2).^2);
% % 
% %         dGdEx_meas2(xx,yy) = dGdEx_meas2(xx,yy) + (G2-G)/delta;
% % 
% % 
% % 
% % 
% %     end
% % end
% % 
% 
% 
% % dGdEx_meas_sum = dGdEx_meas + dGdEx_meas2;
% 
% 
% figure(403)
% clf
% imagesc(dGdEx_meas')
% axis xy image
% 
% figure(40)
% clf
% %imagesc(dGdEx_meas')
% imagesc(x_grid, r_grid, dGdEx_meas');
% hold on
% plot(x_p,y_p,'rx')
% % plot(xv(1:100,1),xv(101:200,1),'k')
% % plot(xv3(1:100,1),xv3(101:200,1),'b')
% %axis xy image
% colorbar
% title('Meas') 
% % figure(42)
% % imagesc(x_grid, y_grid, dGdEx_meas2');
% % axis xy image
% % colorbar
% % title('Meas') 
% % figure(43)
% % imagesc(x_grid, y_grid, (dGdEx_meas' + dGdEx_meas2'));
% % axis xy image
% % colorbar
% % title('Meas') 



%dGdV_xr_row = 
%%
% 
% delta = 1e-3;
% dGdV_meas = 0*U(:,:,1);
% 
% 
% for xx = 1:11%1:length(x_grid)
%     for yy = 1:11%1:length(y_grid)
%           fprintf('%i, %i\n', xx, yy);
%             UV2 = U;
%             
%             accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
%     n_charges, n_masses);
% 
%             UV2(xx,yy,1) = UV2(xx,yy,1) + delta;
%             max(max(UV2(:,:,1)-V(:,:,1)))
%             E_x2         = -centeredDiff(UV2, 1) / d_x;
%             E_r2         = -centeredDiff(UV2, 2) / d_y;
%             E_z2         = -centeredDiff(UV2, 3) / d_z;
%             
%             [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
%             [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));
%             
%             G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
%             G2
%             G-G2
%             max(xv-xv3)
%             max(xv2-xv4)
%             dGdV_meas(xx,yy) = dGdV_meas(xx,yy) + (G2-G)/(delta);
%             
%             
%     end
% end
