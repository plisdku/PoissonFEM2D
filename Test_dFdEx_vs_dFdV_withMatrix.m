%u = femProblem.uCartesian*1e3;
u = UU1;

Nx          = measNxy(1);%size(xCoarse,2);
Ny          = measNxy(2);%size(yCoarse,2);
Nz          = 3;

uu = zeros(Nx,Ny,Nz);
uu(:,:,1) = u;%u2_reshape;
uu(:,:,2) = u;%u2_reshape;

x_grid = linspace(measBox(1),measBox(3),measNxy(1));%xCoarse;
y_grid = linspace(measBox(2),measBox(4),measNxy(2));


z_grid = [-2, 2, 4];


xv0 = [-0.8 -0.8;1.4 -1.4; 2 2;1000000 1000000;0 0;0 0];
x_p = 14e-3;
y_p  = 0;
z_p  = 0;
vx_p = 1;
vy_p = 0;
vz_p = 0;

obj_weights = [1, 1, 0, 0, 0, 0];

%%
U = [uu(:,end:-1:2,:) uu(:,:,:)];
a = fliplr(eye(Nx));
A = [a(1:(end-1),:); eye(Ny)];
U_matrix = (A*u')';


figure(401)
clf
hold on
imagesc(x_grid,r_grid,U(:,:,1)')

figure
clf
hold on
imagesc(x_grid,r_grid,U_matrix')

U(:,:,1) = U_matrix; 
U(:,:,2) = U_matrix;
U(:,:,3) = U_matrix;
%%
%V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
r_grid = [-y_grid(end:-1:2) y_grid];

d_x = x_grid(2) - x_grid(1);
d_r = y_grid(2) - y_grid(1);
d_z = z_grid(2) - z_grid(1);


E_x         = -centeredDiff(U, 1) / d_x;
E_r         = -centeredDiff(U, 2) / d_r;
E_z         = -centeredDiff(U, 3) / d_z;


Nt = 100;

nParticle = 1;
%dGdEx_meas = 0*E_x(:,:,1);
n_charges = 0.1;
n_masses = 1;
%n_charges = 1./elementary_charge;
%n_masses  = 1./atomic_mass;
accelFunc = accelerationFunction( x_grid, r_grid, z_grid, ...
    n_charges, n_masses);
%Nt = 4*Nt;
ts = linspace(0, 0.45e-6, Nt);



[xv, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x, E_r, E_z));
%[xv2, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x, E_r, E_z));


hit_objective = @(x_v) hitObjective3D_wrap(...
            x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);

%G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
G = hit_objective(xv);

[ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
    
    
figure(401)
clf
hold on
imagesc(x_grid,r_grid,U(:,:,2)')
axis xy image
xlabel('x')
ylabel('r')
figure(402)
clf
hold on
imagesc(x_grid,r_grid,E_x(:,:,2)')
axis xy image
xlabel('x')
ylabel('r')


plot(xv(ix_x,1),xv(ix_y,1),'ok', 'LineWidth', 3)
%plot(xv2(ix_x,1),xv2(ix_y,1),'k', 'LineWidth', 3)
plot(xv(1,1),xv(ix_y(1),1),'rx')
%plot(xv2(1,1),xv2(ix_y(1),1),'rx')
plot(x_p,y_p,'bx')
plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
plot(measBox([1,3,3,1,1]), -measBox([2,2,4,4,2]), 'w--');


%%

%nParticle = 2;


% hit_objective = @(x_v) hitObjective3D_wrap(...
%             x_v, x_p, y_p, z_p, vx_p, vy_p, vz_p, obj_weights);
%V_Comsol = zeros(Nx, Ny, Nz, 3);
V_Comsol = cell(3,1);
V_Comsol{1} = E_x;
V_Comsol{2} = E_r; 
V_Comsol{3} = E_z;
%[dGdEx0, dGdEy0, dGdEz0, G0, xv00, DG, xv_dual] = VV_get_dual_E_final(n_charges, n_masses, E_x, E_y, E_z,  x_grid, y_grid, z_grid, xv0, nParticle, hit_objective);
[dGdEx_sum, dGdEy_sum, dGdEz_sum, dGdV, dGdV_xr, G_sum, xv_all, DG, xv_dual, Nt] ...
    = VV_get_dual_E_cylindrical(n_charges, n_masses, V_Comsol, ...
    x_grid, r_grid, z_grid, xv0, nParticle, hit_objective, ts);




figure(4040)
clf
imagesc(x_grid, r_grid, dGdEx_sum(:,:,2)');
axis xy image
hold on
title('dual E_x')


figure(501)
imagesc(x_grid,r_grid,dGdEy_sum(:,:,2)');
axis xy image
title('Dual E_y')

%%

figure(501)
imagesc(x_grid,r_grid,dGdV(:,:,2)')
axis xy image

xlabel('x')
ylabel('y')
title('dual V')


%%
figure(502)
imagesc(x_grid,y_grid,dGdV_xr')
axis xy image

xlabel('x')
ylabel('r')
title('dual V in x,r')

%%
dGdV_measU2 = 0*U(:,:,1); 

%G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
G = hit_objective(xv);
%delta = 1e-3;
UU2_d = (A*UU2')';
UU1_d = (A*UU1')';


for xx = 3:17%1:length(x_grid)
    for yy = 31:35%5:8%12:15%(5:8)%1:length(y_grid)
        
        fprintf('%i, %i\n', xx, yy);
        delta = UU2_d(xx,yy) - UU1_d(xx,yy);%1e-3;%(UU2(xx,yy)-UU1(xx,yy));
        U_meas = U;
        U_meas(xx,yy,2) = U_meas(xx,yy,2) + delta;
        %max(max(Ex2-E_x))
        
        E_x2         = -centeredDiff(U_meas, 1) / d_x;
        E_r2         = -centeredDiff(U_meas, 2) / d_r;
        E_z2         = -centeredDiff(U_meas, 3) / d_z;
        
        [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
        %[xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));

        %G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        %G2 = sqrt(1/nParticle*(hit_objective(xv3).^2));
        G2 = hit_objective(xv3);
        G2
        G-G2
        max(xv-xv3)
        max(xv2-xv4)
        dGdV_measU2(xx,yy) = dGdV_measU2(xx,yy) + (G2-G)/(delta);




    end
end

%%
delta_fields = UU2_d - UU1_d;
delta_F_adj = sum(sum(dGdV(:,:,1) .* delta_fields))
delta_F_meas = sum(sum(dGdV_measU2 .* delta_fields))

%%
UU2_d2(:,:,1) = UU2_d;
UU2_d2(:,:,2) = UU2_d;

E_x_test         = -centeredDiff(UU2_d2, 1) / d_x;
E_r_test         = -centeredDiff(UU2_d2, 2) / d_r;
E_z_test         = -centeredDiff(UU2_d2, 3) / d_z;

[xv_test, accel] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x_test, E_r_test, E_z_test));

G_test = hit_objective(xv_test);

delta_F = G_test - G
%%
ratio = dGdV(:,:,1)'./dGdV_meas';
sum(sum(ratio(5:8, 3:15)))
%%
figure(4004)
clf
imagesc(dGdV_measU2')
axis xy image

figure(500)
imagesc(dGdV(:,:,1)')
axis xy image
%%
figure(5000)
imagesc(dGdV(:,:,1)')
axis xy image

xlabel('x')
ylabel('y')
title('dual V')

%%
figure(5001)
imagesc(dGdV(:,:,1)'./dGdV_measU2')
title('Factor')
axis xy image

%%
dGdV_meas = 0*U(:,:,1); 

G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));

delta = 1e-3;


for xx = 3:15%1:length(x_grid)
    for yy = 5:8%1:length(y_grid)
        
        fprintf('%i, %i\n', xx, yy);

        U2 = U;
        U2(xx,yy,1) = U2(xx,yy,1) + delta;
        %max(max(Ex2-E_x))
        
        E_x2         = -centeredDiff(U2, 1) / d_x;
        E_r2         = -centeredDiff(U2, 2) / d_r;
        E_z2         = -centeredDiff(U2, 3) / d_z;
        
        [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
        [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
        max(xv-xv3)
        max(xv2-xv4)
        dGdV_meas(xx,yy) = dGdV_meas(xx,yy) + (G2-G)/(delta);




    end
end
%% 

figure(4003)
clf
imagesc(dGdV_meas(:,:,1)')
axis xy image



%%
dGdEx_meas = 0*E_x(:,:,1);

G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));
%G
delta = 1e-3;


for xx = 3:15%1:length(x_grid)
    for yy = 5:8%1:length(y_grid)
        
        fprintf('%i, %i\n', xx, yy);

        Ex2 = E_x;
        Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
        %max(max(Ex2-E_x))
        [xv3, accel3] = velocityVerlet3D(ts, xv0(:,1), accelFunc(Ex2, E_r, E_z));
        [xv4, accel4] = velocityVerlet3D(ts, xv0(:,2), accelFunc(Ex2, E_r, E_z));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
%         if G == G2
%             keyboard
%         end 
        max(abs(xv-xv3))
        max(abs(xv2-xv4))
        dGdEx_meas(xx,yy) = dGdEx_meas(xx,yy) + (G2-G)/(delta);




    end
end


%%

figure(403)
clf
hold on
imagesc(x_grid,r_grid,dGdEx_meas')
plot(xv4(ix_x,1),xv4(ix_y,1),'ok', 'LineWidth', 3)

axis xy image

%% 

figure(500)
clf
hold on
plot(ts, accel2(:,1))
plot(ts, accel4(:,1))
plot(ts, accel(:,1),'o')
plot(ts, accel3(:,1),'x')


%%
dGdV_measU2 = 0*U(:,:,1); 

G = sqrt(1/nParticle*(hit_objective(xv).^2 + hit_objective(xv2).^2));

%delta = 1e-3;


for xx = 2:16%1:length(x_grid)
    for yy = 4:9%1:length(y_grid)
        
        fprintf('%i, %i\n', xx, yy);
        delta = 1e-3;%(UU2(xx,yy)-UU1(xx,yy));
        U_meas = U;
        U_meas(xx,yy,1) = U_meas(xx,yy,1) + delta;
        %max(max(Ex2-E_x))
        
        E_x2         = -centeredDiff(U_meas, 1) / d_x;
        E_r2         = -centeredDiff(U_meas, 2) / d_r;
        E_z2         = -centeredDiff(U_meas, 3) / d_z;
        
        [xv3, accel2] = velocityVerlet3D(ts, xv0(:,1), accelFunc(E_x2, E_r2, E_z2));
        [xv4, accel2] = velocityVerlet3D(ts, xv0(:,2), accelFunc(E_x2, E_r2, E_z2));

        G2 = sqrt(1/nParticle*(hit_objective(xv3).^2 + hit_objective(xv4).^2));
        G2
        G-G2
        max(xv-xv3)
        max(xv2-xv4)
        dGdV_measU2(xx,yy) = dGdV_measU2(xx,yy) + (G2-G)/(delta);




    end
end
