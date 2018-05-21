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


s_vec = [4e-3 3e-3 2e-3 1.5e-3 1e-3 0.75e-3 0.5e-3];
%Nx_vec = [100 200 300 400 500];
%Ny_vec = [40 80 120 160 200];

u_cell = cell(7,1);

for i = 1:7
    s = s_vec(i);
    geom2d = ParameterizedGeometry2D();
    
    geom2d.addContour(@(p) [-Lx, 0, Lx, Lx, 0, -Lx], @(p) [0, 0, 0, Ly, Ly, Ly], s*[0.5, 4, 0.5, 4, 4, 4], 1, 1:6);
    geom2d.addContour(@(p) [-Lx+e_d+Wd, -Lx+e_d+Wd+l1, -Lx+e_d+Wd+l1, -Lx+e_d+Wd], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 2, 1:4);
    geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1+l2, -Lx+e_d+Wd+l1+L1], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 3, 1:4);
    geom2d.addContour(@(p) [-Lx+e_d+Wd+l1+L1+l2+L1, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1+l3, -Lx+e_d+Wd+l1+L1+l2+L1], @(p) [D1, D1, Ly-d, Ly-d], s*0.5, 4, 1:4);


    fem = FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
    fem.setNeumann(1, @(p,x,y) 0.0);
    fem.setDirichlet(2, @(p,x,y) 0 + p(1));
    fem.setDirichlet(3, @(p,x,y) Vb + p(2));
    fem.setDirichlet(4, @(p,x,y) 0 + p(3));

    fem.setFreeCharge(@(p,x,y) 0.0);

    disp('starting')
    x0 = zeros(1,4)';
    [~] = fem.instantiateProblem(p);

    [femProblem, dDirichlet_dp, dnx_dp, dny_dp] = fem.adjustProblem(p);

    xyGeomNodes = femProblem.poi.tnMesh.xyNodes;
    
    measBox = [-55e-3, 0, 55e-3, 3e-3]; %[-0.5, -0.5, 0.5, 0.5];
    measNxy = [200, 80];
    
    fprintf('Forward solution... ');
    xCoarse = linspace(-Lx, Lx, 200);
    yCoarse = linspace(0, Ly, 200);
    femProblem.solveCartesian(measBox(1:2), measBox(3:4), measNxy);

    u_cell{i}  = femProblem.uCartesian;   
    
    disp('finished')
    
    figure(100); clf
     
%      xCoarse = linspace(-Lx, Lx, 200);
%     yCoarse = linspace(0, Ly, 200);
    u = femProblem.poi.tnMesh.rasterizeField(femProblem.u, xCoarse, yCoarse);
    imagesc_centered(xCoarse, yCoarse, u'); axis xy image; colorbar
    ax = axis;
    colormap orangecrush
    hold on
    femProblem.poi.tnMesh.plotMesh();
    
    geometry = fem.instantiatedGeom.geometry;
    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
        [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'color', [0.8 0.8 0.8], 'linewidth', 2)
    
%     [ix_x, ix_y, ix_z, ~] =...
%         get_Index3D(Nt);
%     for i = 1:size(xv,2)
%     plot(xv(ix_x,i),xv(ix_y,i),'w', 'LineWidth', 1)
%     plot(xv(ix_x(1),i),xv(ix_y(1),i),'rx')
%     end
    %plot(x_p,y_p,'bx')
%     xlim([-Lx Lx])
%     ylim([0 Ly])


    %id = femProblem.iDirichlet;
%     id = ':';
%     quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'w-', 'linewidth', 2)
%     quiver(xyGeomNodes(id,1), xyGeomNodes(id,2), femProblem.dF_dxy(id,1), femProblem.dF_dxy(id,2), 'g-', 'linewidth', 1)
    plot(measBox([1,3,3,1,1]), measBox([2,2,4,4,2]), 'w--');
    %axis xy image
    axis(ax)
    
    figure(110)
    imagesc(femProblem.uCartesian')
    axis image xy
    colorbar 
    colormap orangecrush 
    
end
    
% fn_handle = @(fem) @(p) return_FEM_function(fem, p, Lx, Ly); 
% % minX = ones(1,20)';
% % minX(1:16) = minX(1:16)*(-5e-3);
% % maxX = ones(1,20)';
% % maxX(1:16) = maxX(1:16)*(10e-3);
% % 
% % minX(17:20) = minX(17:20)*(-25);
% % maxX(17:20) = maxX(17:20)*25;
% maxX = ones(1,4)'*25;
% minX = ones(1,4)'*-25;
% 
% 
% [x, fval, iter] = extremize12(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', @ plotfunc);
% 

%%
i_end = length(u_cell);
u_diff_sum = zeros(1,i_end);
RMS_err = u_diff_sum;
L1_err = u_diff_sum;
for i = 1:i_end
    
    u_diff = (u_cell{i_end} - u_cell{i});
    u_1 = u_cell{i_end};
    u_2 = u_cell{i};
    u_diff_sum(i) = norm(u_1(:) - u_2(:)) ;
    
    figure(40+i)
    imagesc(u_diff)
    axis xy image 
    RMS_err(i) = norm(u_diff(:))./norm(u_1(:));
    L1_err(i) = max(abs(u_diff(:)));
    
end
%%
figure(50)
clf
loglog(s_vec,u_diff_sum,'o-')
hold on

l = linspace(1,i_end,i_end);
loglog(s_vec,10e4*s_vec.^(1))
loglog(s_vec,10e4*s_vec.^2)
loglog(s_vec,10e4*s_vec.^3)
xlabel('s-mesh')
ylabel('Diff. Potential')
[a, ~]=polyfit(log2(s_vec(1:6)),log2(u_diff_sum(1:6)),1);
title(sprintf('Convergance Rate %0.2f', a(1)))

figure(60)
loglog(s_vec, RMS_err,'o-')
xlabel('s-mesh')
ylabel('RMS rel. err.')
axis tight

[a, ~]=polyfit(log2(s_vec(1:6)),log2(RMS_err(1:6)),1);

title(sprintf('N_{field} = %i Convergance Rate %0.2f', N_field, a(1)))



figure(61)
loglog(s_vec,L1_err,'o-')
xlabel('s-mesh')
ylabel('L_{inf} Err.')
axis tight

[a, ~]=polyfit(log2(s_vec(1:6)),log2(L1_err(1:6)),1);

title(sprintf('N_{field} = %i Convergance Rate %0.2f', N_field, a(1)))

%title(sprintf('N_{field} = %i', N_field))



%%

% %%
%  u_cell_1 = u_cell;
%  save('u_cell_first', 'u_cell_2')





