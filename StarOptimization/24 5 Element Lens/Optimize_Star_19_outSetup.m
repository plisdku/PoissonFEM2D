

import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;
geom2d = PoissonFEM2D.ParameterizedGeometry2D();


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

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

%L_vec = [L2, L1, L1, L1, L1, L1, L1];
D_vec = [D12, D12, D12, D12, D12, D12, D12];

r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 5*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);



% angles_start = cell(1,N_elements);
% r_start     = angles_start;
% s_vecs      = angles_start;
% 
% center_start = cell(1,N_elements);
% 
% for ii = 1:N_elements
%    
%     [angles_start{ii}, r_start{ii}, s_vecs{ii}] = ...
%         getInitialVectorsStarBox(L_vec(ii), r_i, r_o, N_moveable(ii));
%     
% 
% end
% 
% element_height = ((r_o-r_i)*0.5)+r_i;
% 
% center_start{2} = [-L_vec(2)/2-D_vec(1)-L_vec(1)/2 element_height];
% center_start{1} = [0 element_height];
% center_start{3} = [L_vec(2)/2+D_vec(2)+L_vec(3)/2 element_height];
% 
% for ii = 4:2:(N_elements-1)
%     
%     
%     center_start{ii} = [center_start{ii-2}(1)-L_vec(ii-2)/2-D_vec(ii-1)-L_vec(ii)/2 ...
%         element_height];
%     center_start{ii+1} = [center_start{ii-1}(1)+L_vec(ii-1)/2+D_vec(ii)+L_vec(ii+1)/2 ...
%         element_height];
%     
% end
% % 
% % 
% % [angles1_start, r1_start, s1_vec] = getInitialVectorsStarBox(L1,r_i,r_o,5);
% % [angles2_start, r2_start, s2_vec] = getInitialVectorsStarBox(L2,r_i,r_o,5);
% % [angles3_start, r3_start, s3_vec] = getInitialVectorsStarBox(L3,r_i,r_o,5);
% % [angles4_start, r4_start, s4_vec] = getInitialVectorsStarBox(L4,r_i,r_o,5);
% % [angles5_start, r5_start, s5_vec] = getInitialVectorsStarBox(L5,r_i,r_o,5);
% % [angles6_start, r6_start, s6_vec] = getInitialVectorsStarBox(L6,r_i,r_o,5);
% % [angles7_start, r7_start, s7_vec] = getInitialVectorsStarBox(L7,r_i,r_o,5);
% % 
% % center1_start = [-L2/2-D12-L1/2 ((r_o-r_i)*0.5)+r_i];
% % center2_start = [0 ((r_o-r_i)*0.5)+r_i];
% % center3_start = [L2/2+D23+L3/2  ((r_o-r_i)*0.5)+r_i];
% % center4_start = [-L1-D12-L2/2-D14-L4/2  ((r_o-r_i)*0.5)+r_i];
% % center5_start = [L2/2+L3+D23+L5/2+D35 ((r_o-r_i)*0.5)+r_i];
% % center6_start = [-L1-D12-L2/2-D14-L4-D46-L6/2  ((r_o-r_i)*0.5)+r_i];
% % center7_start = [L2/2+L3+D23+L5+D35+D57+L7/2 ((r_o-r_i)*0.5)+r_i];
% 
% 
% 
% x_star = cell(1,N_elements);
% y_star = cell(1,N_elements);
% l = cell(1,N_elements);
% end_p = cell(1,N_elements);
% 
% for ii = 1:N_elements
%     
%     if ii ==1 
%         
%         start_p = 0;
%     else
%         start_p = end_p{ii-1};
%     end
%         
%     [x_star{ii}, y_star{ii}, l{ii}, end_p{ii}] = xy_star(N_moveable(ii)*4+4,...
%         start_p, r_start{ii}, center_start{ii}, angles_start{ii});
%     
% end 
% % 
% % 
% % 
% % [x1_star, y1_star, l1, end_p1] = xy_star(length(angles1_start), 0, r1_start,...
% %     center1_start, angles1_start);
% % [x2_star, y2_star, l2, end_p2] = xy_star(length(angles1_start), end_p1, r2_start,...
% %     center2_start', angles2_start);
% % [x3_star, y3_star, l3, end_p3] = xy_star(length(angles1_start), end_p2, r3_start,...
% %     center3_start, angles3_start);
% % [x4_star, y4_star, l4, end_p4] = xy_star(length(angles1_start), end_p3, r4_start,...
% %     center4_start, angles4_start);
% % [x5_star, y5_star, l5, end_p5] = xy_star(length(angles1_start), end_p4, r5_start,...
% %     center5_start, angles5_start);
% % [x6_star, y6_star, l6, end_p6] = xy_star(length(angles1_start), end_p5, r6_start,...
% %     center6_start, angles6_start);
% % [x7_star, y7_star, l7, end_p7] = xy_star(length(angles1_start), end_p6, r7_start,...
% %     center7_start, angles7_start);
% 
% 
% %%
% 
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
% 
% 
% for ii = 1:N_elements
%    
%     geom2d.addContour(x_star{ii}, y_star{ii}, s*ratio*s_vecs{ii}, ii+1, 1:l{ii});
%     %geom2d.addContour(x_star{1}, y_star{1}, s*ratio*s_vecs{1}, 2, 1:l{2});
%   %  geom2d.addContour(x1_star, y1_star, s*ratio*s1_vec, 2, 1:l1);
%    % geom2d.addContour(x2_star, y2_star, s*ratio*s2_vec, 3, 1:l2);
% % 
%     % geom2d.addContour(x_star{3}, y_star{3}, s*ratio*s_vecs{3}, 4, 1:l{3});
% % 
%     % geom2d.addContour(x4_star, y4_star, s*ratio*s4_vec, 5, 1:l4);
% % 
% %     geom2d.addContour(x5_star, y5_star, s*ratio*s5_vec, 6, 1:l5);
% % 
% %     geom2d.addContour(x6_star, y6_star, s*ratio*s6_vec, 6, 1:l6);
% % 
% %     geom2d.addContour(x7_star, y7_star, s*ratio*s7_vec, 6, 1:l7);
% 
% end
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);

fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
fem.setDirichlet(4, @(p,x,y) 0 );
fem.setDirichlet(5, @(p,x,y) 0 );
fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);


%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p7, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 2);
toc(t1)
save('StarOpti_fiveElement_shapeOnlychromtighter' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')




