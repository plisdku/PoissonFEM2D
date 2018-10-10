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

[x1_star, y1_star, l1, end_p1] = xy_star_shapeonly(length(angles1_start), 0, r1_start,...
    center1_start, angles1_start);
[x2_star, y2_star, l2, end_p2] = xy_star_shapeonly(length(angles1_start), end_p1, r2_start,...
    center2_start', angles2_start);
[x3_star, y3_star, l3, end_p3] = xy_star_shapeonly(length(angles1_start), end_p2, r3_start,...
    center3_start, angles3_start);

geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);

geom2d.addContour(x1_star, y1_star, s*ratio*s1_vec, 2, 1:l1);

geom2d.addContour(x2_star, y2_star, s*ratio*s2_vec, 3, 1:l2);

geom2d.addContour(x3_star, y3_star, s*ratio*s3_vec, 4, 1:l3);

safetyplotGeom(geom2d, zeros(1,end_p3))


%%
% if isrow(p) == 1 
%         p = p';
% end
marksize = 10;
colors = jet(6); 
geometry = geom2d.evaluateGeometry(xHist(:,1));
figure(77)
clf
hold on
   plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
[geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o-','MarkerSize', marksize ,'color', colors(5,:), 'linewidth', 2)
geometry = geom2d.evaluateGeometry(xHist(:,5));
   plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
[geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o--', 'MarkerSize', marksize ,'color', colors(2,:), 'linewidth', 2)
geometry = geom2d.evaluateGeometry(xHist(:,10));
   plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
[geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o-.', 'MarkerSize', marksize ,'color', colors(6,:), 'linewidth', 2)
% geometry = geom2d.evaluateGeometry(xHist(:,19));
%    plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
% [geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o-', 'MarkerSize', marksize ,'color', colors(4,:), 'linewidth', 2)
geometry = geom2d.evaluateGeometry(xHist(:,25));
   plot([geometry.vertices(geometry.lines(:,1),1), geometry.vertices(geometry.lines(:,2),1)]', ...
[geometry.vertices(geometry.lines(:,1),2), geometry.vertices(geometry.lines(:,2),2)]', 'o:', 'MarkerSize', marksize ,'color', colors(1,:), 'linewidth', 2)

% 
%%
x1_p = @(p) x1_star(xHist(:,25))';
x2_p = @(p) x2_star(xHist(:,25))';
x3_p = @(p) x3_star(xHist(:,25))';

y1_p = @(p) y1_star(xHist(:,25))';
y2_p = @(p) y2_star(xHist(:,25))';
y3_p = @(p) y3_star(xHist(:,25))';
% 
% figure()
% hold on
% plot(x1_p(1:8), y1_p(1:8))
% plot(x2_p, y2_p)
% plot(x3_p, y3_p)
% 
% 
% [X Y Z] = cylinder(y1_p(1:8));
% 
% figure()
% surf(X,Y,Z)

%% paths
N = 1001;
phi = linspace(0,2*pi, N);
path1 = @(p,t) [0*cos(phi); 0*sin(phi); 0*ones(1,N)];
path2 = @(p,t) [center2_start(2)*cos(phi); center2_start(2)*sin(phi); center2_start(1)*ones(1,N)];
path3 = @(p,t) [center3_start(2)*cos(phi); center3_start(2)*sin(phi); center3_start(1)*ones(1,N)];
p1 = path1();
p2 = path2();
p3 = path3();
figure()
hold on 
plot3(p1(1,:), p1(2,:), p1(3,:))
plot3(p2(1,:), p2(2,:), p2(3,:))
plot3(p3(1,:), p3(2,:), p3(3,:))

%% stupid up and down arrows 
m = [ 0 1 0; 0 0 1; 1 0 0];
up = @(p) m*[cos(phi); sin(phi); zeros(1,N)];
right = @(p) m*[zeros(1,N); zeros(1,N); ones(1,N)];

tube1 = dmodel.ExtrudePath('Path', path1, 'X', x1_p, 'Y', y1_p, 'U', right, 'V', up, 'Closed', true);
tube2 = dmodel.ExtrudePath('Path', path1, 'X', x2_p, 'Y', y2_p, 'U', right, 'V', up,  'Closed', true);
tube3 = dmodel.ExtrudePath('Path', path1, 'X', x3_p, 'Y', y3_p, 'U', right, 'V', up,  'Closed', true);

figure(101); clf
hold on
b1 = tube1.mesh([]);
b2 = tube2.mesh([]);
b3 = tube3.mesh([]);
elec_color = [0.6 0.6 0.6];
h = flatPatch('Faces', b1.faces, 'Vertices', b1.patchVertices, 'FaceColor', elec_color, 'FaceAlpha', 0.8, 'EdgeAlpha', 0);
flatPatch('Faces', b2.faces, 'Vertices', b2.patchVertices, 'FaceColor', elec_color, 'FaceAlpha', 0.8, 'EdgeAlpha', 0);
flatPatch('Faces', b3.faces, 'Vertices', b3.patchVertices, 'FaceColor', elec_color, 'FaceAlpha', 0.8, 'EdgeAlpha', 0);

%set(gca,'View',[0 90])
xlabel('x')
ylabel('z')
zlabel('y')
camlight right
view(3)
camroll(90)
set(gca, 'FontSize', 16)
%camup([1 0 0])

