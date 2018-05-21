tt = linspace(0,2*pi,101);
%tt = tt(1:end-1);

R = 1;

x = geometry.vertices(:,1);
y = geometry.vertices(:,2);
x_n = [1 2 2 1]';
y_n = [1 1 2 2]';
%pathFn = @(p) [0*tt; R*sin(tt); -R*cos(tt)];
pathFn = @(p) repmat([0 0 0]',1, length(tt)); 
uFn = @(p) repmat([1 0 0]', 1, length(tt));
vFn = @(p) [0*tt; R*sin(tt); -R*cos(tt)];
%vFn = @(p) repmat([0 1 0]',1, length(tt));
%vFn = @(p) [R*cos(tt), R*sin(tt), 0*tt];
%tube1 = dmodel.ExtrudePath('X', @(p) x_n, 'Y',@(p) y_n, 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
tube1 = dmodel.ExtrudePath('X', @(p) x(7:10), 'Y',@(p) y(7:10), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
outer_tube = dmodel.ExtrudePath('X', @(p) x([1,3,4,6]), 'Y',@(p) y([1,3,4,6]), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
tube2 = dmodel.ExtrudePath('X', @(p) x(11:14), 'Y',@(p) y(11:14), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
tube3 = dmodel.ExtrudePath('X', @(p) x(15:18), 'Y',@(p) y(15:18), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
tube4 = dmodel.ExtrudePath('X', @(p) x(19:22), 'Y',@(p) y(19:22), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);
tube5 = dmodel.ExtrudePath('X', @(p) x(23:26), 'Y',@(p) y(23:26), 'U', uFn, 'V', vFn, 'Path', pathFn, 'Closed', true);

figure(101); clf
hold on
b1 = tube1.mesh([]);
b2 = tube2.mesh([]);
b3 = tube3.mesh([]);
b4 = tube4.mesh([]);
b5 = tube5.mesh([]);
bout = outer_tube.mesh([]);

flatPatch('Faces', b1.faces, 'Vertices', b1.patchVertices, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeAlpha', 0);
flatPatch('Faces', b2.faces, 'Vertices', b2.patchVertices, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeAlpha', 0);
flatPatch('Faces', b3.faces, 'Vertices', b3.patchVertices, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeAlpha', 0);
flatPatch('Faces', b4.faces, 'Vertices', b4.patchVertices, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeAlpha', 0);
flatPatch('Faces', b5.faces, 'Vertices', b5.patchVertices, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeAlpha', 0);

flatPatch('Faces', bout.faces, 'Vertices', bout.patchVertices, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

xlabel('x')
ylabel('y')
zlabel('z')
camlight right