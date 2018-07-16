%% Test cases for BoxVector Function 

V1 = [0, 0];
V2 = [2, 0];
V3 = [2, 4];
V4 = [0, 4];

Verticies = [V1; V2; V3; V4];
Verticies_P = [V1; V2; V3; V4; V1];

figure()
hold on
plot(Verticies_P(:,1), Verticies_P(:,2),'-o')

xlim([min(Verticies(:,1)) max(Verticies(:,1))] + [-2 2])
ylim([min(Verticies(:,2)) max(Verticies(:,2))] + [-2 2])

%%

[Box_x, Box_y] = BoxVector(Verticies, [2,2,2,2]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%%

[Box_x, Box_y] = BoxVector(Verticies, [3,3,3,3]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%%
[Box_x, Box_y] = BoxVector(Verticies, [5,5,5,5]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%%



[Box_x, Box_y] = BoxVector(Verticies, [5,2,5,8]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%%
V1 = [0, 0];
V2 = [2, 2];
V3 = [2, 4];
V4 = [0, 4];

Verticies = [V1; V2; V3; V4];
[Box_x, Box_y] = BoxVector(Verticies, [5,3,8,12]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%%

V1 = [0, 0];
V2 = [2, 2];
V3 = [3, 5];
V4 = [0, 4];

Verticies = [V1; V2; V3; V4];
[Box_x, Box_y] = BoxVector(Verticies, [5,3,8,12]);

figure()
plot([Box_x Box_x(1)],[Box_y Box_y(1)],'-o')

xlim([min(Box_x(:)) max(Box_x(:))] + [-2 2])
ylim([min(Box_y(:)) max(Box_y(:))] + [-2 2])

%% Test with function handle for s-vector

D_outer = Ly - d;

Ly = 45e-3;

D12 = 5e-3;
D23 = D12;
L1 = 2e-3;
L2 = 5e-3;
L3 = L1;

r_i = 3e-3;
r_o = 11e-3;

E2_V1 = [-L2/2, r_i];
E2_V2 = [L2/2, r_i];
E2_V3 = [L2/2, r_o];
E2_V4 = [-L2/2, r_o];
E2 = [E2_V1; E2_V2; E2_V3; E2_V4];
[E2_x, E2_y, s1_vec] = BoxVector(E2, [12 12 12 12], @(N) s_vec_lin(N,1,2));
l2 = length(E2_x);
assert(isequal(size(s1_vec),[1 l2]), 's_vec not the right length')

%%
[E2_x, E2_y, s1_vec] = BoxVector(E2, [2 2 2 2], @(N) s_vec_lin(N,1,2));
l2 = length(E2_x);
assert(isequal(s1_vec,[2 2 1 1]), 's_vec not correct (2 points per contour)')