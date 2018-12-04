
tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 10*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')



tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 15*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')




tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 20*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')





tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 30*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')





tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 50*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')





tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 75*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')



tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 100*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')






tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 200*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')




tic
import PoissonFEM2D.*
%%

Lx_outer = 3*60e-3;
isAxisymmetric = 1;

N_field = 5;
N_geom = 2;
N_quad = N_field + isAxisymmetric;

s = 1.4e-3;%4e-3; %1.3e-3; % mesh scale
ratio = 0.22;


Vb = 29e3;
D1 = 3e-3;

Ly = 3*45e-3;

D12 = 5e-3;
D23 = D12;

L1 = 2e-3;
L2 = 5e-3;
L3 = L2;


r_i = 3e-3;
r_o = 11e-3;

N_elements = 33;
N_moveable = 400*ones(1,N_elements);

L_vec = L1*ones(1,N_elements);
L_vec(1) = L2;
D_vec = D12*ones(1,N_elements-1);

geom2d = PoissonFEM2D.ParameterizedGeometry2D();
geom2d.addContour(@(p) [-Lx_outer -L2/2-D12-L1 0 L2/2+D23+L3 Lx_outer Lx_outer 0 -Lx_outer], @(p) [0, 0, 0, 0, 0, Ly, Ly, Ly], s*[4.5, ratio, ratio, ratio, 4.5, 9, 9, 9], 1, 1:8);
[geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio);

safetyplotGeom(geom2d, zeros(1,end_p{N_elements}))

fem = PoissonFEM2D.FEMInterface(geom2d, N_field, N_geom, N_quad, isAxisymmetric);
fem.setNeumann(1, @(p,x,y) 0.0);
fem.setDirichlet(2, @(p,x,y) 0 );
fem.setDirichlet(3, @(p,x,y) Vb);
for ii = 4:N_elements

    fem.setDirichlet(ii, @(p,x,y) 0 );
    
end
%fem.setDirichlet(5, @(p,x,y) 0 );
%fem.setDirichlet(6, @(p,x,y) 0 );

fem.setFreeCharge(@(p,x,y) 0.0);
fprintf('Number of moveable parameters is %d \n', end_p{N_elements})

%%


x0 = zeros(1,end_p{N_elements})';
%x0(end) = 1;


fn_handle = @(fem) @(p) return_FEM_function_inPlace_timer(fem, p, Lx_outer, Ly); 
%[minX, maxX] = getBoundsStar_shape_7elem([5 5 5], r1_start, r2_start, r3_start, r_i);

t1 = tic;
plot_func = @(xHistory,fHistory,DfHistory,max_x,alpha) plotfunc_star_alpha44(xHistory, fHistory, DfHistory, max_x, alpha, end_p{N_elements}, [5 5 5], 1);


%plot_func = @(xHistory,fHistory,DfHistory,max_x, alpha) plotfunc_star_voltage_alpha(xHistory, fHistory, DfHistory, max_x, alpha, end_p3, [5 5 5], 5);
disp('Setup time..')
toc
fprintf('\n')

[x, fval, iter, xHist, fHist, DfHist, times, alphas] = extremize_lars(fn_handle(fem), x0, 'Bounds', [minX, maxX], 'Callback', plot_func, 'MaxIter', 1);
toc(t1)
save('StarOpti_15Element_shapeOnly' , 'x' , 'fval', 'iter', 'xHist', 'fHist', 'DfHist', 'alphas')









