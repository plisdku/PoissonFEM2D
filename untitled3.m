figure(); hold on

Nt = 300;
[ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
    
plot(xv_1(ix_x,:), xv_1(ix_y,:), 'k', 'LineWidth', 1)
plot(xv_2(ix_x,:), xv_2(ix_y,:), 'b', 'LineWidth', 1)
plot(xv_3(ix_x,:), xv_3(ix_y,:), 'y', 'LineWidth', 1)

Nt = 200;
[ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
    
plot(xv_4(ix_x,:), xv_4(ix_y,:), 'r', 'LineWidth', 1)

Nt = 400;
[ix_x, ix_y, ix_z, ~] =...
        get_Index3D(Nt);
plot(xv_5(ix_x,:), xv_5(ix_y,:), 'g', 'LineWidth', 1)
