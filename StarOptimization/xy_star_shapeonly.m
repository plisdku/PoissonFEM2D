function [x_star_func, y_star_func, l, end_p]= xy_star_shapeonly(N_points,...
    start_p, r_start, center_start, angles)

    pxr_end = start_p + N_points;
    pxc_end = pxr_end;
    px_end = pxc_end + N_points;
%     pyr_end = px_end + N_points; 
%     pyc_end = pyr_end + 2;
%     py_end = pyc_end + N_points;
    end_p = px_end;
    
    x_star_func = @(p) getStarContourVertex_x_warp(N_points,...
        p((start_p+1):pxr_end)' + r_start,...
        center_start,...
        angles + p((pxc_end+1):px_end)');
    
    y_star_func = @(p) getStarContourVertex_y_warp(N_points,...
        p((start_p+1):pxr_end)' + r_start,...
        center_start,...
        angles + p((pxc_end+1):px_end)');

    l = N_points; 

end