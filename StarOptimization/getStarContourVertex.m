function [x_star, y_star] = getStarContourVertex(N_points, r_vec, center, angles) 

    assert(length(r_vec) == N_points,...
        'Number of Radii and number of points must agree')

    %alpha_vec = -3*pi/4 + linspace(0,2*pi,N_points+1);
    %alpha_vec = alpha_vec(1:(end-1));
    
    x_vec = r_vec .* cos(angles);
    y_vec = r_vec .* sin(angles);
    
    x_star = x_vec + center(1);
    y_star = y_vec + center(2);

end