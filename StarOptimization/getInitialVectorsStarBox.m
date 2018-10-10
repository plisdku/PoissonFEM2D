function [phi, r_vec, s_vec] = getInitialVectorsStarBox(L,r_i,r_o, N_vec)

   
    if length(N_vec) == 1
        
        N_vec = [N_vec N_vec N_vec N_vec];
        
    end
    
    R = (r_o - r_i)*0.5;
    alpha = atan(R/(0.5*L));


    bottom_points = linspace(-L/2,L/2,N_vec(1)+2);
    right_points = linspace(-R,R,N_vec(2)+2);
    top_points = linspace(L/2,-L/2,N_vec(3)+2);
    left_points = linspace(R,-R,N_vec(4)+2);

    bottom_angles   = 1.5*pi + atan(bottom_points/R);
    bottom         = bottom_angles(2:(end-1));
    right_angles    = atan(right_points/(0.5*L));
    right           = right_angles(2:(end-1));
    top_angles      = 0.5*pi + atan(top_points/R);
    top             = top_angles((end-1):-1:2);
    left_angles     = pi + atan(left_points/(0.5*L));
    left            = left_angles((end-1):-1:2);

    phi = [pi+alpha bottom -alpha right alpha top pi-alpha left]; 

    r_corner = sqrt( (0.5*L)^2 + R^2);
    r_bottom = R*sec(bottom-1.5*pi);
    r_right = 0.5*L*sec(right); 
    r_top = R*sec(top-0.5*pi);
    r_left = 0.5*L*sec(left-pi);

    r_vec = [r_corner r_bottom r_corner r_right r_corner r_top ...
        r_corner r_left];
    
    s_vec = s_vec_lin((N_vec + 2),2,1);

end

% x_vec = r_vec .* cos(phi); 
% y_vec = r_vec .* sin(phi); 
% 
% figure(22112)
% plot(x_vec, y_vec, 'xr')
% 
% figure(2212222)
% plot(cos(phi), sin(phi), 'xk')