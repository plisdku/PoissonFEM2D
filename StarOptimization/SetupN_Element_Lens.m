function [geom2d, end_p, minX, maxX] = SetupN_Element_Lens(N_elements, N_moveable, L_vec, D_vec, r_i, r_o, geom2d, s, ratio)

angles_start = cell(1,N_elements);
r_start     = angles_start;
s_vecs      = angles_start;

center_start = cell(1,N_elements);

for ii = 1:N_elements
   
    [angles_start{ii}, r_start{ii}, s_vecs{ii}] = ...
        getInitialVectorsStarBox(L_vec(ii), r_i, r_o, N_moveable(ii));
    

end

element_height = ((r_o-r_i)*0.5)+r_i;

center_start{2} = [-L_vec(1)/2-D_vec(1)-L_vec(2)/2 element_height];
center_start{1} = [0 element_height];
center_start{3} = [L_vec(1)/2+D_vec(2)+L_vec(3)/2 element_height];

for ii = 4:2:(N_elements-1)


        center_start{ii} = [center_start{ii-2}(1)-L_vec(ii-2)/2-D_vec(ii-1)-L_vec(ii)/2 ...
            element_height];
        center_start{ii+1} = [center_start{ii-1}(1)+L_vec(ii-1)/2+D_vec(ii)+L_vec(ii+1)/2 ...
            element_height];

end

if mod(N_elements,2) == 0
    center_start{N_elements} = [center_start{N_elements-2}(1)-L_vec(N_elements-2)/2-D_vec(N_elements-1)-L_vec(N_elements)/2 ...
            element_height];
end



x_star = cell(1,N_elements);
y_star = cell(1,N_elements);
l = cell(1,N_elements);
end_p = cell(1,N_elements);

for ii = 1:N_elements
    
    if ii ==1 
        
        start_p = 0;
    else
        start_p = end_p{ii-1};
    end
        
    [x_star{ii}, y_star{ii}, l{ii}, end_p{ii}] = xy_star(N_moveable(ii)*4+4,...
        start_p, r_start{ii}, center_start{ii}, angles_start{ii});
    
end 

minX = [];
maxX = [];

for ii = 1:N_elements
   
    geom2d.addContour(x_star{ii}, y_star{ii}, s*ratio*s_vecs{ii}, ii+1, 1:l{ii});
    minX = [minX; ones(1,(N_moveable(ii)*4+4))'*(-r_start{ii}(1)); -5e-3; -r_i; -pi/(N_moveable(ii)*4+4)*ones(1,(N_moveable(ii)*4+4))'];
    maxX = [maxX; ones(1,(N_moveable(ii)*4+4))'*(r_start{ii}(1)); 5e-3; 10*r_i; pi/(N_moveable(ii)*4+4)*ones(1,(N_moveable(ii)*4+4))'];
end



% 
% 
% 
%  minX = [ones(1,N_vec(1)*4+4)'*(-r1_start(1)); -5e-3; -r_i; -pi/(N_vec(1)*4+4)*ones(1,(N_vec(1)*4+4))';...
%         ones(1,(N_vec(2)*4+4))'*(-r2_start(1)); -5e-3; -r_i; -pi/(N_vec(2)*4+4)*ones(1,(N_vec(2)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))'];%...
%        % -50e3; -50e3;-50e3;-50e3;-50e3];
%     maxX = [ones(1,(N_vec(1)*4+4))'*(r1_start(1)); 5e-3; 10*r_i; pi/(N_vec(1)*4+4)*ones(1,(N_vec(1)*4+4))';...
%         ones(1,(N_vec(2)*4+4))'*(r2_start(1)); 5e-3; 10*r_i; pi/(N_vec(2)*4+4)*ones(1,(N_vec(2)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))';...
%         ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))'];%...
%         %50e3; 50e3; 50e3; 50e3; 50e3];






    

end 