function [EySpline, l, s_vec] =  yVectorBoxSpline(Ncontrol, Ex_func,...
    Ey_func, Nspine)
    
    EySpline = @(p) createBoxSpline_y_warp(Ncontrol, Ex_func,...
        Ey_func, p, Nspine);

    
    l = 4*Nspine -4;
    
    s_vec = s_vec_lin([Nspine Nspine Nspine Nspine] ,2,1);

end 