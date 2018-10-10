function [ExSpline, l, s_vec] =  xVectorBoxSpline(Ncontrol, Ex_func,...
    Ey_func, Nspine)

    ExSpline = @(p) createBoxSpline_x_warp(Ncontrol, Ex_func,...
        Ey_func, p, Nspine);
    
    l = 4*Nspine -4;
    
    s_vec = s_vec_lin([Nspine Nspine Nspine Nspine] ,2,1);
    
    
end 