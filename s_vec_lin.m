function s_func = s_vec_lin(N, top, bottom)

    contour4 = linspace(top,bottom,N(3));
    s_func = [bottom*ones(1,N(1)-1) linspace(bottom,top,N(2)) top*ones(1,N(3)-2) contour4(1:(end-1))];

end