function [minX, maxX] = getBoundsStar(N_vec, r1_start, r2_start, r3_start, r_i)
 

    minX = [ones(1,N_vec(1)*4+4)'*(-r1_start(1)); -5e-3; -r_i; -pi/(N_vec(1)*4+4)*ones(1,(N_vec(1)*4+4))';ones(1,(N_vec(2)*4+4))'*(-r2_start(1)); -5e-3; -r_i; -pi/(N_vec(2)*4+4)*ones(1,(N_vec(2)*4+4))';ones(1,(N_vec(3)*4+4))'*(-r3_start(1)); -5e-3; -r_i; -pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))'];
    maxX = [ones(1,(N_vec(1)*4+4))'*(r1_start(1)); 5e-3; 10*r_i; pi/(N_vec(1)*4+4)*ones(1,(N_vec(1)*4+4))';ones(1,(N_vec(2)*4+4))'*(r2_start(1)); 5e-3; 10*r_i; pi/(N_vec(2)*4+4)*ones(1,(N_vec(2)*4+4))';ones(1,(N_vec(3)*4+4))'*(r3_start(1)); 5e-3; 10*r_i; pi/(N_vec(3)*4+4)*ones(1,(N_vec(3)*4+4))'];

end