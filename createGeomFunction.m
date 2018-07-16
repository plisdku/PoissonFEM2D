%% Function to return function handles to a function  that is able to move
% all points of the given input vector. 


% ATTENTION: if changes in this function are made/variants produce, ensure
% the right orientation of p-vector. It gets returned as a row vector. 


function [x_func, y_func, end_p] = createGeomFunction(x_vec, y_vec, last_p)

    x_func = @(p) x_vec + p((last_p+1):(last_p+length(x_vec)))';
    y_func = @(p) y_vec +...
        p((last_p+1+length(x_vec)):(last_p+length(x_vec)+length(y_vec)))';
    end_p = last_p+length(x_vec)+length(y_vec);

end