%% Test the sensitivity of the horrible determinant formula


f = @(A) sqrt(det(A'*A));

dfdA = @(A) f(A)*A*inv(A'*A);

%%

A = [1,2; 3,4];
dfda_calc = A;

delta = 1e-6;

for ii = 1:numel(A)
    
    A2 = A;
    A2(ii) = A2(ii) + delta;
    
    dfda_calc(ii) = (f(A2)-f(A))/delta;
    
end

fprintf('Numerical dfdA:\n');
disp(dfda_calc);
fprintf('Analytical dfdA:\n');
disp(dfdA(A));