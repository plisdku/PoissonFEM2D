import support.*

isClose_tol = @(a,b,tol) abs(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Test differentiation with diffPoly

assert(isequal(diffPoly([1 0]), 1));
assert(isequal(diffPoly([2 1 1]), [4 1]));

disp('Polynomial differentiation PASSED');

%% Test indefinite integration with intPoly

assert(isequal(intPoly(1), [1 0]));
assert(isequal(intPoly([4 1]), [2 1 0]));

%% Test definite integration with intPoly

assert(isClose( intPoly([4 1], 0, 1), ...
    polyval(intPoly([4 1]), 1) - polyval(intPoly([4 1]), 0)));

%% Done.

disp('Polynomial integration PASSED')

