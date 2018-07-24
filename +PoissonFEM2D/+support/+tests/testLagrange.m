import PoissonFEM2D.support.*

isClose_tol = @(a,b,tol) abs(a-b) < tol;
isClose = @(a,b) isClose_tol(a,b,1e-12);

%% Single node (degenerate case).
% Only one polynomial, == 1.

lp = lagrangePolynomials(0);
assert(isequal(size(lp), [1,1]));

assert(polyval(lp(1,:), 0) == 1);
assert(polyval(lp(1,:), 1) == 1);

%% Several nodes

nodes = [0 2 3 3.5];

lp = lagrangePolynomials(nodes);
assert(isequal(size(lp), [4, 4]));

for whichPoly = 1:4
    for whichNode = 1:4
        
        %fprintf('lp(%i) at %2.2f\n', whichPoly, nodes(whichNode));
        
        if whichPoly == whichNode
            expectation = 1.0;
        else
            expectation = 0.0;
        end
        
        assert(isClose(polyval(lp(whichPoly, :), nodes(whichNode)), expectation));
    end
end

%%

disp('Lagrange polynomial test PASSED');

