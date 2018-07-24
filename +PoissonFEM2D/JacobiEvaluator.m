classdef JacobiEvaluator
    
    properties
        N;      % NUMBER OF polynomials, order 0 through N-1
        coefficientArray;   % Indexed [n+1, coeff, alpha+1, beta+1]
        derivCoefficientArray; % Indexed [n+1, coeff, alpha+1, beta+1]
    end
    
    
    methods
        
        function obj = JacobiEvaluator(inN)
            obj.N = inN;
            obj.coefficientArray = zeros([inN, inN, inN, inN]);
            
            for ii = 0:inN-1
                for jj = 0:inN-1
                    obj.coefficientArray(:,:,ii+1,jj+1) = PoissonFEM2D.support.jacobiPolynomialsNormalized(inN,ii,jj);
                end
            end
            
            for ii = 0:inN-1
                for jj = 0:inN-1
                    for kk = 1:inN
                        obj.derivCoefficientArray(kk,:,ii+1,jj+1) = PoissonFEM2D.support.diffPoly(obj.coefficientArray(kk,:,ii+1,jj+1));
                    end
                end
            end
        end
        
        function y = evaluate(obj, nn, ii, jj, x)
            coeffs = obj.coefficientArray(nn+1,end-nn:end,ii+1,jj+1); % if I strip the leading zeros
            %coeffs = obj.coefficientArray(nn+1,:,ii+1,jj+1); % if I keep the leading zeros            
            y = PoissonFEM2D.paulyval(coeffs, x);
        end
        
        function dydx = evaluateDerivative(obj, nn, ii, jj, x)
            coeffs = obj.derivCoefficientArray(nn+1,end-nn+1:end,ii+1,jj+1); % strip the leading zeros
            dydx = PoissonFEM2D.paulyval(coeffs, x);
        end
    end
    
end