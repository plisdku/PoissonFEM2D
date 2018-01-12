classdef BasisEvaluator1d
    
    properties
        jac1d@JacobiEvaluator;
        N;
    end
    
    methods
        
        function obj = BasisEvaluator1d(N)
            obj.jac1d = JacobiEvaluator(N);
            obj.N = N;
        end
        
        
        function y = evaluate(obj, ii, rr)
            y = obj.jac1d.evaluate(ii, 0, 0, rr);
        end
        
        function dydr = evaluateDerivative(obj, ii, rr)
            dydr = obj.jac1d.evaluateDerivative(ii, 0, 0, rr);
        end
        
        function V = vandermonde(obj, rr)
            
            if nargin == 1
                rr = support.nodes1d(obj.N);
            end
            
            V = zeros(numel(rr), obj.N);
            
            for cc = 1:obj.N
                V(:,cc) = obj.evaluate(cc-1, rr);
            end
        end
        
        function dVdr = vandermondeDerivative(obj, rr)
            if nargin == 1
                rr = support.nodes1d(obj.N);
            end
            
            dVdr = zeros(numel(rr), obj.N);
            
            for cc = 1:obj.N
                dVdr(:,cc) = obj.evaluateDerivative(cc-1, rr);
            end
        end
        
    end
    
    
end