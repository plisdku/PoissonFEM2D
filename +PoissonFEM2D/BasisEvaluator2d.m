classdef BasisEvaluator2d
    
    properties
        jac1d%@JacobiEvaluator;
        N;
    end
    
    methods
        
        function obj = BasisEvaluator2d(N)
            obj.jac1d = PoissonFEM2D.JacobiEvaluator(2*N);
            obj.N = N;
        end
        
        
        function y = evaluate(obj, ii, jj, rr, ss)
            
            a = 2*(1+rr)./(1-ss) - 1;
            a(ss == 1) = 1.0;
            b = ss;
            
            p_x = obj.jac1d.evaluate(ii, 0, 0, a);
            p_y = obj.jac1d.evaluate(jj, 2*ii+1, 0, b);
            extraFactor = (1-b).^ii;
            
            y = sqrt(2.0)*p_x.*p_y.*extraFactor;
            
        end
        
        function [dydr, dyds, y] = evaluateDerivative(obj, ii, jj, rr, ss)
            
            ss(ss == 1.0) = 1 - 1e-9;
            
            a = 2*(1+rr)./(1-ss) - 1;
            b = ss;
            
            da_dr = 2./(1-ss);
            da_ds = 2*(1+rr)./(1-ss).^2;
            %db_dr = 0*rr;
            %db_ds = ones(size(s));
            
            p_r = obj.jac1d.evaluate(ii, 0, 0, a);
            p_s = obj.jac1d.evaluate(jj, 2*ii+1, 0, b);
            extraFactor = (1-b).^ii;
            
            dpr_da = obj.jac1d.evaluateDerivative(ii, 0, 0, a);
            dpr_dr = dpr_da .* da_dr;
            dpr_ds = dpr_da .* da_ds;
            
            dps_ds = obj.jac1d.evaluateDerivative(jj, 2*ii+1, 0, b); % .* db_ds;
            dextra_ds = -ii*(1-b).^(ii-1); % .* db_ds
            
            y = sqrt(2.0)*p_r.*p_s.*extraFactor;
            
            dydr = sqrt(2.0) * p_s .* extraFactor .* dpr_dr;
            dyds = sqrt(2.0) * (dpr_ds .* p_s .* extraFactor + ...
                p_r .* dps_ds .* extraFactor + ...
                p_r .* p_s .* dextra_ds);
        end
        
        function V = vandermonde(obj, rr, ss)
           
            [is,js] = ndgrid(0:obj.N-1,0:obj.N-1);
            pickThese = (is+js < obj.N);
            is = is(pickThese);
            js = js(pickThese);
            
            if nargin < 3
                [rr,ss] = PoissonFEM2D.support2d.nodes2d(obj.N);
            end
            
            V = zeros(numel(rr), numel(is));
            
            for nn = 1:length(is)
                V(:,nn) = obj.evaluate(is(nn), js(nn), rr, ss);
            end
        end
        
        function [dVdr, dVds] = vandermondeDerivative(obj, rr, ss)
            
            [is,js] = ndgrid(0:obj.N-1,0:obj.N-1);
            pickThese = (is+js < obj.N);
            is = is(pickThese);
            js = js(pickThese);
            
            if nargin < 3
                [rr,ss] = PoissonFEM2D.support2d.nodes2d(obj.N);
            end
            
            dVdr = zeros(numel(rr), numel(is));
            dVds = zeros(numel(rr), numel(is));
            
            for nn = 1:length(is)
                [dVdr(:,nn), dVds(:,nn)] = obj.evaluateDerivative(is(nn), js(nn), rr, ss);
            end
        end
        
    end
    
    
end