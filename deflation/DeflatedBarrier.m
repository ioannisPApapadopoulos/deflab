classdef DeflatedBarrier < handle
    
    properties
        muStart = 1
        muEnd = 1e-10
        lbTrue
        ubTrue
        lbEnlarged
        ubEnlarged
        constrainedSet
        complementSet
        kmu = 0.5
        theta = 2
    end
    methods
        function dab = DeflatedBarrier(muStart, muEnd, lb, ub)
            dab.muStart = muStart;
            dab.muEnd = muEnd;
            dab.lbTrue = lb;
            dab.ubTrue = ub;
            dab.lbEnlarged = lb-(1e-5);
            dab.ubEnlarged = ub+(1e-5);           
        end
        
        function solution = solve(dab, x, residual, jacobian)
            mu = dab.muStart;
            tol = 1e10;
            indices = 1:length(x);
            dab.constrainedSet = union(find(dab.lbEnlarged(dab.lbEnlarged>-tol)),...
                find(dab.ubEnlarged(dab.ubEnlarged<tol)));
            dab.complementSet = setdiff(indices, dab.constrainedSet);
            solution = dab.subproblem(mu, jacobian, residual, x);
        end
        
        function solution = subproblem(dab, mu, x, residual, jacobian)
            solver = NonlinearSolver;
            while mu > dab.muEnd
                mubarrierresidual = @(x) dab.barrierResidual(mu, residual, x);
                mubarrierjacobian = @(x) dab.barrierJacobian(mu, jacobian, x);
                fprintf("\nConsidering mu = %e\n", mu)
                x = solver.ssls(mubarrierjacobian, mubarrierresidual, x, dab.lbTrue, dab.ubTrue);
                mu = dab.updateMu(mu);
            end
            fprintf("\nConsidering mu = 0\n")
            x = solver.ssls(jacobian, residual, x, dab.lbTrue, dab.ubTrue);
            solution = x;
        end
        
        function out = barrierResidual(dab, mu, residual, x)  
            out = zeros(length(x),1);
            evaluatedResidual = residual(x);
            out(dab.complementSet) = evaluatedResidual(dab.complementSet);
            out(dab.constrainedSet) = evaluatedResidual(dab.constrainedSet)...
                -mu*1./(x(dab.constrainedSet) - dab.lbEnlarged(dab.constrainedSet))...
                +mu*1./(dab.ubEnlarged(dab.constrainedSet) - x(dab.constrainedSet));
        end
 
        function out = barrierJacobian(dab, mu, jacobian, x)  
            out = zeros(length(x));
            evaluatedJacobian = jacobian(x);
            out(dab.complementSet,dab.complementSet) = evaluatedJacobian(dab.complementSet,dab.complementSet);
            out(dab.constrainedSet,dab.constrainedSet) = evaluatedJacobian(dab.constrainedSet,dab.constrainedSet)...
                +mu*diag(1./(x(dab.constrainedSet) - dab.lbEnlarged(dab.constrainedSet)).^2)...
                +mu*diag(1./(dab.ubEnlarged(dab.constrainedSet) - x(dab.constrainedSet)).^2);
        end
        
        function outmu = updateMu(dab, mu)
            outmu = min(dab.kmu * mu, mu^dab.theta);
        end
    end
    
end
