classdef NonlinearSolver < handle
    
    properties
        LinearSolver
        LineSearch
        max_iter = 100;
        deflation = []
        damping = 1;
        tol = 1e-9;
    end
    
    methods
        
        function [state, residual] = checkArguments(~, state, residual)
            chk = size(state);
            
            if chk(1) > 1 && chk(2) > 1
                error('Initial guess should be a column vector not a matrix')
            elseif chk(2) > 1
                warning('Primal variables should be represented by a column vector not a row vector')
                state = state';
            end
            
            chk = size(residual(state));
            if chk(1) > 1 && chk(2) > 1
                error('System function handle should output a column vector not a matrix')
            elseif chk(2) > 1
                warning('System function handle should output a column vector not a row vector')
                residual = @(x) residual(x)';
            end
                
            
        end
        
        
        function y = linesearchBasic(~, x, update, damping)
            y = x + damping * update;
        end
             
        function nls = NonlinearSolver(varargin)
            
            defaultLinearSolver = 'backslash';
            defaultLinesearch = 'basic';
            defaultDamping = 1;
            defaultMaxIter = 100;
            defaultDeflation = [];
            defaultTol = 1e-9;
            
            
            expectedLinesearch = {'basic'};
            validDeflation = @(x) isequal(class(x), 'DeflationOperator');
            validpos = @(x) isnumeric(x) && x > 0;
            
            p = inputParser;
            addParameter(p, 'linearsolver', defaultLinearSolver, @isstring);
            addParameter(p, 'deflation', defaultDeflation, validDeflation);
            addParameter(p, 'linesearch', defaultLinesearch,...
                @(x) any(validatestring(x,expectedLinesearch)));
            addParameter(p, 'damping', defaultDamping, validpos);
            addParameter(p, 'max_iter', defaultMaxIter, validpos);
            addParameter(p, 'tol', defaultTol, validpos);

            parse(p,varargin{:});
            

            if isequal(p.Results.linearsolver, 'backslash')
                nls.LinearSolver = MatlabBackslash;
            end
            
            if isequal(p.Results.linesearch, 'basic')
                nls.LineSearch = BasicLinesearch;
            end
            
            nls.deflation = p.Results.deflation;
            nls.damping = p.Results.damping;
            nls.max_iter = floor(p.Results.max_iter);
            nls.tol = p.Results.tol;
        end
        
        function root = newton(nls, state, residual, jacobian)
            [state, residual] = nls.checkArguments(state, residual);
            iter = 0;
            x = state;
            evaluatedResidual = residual(x);
            evaluatedJacobian = jacobian(x);
            normEvaluatedResidual = norm(evaluatedResidual,2);
            fprintf('Iteration %i, residual norm = %e\n', iter, normEvaluatedResidual);
            while normEvaluatedResidual > nls.tol  && iter < nls.max_iter
                update = nls.LinearSolver.solve(evaluatedJacobian, evaluatedResidual);
               
                
                
                if ~isempty(nls.deflation)
                    stepadjustment = nls.deflation.deflationStepAdjustment(x, update);
                    update = update * stepadjustment;
                end

                update = nls.LineSearch.adjust(x, update, nls.damping);
                x = x + update;
                evaluatedResidual = residual(x);
                evaluatedJacobian = jacobian(x);
                normEvaluatedResidual = norm(evaluatedResidual,2);
                iter = iter+1;
                fprintf('Iteration %i, residual norm = %e\n',iter, normEvaluatedResidual);
            end
            if iter == nls.max_iter
                disp("Iteration max reached")
            end
            root = x;
        end
        
        function root = bensonmunson(nls, state, lb, ub, residual, jacobian)
            
            [state, residual] = nls.checkArguments(state, residual);
            index = 1:size(state,1);
            iter = 0;
            inactive = index;
            
            x = state;
            x = nls.project(x,lb,ub);
            evaluatedJacobian = jacobian(x);
            evaluatedResidual = residual(x);
            
            normResidualOmega = norm(nls.reducedResidual(evaluatedResidual, x, lb, ub));
            fprintf('Iteration %i, residual norm = %e\n',iter, normResidualOmega);
            
            n = length(x);
            while (normResidualOmega) > 1e-7 && (iter < nls.max_iter)
                
                update = zeros(n, 1);
                update(inactive) = nls.LinearSolver.solve(evaluatedJacobian(inactive, inactive),...
                    evaluatedResidual(inactive));
                update = nls.project(update,lb-x,ub-x);
                
                if norm(update) < 1e-10
                    update = -evaluatedResidual;
                end
                
                if ~isempty(nls.deflation)
                    stepadjustment = nls.deflation.deflationStepAdjustment(x, update);
                    update = update * stepadjustment;
                end
                
                update = nls.LineSearch.adjust(x, update, nls.damping);
                x = x + update;
                x = nls.project(x,lb,ub);
                
                evaluatedResidual = residual(x);
                evaluatedJacobian = jacobian(x);
                normResidualOmega = norm(nls.reducedResidual(evaluatedResidual, x, lb, ub));
                
                active_lb = index(x<=lb);
                active_lb = intersect(active_lb,index(evaluatedResidual>0));
                active_ub = index(x>=ub);
                active_ub = intersect(active_ub,index(evaluatedResidual<0));
                active = [active_lb, active_ub];
                index2 = index;
                index2(active) = 0;
                inactive  = find(index2);
                
                iter = iter + 1;
                fprintf('Iteration %i, residual norm = %e\n',iter, normResidualOmega);
            end
            
            if iter == nls.max_iter
                disp("Iteration max reached")
                disp(normResidualOmega)
            end
            root = x;
        end
        
        function out = reducedResidual(nls, evaluatedResidual, x, lb, ub)
            out = evaluatedResidual;
            out(x<=lb) = min(out(x<=lb),0);
            out(x>=ub) = max(out(x>=ub),0);
        end
        
        function b = project(nls, x, lb, ub)
            b = x;
            b(x<lb) = lb(x<lb);
            b(x>ub) = ub(x>ub);
        end
        
        function root = hik(nls, state, lb, ub, residual, jacobian)
            [state, residual] = nls.checkArguments(state, residual);
            index = 1:size(state,1);
            iter = 0;
            inactive = index;
            active_lb = [];
            active_ub = [];
            active    = [];
            
            x = state;
            x = nls.project(x,lb,ub);
            n = length(x);
            
            evaluatedJacobian = jacobian(x);
            evaluatedResidual = residual(x);
            dual = zeros(n,1);
            
            normResidualOmega = norm(nls.reducedResidual(evaluatedResidual, x, lb, ub));
            %normResidualOmega = norm(dual - max(zeros(n,1), dual+x-ub) - max(zeros(n,1), dual-x-lb));
            %normResidualOmega = normResidualOmega + norm(evaluatedResidual+dual);
            fprintf('Iteration %i, residual norm = %e\n',iter, normResidualOmega);
            
            
            while (normResidualOmega) > 1e-7 && (iter < nls.max_iter)
                
                update = zeros(n, 1);
                update(active_lb) = lb(active_lb) - x(active_lb);
                update(active_ub) = ub(active_ub) - x(active_ub);
                dual(inactive) = 0;
                if ~isempty(active)
                    correctedResidual = evaluatedResidual(inactive) + evaluatedJacobian(inactive, active)*update(active);
                else
                    correctedResidual = evaluatedResidual(inactive);
                end
                update(inactive) = nls.LinearSolver.solve(evaluatedJacobian(inactive, inactive),...
                    correctedResidual);
                
                if ~isempty(nls.deflation)
                    stepadjustment = nls.deflation.deflationStepAdjustment(x, update);
                    update = update * stepadjustment;
                end
                
                update = nls.LineSearch.adjust(x, update, nls.damping);
                x = x + update;
                
                evaluatedResidual = residual(x);
                evaluatedJacobian = jacobian(x);

                % normResidualOmega = norm(nls.reducedResidual(evaluatedResidual, x, lb, ub));
                
                % which way should the sign be?
                dual(active_lb) = evaluatedResidual(active_lb);
                dual(active_ub) = -evaluatedResidual(active_ub);
                
                active_lb = index((dual-x+lb)>0);
                active_ub = index((dual-ub+x)>0);
                active = [active_lb, active_ub];
                tmp_index = index;
                tmp_index(active) = 0;
                inactive  = find(tmp_index);
                
                %normResidualOmega = norm(dual - max(zeros(n,1), dual+x-ub) - max(zeros(n,1), dual-x+lb));
                %normResidualOmega = normResidualOmega + norm(evaluatedResidual+dual);
                y = nls.project(x,lb,ub);
                evaluatedResidualy = residual(y);
                normResidualOmega = norm(nls.reducedResidual(evaluatedResidualy, y, lb, ub));
                iter = iter + 1;
                fprintf('Iteration %i, residual norm = %e\n',iter, normResidualOmega);
            end
            
            if iter == nls.max_iter
                disp("Iteration max reached")
                disp(normResidualOmega)
            end
            root = x;
        end
        
        function root = ssls(nls, state, lb, ub, residual, jacobian)
            [state, residual] = nls.checkArguments(state, residual);
            iter = 0;
            
            x = state;
            x = nls.project(x,lb,ub);
            
            tol = 1e10;
            wherenoconstraint = intersect(find(lb<-tol), find(ub>tol));
            wherelbconstraint = intersect(find(lb>=-tol), find(ub>tol));
            whereubconstraint = intersect(find(lb<-tol), find(ub<=tol));
            whereequalconstraint = find(lb==ub);
            wherebothconstraint = intersect(find(lb>=-tol), find(ub<=tol));
            
            evaluatedJacobian = jacobian(x);
            evaluatedResidual = residual(x);
                     
            fb = nls.FB(x, evaluatedResidual, lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint);
            normFB = norm(fb);
            fprintf('Iteration %i, residual norm = %e\n',iter, normFB);
            
            
            while (normFB) > 1e-7 && (iter < nls.max_iter)
                
                [dshift, dscale] = nls.computeScaleAndShift(x, evaluatedResidual,lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint);
                shiftedJacobian = diag(dshift) + diag(dscale) * evaluatedJacobian;

                update = nls.LinearSolver.solve(shiftedJacobian,fb);
                
                if ~isempty(nls.deflation)
                    stepadjustment = nls.deflation.deflationStepAdjustment(x, update);
                    update = update * stepadjustment;
                end
                
                update = nls.LineSearch.adjust(x, update, nls.damping);
                x = x + update;
                
                evaluatedResidual = residual(x);
                evaluatedJacobian = jacobian(x);                
              
                fb = nls.FB(x, evaluatedResidual, lb, ub, wherenoconstraint,wherelbconstraint,whereubconstraint,whereequalconstraint,wherebothconstraint);
                normFB = norm(fb);                
                iter = iter + 1;
                fprintf('Iteration %i, residual norm = %e\n',iter, normFB);
            end
            
            if iter == nls.max_iter
                disp("Iteration max reached")
            end
            root = x;
        end
        
        function out = Phi(nls, a, b)
            out = a + b - sqrt(a.^2 + b.^2);      
        end
        
        function out = dPhi(nls, a, b)
            if any(abs(a) > 1e-6) || any(abs(b) > 1e-6)
                out =  1 - a./sqrt(a.^2 + b.^2);
            else
                out = 0.5;
            end
        end
        
        function out = FB(nls, state, residual, lb, ub, wherenoconstraint, wherelbconstraint, whereubconstraint, whereequalconstraint, wherebothconstraint)
            out = zeros(length(state),1);
            % FIXME add a check for all indices here
            
            out(wherenoconstraint) = residual(wherenoconstraint);
            
            idx = whereubconstraint;
            out(idx) = nls.Phi(ub(idx) - state(idx), -residual(idx));
            
            idx = wherelbconstraint;
            out(idx) = nls.Phi(state(idx)-lb(idx),residual(idx));
            
            idx = wherebothconstraint;
            out(idx) = nls.Phi(state(idx)-lb(idx),-nls.Phi(ub(idx)-state(idx),-residual(idx)));
            
            idx = whereequalconstraint;
            out(idx) = lb(idx) - state(idx);
        end

        function [dshift, dscale] = computeScaleAndShift(nls, x, residual, lb, ub, wherenoconstraint, wherelbconstraint, whereubconstraint, whereequalconstraint, wherebothconstraint)
            n = length(x);
            dshift = ones(n, 1);
            dscale = ones(n, 1);
            
            dshift(wherenoconstraint) = 0;
            dscale(wherenoconstraint) = 1;
            
            idx = whereubconstraint;
            dshift(idx) = nls.dPhi(ub(idx) - x(idx), -residual(idx));
            dscale(idx) = nls.dPhi(-residual(idx),ub(idx)-x(idx));
            
            idx = wherelbconstraint;
            dshift(idx) = nls.dPhi(x(idx) - lb(idx), residual(idx));
            dscale(idx) = nls.dPhi(residual(idx), x(idx) - lb(idx));
           
            
            idx = wherebothconstraint;
            dshift1 = nls.dPhi(x(idx) - lb(idx), -nls.Phi(ub(idx) - x(idx), -residual(idx)));
            dscale1 = nls.dPhi(-nls.Phi(ub(idx)-x(idx),-residual(idx)), x(idx) - lb(idx));
            dshift2 = nls.dPhi(ub(idx)-x(idx),-residual(idx));
            dscale2 = nls.dPhi(-residual(idx),ub(idx) - x(idx));
            dshift(idx) = dshift1 + dscale1.*dshift2;
            dscale(idx) = dscale1.*dscale2;
            
            idx = whereequalconstraint;
            dshift(idx) = 1;
            dscale(idx) = 0;
   
        end
    end
    
end
