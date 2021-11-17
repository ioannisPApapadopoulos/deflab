% Example using ShiftedDeflation vs. ExponentialDeflation
%
% We are solving (x-2)^4 (x-2) = 0 for x. This has the solutions
% x = 1 with multiplicity 4 and x = 2 with multiplicty 1. 
%
% The ShiftedDeflation operator will only decrease the multplicity of the 
% root 1 and the nonlinear solver still converges to 1 after it has been
% deflated away. Whereas the ExponentialDeflation operator removes the root
% 1 and the next solve converges to the other root 2.

addpath(genpath('../deflation'));

% Intialise starting value
x = [0];
dim = length(x);

% Initialise ShiftedDeflation operator, currently holding no known solutions. 
% Since this is a finite-dimensional problem, the inner-product is defined 
% by the identity matrix
deflation = DeflationOperator({}, speye(dim));

% Initialise the solver, and pass the deflation operator into it
nls = NonlinearSolver('deflation', deflation, 'tol', 1e-20);

% Use Newton's method to find the first solution which is 1.
root = nls.newton(x, @F, @jacobian);
disp(root')

% Update found solution list in deflation operator and solve again from
% the same initial guess. We see that as the multiplicity is >1, we
% converge to the same root. 
nls.deflation.updateFoundSolutions(root);
root = nls.newton(x, @F, @jacobian);
disp(root')

% Now initialize a ExponentialDeflation operator and pass into the
% nonlinear solver.
deflation = DeflationOperator({}, speye(dim), 'DeflationOperator', 'ExponentialDeflation');
nls = NonlinearSolver('deflation', deflation, 'tol', 1e-20);

% Again use Newton's method to find the first solution which is 1.
root = nls.newton(x, @F, @jacobian);
disp(root')

% Update found solution list in deflation operator and solve again from
% the same initial guess. The ExponentialOperator completely removes the
% root 1 and we converge to the second root 2. 
nls.deflation.updateFoundSolutions(root);
root = nls.newton(x, @F, @jacobian);
disp(root')

function out = F(x)
    % function handle nonlinear system
    out = [(x-1)^4*(x-2)];
end
        
function out = jacobian(x)
    % Jacobian of nonlinear system
    out = 4*(x-1)^3*(x-2) + (x-1)^4;
end