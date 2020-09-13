% Example using Newton's method & deflation
%
%
% This is a nonlinear system representing the intersection between a straight
% line and a parabola. It consists of two unknowns and has the solutions 
% (1,2) and (0,1). 
%
% x - y = -1
% y = x^2 + 1

addpath(genpath('deflatedbarrier'));

% Intialise starting value
x = [3,3]';
dim = length(x);

% Initialise deflation operator, currently holding no known solutions. 
% Since this is a finite-dimensional problem, the inner-product is defined 
% by the identity matrix
deflation = DeflationOperator({}, speye(dim));

% Initialise the solver, and pass the deflation operator into it
nls = NonlinearSolver('deflation', deflation);

% Use Newton's method to find the first solution
root = nls.newton(x, @F, @jacobian);
disp(root')

% Update found solution list in deflation operator and solve again from
% the same initial guess and using the same solver to find a new solution
nls.deflation.updateFoundSolutions(root);
root = nls.newton(x, @F, @jacobian);
disp(root')

function out = F(xvec)
    % function handle nonlinear system
    x1 = xvec(1);
    x2 = xvec(2);
               
    out =[x1-x2+1;x2-x1^2-1];
end
        
function out = jacobian(xvec)
    x1 = xvec(1);
    % Jacobian of nonlinear system
    out = [1,     -1;...
           -2*x1,  1];
end