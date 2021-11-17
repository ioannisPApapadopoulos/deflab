% Example using a semismooth Newton method & deflation
%
%
% The following description is taken from Deflation for semismooth
% equations by Farrell, Croci and Surowiec (2019) doi: 10.1080/10556788.2019.1613655 
%
% This is a nonconvex quadratic programming problem with linear constraints 
% suggestedby N. I. M. Gould in an invited lecture to the 19th biennial 
% conference on NA. It is a quadratic minimization problem with an indefinite 
% Hessian of the form min_x f(x) = −2(x1−1/4)^2 + 2(x2−1/2)^2
% s.t. x1+x2<=1, 6*x1+2*x2<=3, x1,x2>=0.

% The first order Karush–Kuhn–Tucker optimality conditions yield an NCP.
% The nonconvexity of the function makes this problem difficult; 
% it attains two minima with similar functional values and has a saddle point 
% at x = [1/4,1/2]^T. The central path to be followed by an interior point method 
% is pathological, with different paths converging to the different minima.

addpath(genpath('../deflation'));

% Intialise starting value
x = [0.2,0.2,0,0]';

% Lower and upper bounds
lb = [0,0,0,0]';
ub = [1e100,1e100,1e100,1e100]';
dim = length(x);

% Initialise deflation operator, currently holding no known solutions. 
% Since this is a finite-dimensional problem, the inner-product is defined 
% by the identity matrix
deflation = DeflationOperator({}, speye(dim));

% Initialise the solver, and pass the deflation operator into it
nls = NonlinearSolver('deflation', deflation);

% Use Benson and Munson SSLS solver to find the first solution
root = nls.ssls(x, lb, ub, @F, @jacobian);
disp(root')

% Update found solution list in deflation operator and solve again from
% the same initial guess and using the same solver to find a new solution
nls.deflation.updateFoundSolutions(root);
root = nls.ssls(x, lb, ub, @F, @jacobian);
disp(root')

% Repeat again to find the final solution
nls.deflation.updateFoundSolutions(root);
root = nls.ssls(x, lb, ub, @F, @jacobian);
disp(root')


function out = F(xvec)
    % function handle for the KKT conditions
    x1 = xvec(1);
    x2 = xvec(2);
    l1 = xvec(3);
    l2 = xvec(4);
               
    out =[-4*(x1-1/4)+3*l1+l2;4*(x2-1/2)+l1+l2;3-6*x1-2*x2;1-x1-x2];
end
        
function out = jacobian(~)
    % Jacobian of the KKT conditions
    out = [-4, 0, 3, 1;...
            0, 4, 1, 1;...
           -6, -2 ,0, 0;...
           -1, -1, 0, 0];
end
