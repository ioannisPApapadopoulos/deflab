classdef MatlabBackslash
    
    methods
        
        function update = solve(self, jacobian, residual)
            update = -jacobian\residual;
        end
        
    end
end
