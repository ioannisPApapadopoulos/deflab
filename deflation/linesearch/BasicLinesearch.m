classdef BasicLinesearch
    
    methods
        
        function update = adjust(~, ~, update, damping)
            update = update * damping;
        end
        
    end
end
