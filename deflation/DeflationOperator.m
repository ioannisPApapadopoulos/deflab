classdef DeflationOperator < handle

    properties
        foundRoots
        deflation
        innerProductMatrix
    end

    methods

        function self = DeflationOperator(foundRoots, innerProductMatrix, varargin)

            defaultDeflationOperator = 'ShiftedDeflation';
            p = inputParser;
            addParameter(p, 'DeflationOperator', defaultDeflationOperator);
            parse(p,varargin{:});

            for iter = 1:length(foundRoots)
                if isrow(foundRoots{iter})
                    foundRoots{iter} = foundRoots{iter}';
                end
            end
            self.foundRoots = foundRoots;
            self.innerProductMatrix = innerProductMatrix;
            
            if isequal(p.Results.DeflationOperator, 'ShiftedDeflation')
                self.deflation = ShiftedDeflation(2,1,foundRoots,self.innerProductMatrix,[]);
            elseif isequal(p.Results.DeflationOperator, 'ExponentialDeflation')
                self.deflation = ExponentialDeflation(1,0,foundRoots,self.innerProductMatrix,[]);
            else
                error('DeflationOperator should be of type ShiftedDeflation or ExponentialDeflation')
            end

        end

        function updateFoundSolutions(self, newroot)
            if isrow(newroot) && length(newroot) > 1
                warning('Discovered roots should be given as a column vector')
                newroot = newroot';
            end
            self.foundRoots{end+1} = newroot;
            self.deflation.roots{end+1} = newroot;
        end

        function stepadjustment = deflationStepAdjustment(self, state, update)
            if isempty(self.foundRoots)
                stepadjustment = 1;
            else
                dMy = self.getdMy(state, update);
                minv = 1.0 / self.deflation.evaluate(state);
                stepadjustment = 1/(1 - minv*dMy);
            end
        end

        function out = getdMy(self, state, update)
            deriv = self.deflation.derivative(state);
            % defcon has a minus sign here, but that's because PETSc
            % calculates the update so that state = state - update rather
            % than state = state + update
            out = deriv*self.innerProductMatrix*update;
        end

    end
end
