classdef GradientMethod
    % When you define any abstract methods or properties,
    % MATLAB automatically sets the class Abstract attribute to true.
    % If you set Access = private, the desendant of this class can not visit the property
    properties (Access = protected)
        tol (1, 1) double = 1e-5
        iter_max (1, 1) int32 = 1e4
    end

    properties (Abstract)
        name (1,1) string
    end

    methods

        function self = GradientMethod(tol, iter_max)

            arguments
                tol (1, 1) double = 1e-5
                iter_max (1, 1) int32 = 1e4
            end

            self.tol = tol;
            self.iter_max = iter_max;
        end
    end

    methods (Abstract)
        [x1_p, x2_p, value_p, iter, t, success, steps] = descend(self, obj, init_x1, init_x2, record)
        str = self_info(self)
    end

end
