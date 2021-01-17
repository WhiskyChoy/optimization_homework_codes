classdef ObjectiveFunctionAx < ObjectiveFunction
    properties (Access = private)
        obj
    end

    properties (Dependent = true)
        A
    end

    methods
        function self = ObjectiveFunctionAx(A)
            self.obj = get_objective_function_ax(A);
        end

        function res = f2d(self, x1, x2)
            res = self.obj.f(x1, x2);
        end

        function result = f(self, x)
            params = num2cell(x);
            result = self.obj.f(params{:});
        end

        function result = g(self, x)
            params = num2cell(x);
            result = self.obj.g(params{:});
        end

        function result = h(self, x)
            params = num2cell(x);
            result = self.obj.h(params{:});
        end

        function result = get.A(self)
            result = self.obj.A;
        end
    end
end

function obj = get_objective_function_ax(A)
    n = size(A, 2);
    x = sym('x', [n 1]);
    % x.^2, not x^.2
    % sum(A) or sum(A, 1) will add the sum in columns
    obj_func_exp = -sum(log(1 - A * x))-sum(log(1-x.^2));
    g_exp = gradient(obj_func_exp);
    h_exp = jacobian(g_exp);

    obj.f = matlabFunction(obj_func_exp);
    obj.g = matlabFunction(g_exp);
    obj.h = matlabFunction(h_exp);
    obj.A = A;
end