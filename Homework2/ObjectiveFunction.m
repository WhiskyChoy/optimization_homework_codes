classdef ObjectiveFunction

    properties (Access = private)
        obj
    end

    properties (Dependent)
        solutions
    end

    methods

        function self = ObjectiveFunction(obj_func_exp, show_msg)

            arguments
                obj_func_exp (1, 1) sym {quadratic_check} = get_default_function_exp()
                show_msg (1, 1) logical = true
            end

            self.obj = get_objective_function(obj_func_exp, show_msg);
        end

        function res = f(self, x1, x2)
            res = self.obj.f(x1, x2);
        end

        function res = g(self, x1, x2)
            res = self.obj.g(x1, x2);
        end

        function res = h(self, x1, x2)
            res = self.obj.h(x1, x2);
        end

        function res = get.solutions(self)
            res = self.obj.solutions;
        end

    end

end

function obj = get_objective_function(obj_func_exp, show_msg)

    g_exp = gradient(obj_func_exp);
    % We can also use "d_x1_exp = diff(obj_func_exp, x1);"
    %             and "d_x2_exp = diff(obj_func_exp, x2);"
    %             to calculate the partial derivative!
    d_x1_exp = g_exp(1);
    d_x2_exp = g_exp(2);
    % we can use the "hess" function to calculate hessian
    % but it's a waste of computational resource since we already have the gradient
    h_exp = jacobian(g_exp);

    obj.f = get_matlab_func_with_original_param(obj_func_exp);
    obj.g = get_matlab_func_with_original_param(g_exp);
    obj.h = get_matlab_func_with_original_param(h_exp);

    % Solve should be use after using matlabFunction in vscode.
    % It's quite strange that we don't need to do so in matlab.
    sol_struct = solve(d_x1_exp, d_x2_exp);
    sol_length = length(sol_struct.x1);
    solutions = cell(1, sol_length);

    % point_generate = @(x1, x2, value, float_x1, float_x2, float_value, s_type) ...
    %     (struct('x1', x1, 'x2', x2, 'value', value, ...
    %     'float_x1', float_x1, 'float_x2', float_x2, ...
    %     'float_value', float_value, 'type', s_type));

    get_precise_value = @(x1, x2)(subs(obj_func_exp, {'x1', 'x2'}, {x1, x2}));

    [max_func_value, min_func_value] = deal(obj.f(double(sol_struct.x1(1)), double(sol_struct.x2(1))));

    for i = 1:sol_length
        item_x1 = sol_struct.x1(i);
        item_x2 = sol_struct.x2(i);
        func_value = get_precise_value(item_x1, item_x2);

        float_func_value = double(func_value);
        float_x1 = double(item_x1);
        float_x2 = double(item_x2);

        if (float_func_value > max_func_value)
            max_func_value = float_func_value;
        elseif (float_func_value < min_func_value)
            min_func_value = float_func_value;
        end

        wrap_item = StationaryPoint(item_x1, item_x2, func_value, float_x1, float_x2, float_func_value, StationaryType.Unknown);
        solutions{i} = wrap_item;
    end

    for i = 1:sol_length
        % This will chnage the result. If we use item = solutions{i} and item = ... ,
        % the item is temporal and will not change the result of solutions{i}
        solutions{i}.type = stationary_check(obj, solutions{i}, max_func_value, min_func_value);
    end

    obj.solutions = solutions;

    if (show_msg)
        fprintf('The objective function is: %s\n', obj_func_exp);
        fprintf('The partial derivative of x1 for the objective function is: %s\n', d_x1_exp);
        fprintf('The partial derivative of x2 for the objective function is: %s\n', d_x2_exp);

        fprintf('\n')

        for i = 1:length(solutions)
            solution = solutions{i};
            fprintf('(Symbolic)  The stationary point with index %d is: (%s, %s) with the function value %s.\n', i, solution.x1, solution.x2, solution.value);
            fprintf('(Numerical) The stationary point with index %d is: (%.2f, %.2f) with the function value %.2f.\n', i, solution.float_x1, solution.float_x2, solution.float_value')
            fprintf('The type description for this stationary point is "%s".\n\n', solution.type.value);
        end

    end

end

function result_func = get_matlab_func_with_original_param(exp)
    syms x1 x2
    symvars = symvar(exp);
    matlab_func = matlabFunction(exp);

    if (isempty(symvars))
        result_func = @(~, ~)(matlab_func());
    elseif (length(symvars) == 1)
        sym_1 = symvars(1);

        if (sym_1 == x1)
            result_func = @(input, ~)(matlab_func(input));
        elseif (sym_1 == x2)
            result_func = @(~, input)(matlab_func(input));
        else
            eid = 'Symvars:wrongArgs';
            msg = sprintf('Symbolic variable "%s" encountered, but we expected "%s" or "%s"', sym_1, x1, x2);
            throwAsCaller(MException(eid, msg));
        end

    elseif (length(symvars) == 2)
        sym_1 = symvars(1);
        sym_2 = symvars(2);

        if (sym_1 == x1 && sym_2 == x2)
            result_func = matlab_func;
        else
            eid = 'Symvars:wrongArgs';
            msg = sprintf('Symbolic variable "%s" and "%s" encountered, but we expected "%s" and "%s"', sym_1, sym_2, x1, x2);
            throwAsCaller(MException(eid, msg));
        end

    else
        eid = 'Symvars:wrongNum';
        msg = 'The number of symvars in the function expression should be 0, 1 or 2 here.';
        throwAsCaller(MException(eid, msg));
    end

end

function quadratic_check(obj_func_exp)
    symvars = symvar(obj_func_exp);
    syms x1 x2

    if (length(symvars) ~= 2)
        eid = 'Symvars:wrongNum';
        msg = 'The number of symvars in the function expression is not 2, which indicates that it''s not a quadratic function.';
        throwAsCaller(MException(eid, msg));
    elseif (symvars(1) ~= x1 || symvars(2) ~= x2)
        eid = 'Symvars:wrongArgs';
        msg = sprintf('The symbolic variable you use is "%s" and "%s". Please use "x1" and "x2" as the symbolic variables.', symvars(1), symvars(2));
        throwAsCaller(MException(eid, msg));
    end

end

function obj_func_exp = get_default_function_exp
    syms x1 x2
    f1 = -1 + x1 + ((5 - x2) * x2 - 2) * x2;
    f2 = -1 + x1 + ((x2 + 1) * x2 - 10) * x2;
    obj_func_exp = f1^2 +f2^2;
end

function check_result = stationary_check(obj, point, max_val, min_val)
    x1 = point.float_x1;
    x2 = point.float_x2;
    value = point.float_value;
    h = obj.h(x1, x2);
    h_det = det(h);
    h_a = h(1, 1);
    almost_equal = @(a, b, tol)(abs(a - b) < tol);

    if (h_det > 0)

        if (h_a > 0)

            if (almost_equal(value, min_val, 1e-5))
                check_result = StationaryType.GlobalMin;
            else
                check_result = StationaryType.LocalMin;
            end

        elseif (h_a < 0)

            if (almost_equal(value, max_val, 1e-5))
                check_result = StationaryType.GlobalMax;
            else
                check_result = StationaryType.LocalMax;
            end

        else
            check_result = StationaryType.Error;
        end

    elseif (h_det < 0)
        check_result = StationaryType.SaddlePoint;
    else
        check_result = StationaryType.NeedAly;
    end

end
