classdef Util

    methods (Static)

        function check_zero_to_one(param)

            if (param < 0 || param > 1)
                eid = 'Param:wrongRange';
                msg = sprintf('The parameter %.f is outside of the range (0, 1)', param);
                throwAsCaller(MException(eid, msg))
            end

        end

        function test_all_x_inits(obj_func, gradient_method, x_inits)
            fprintf(gradient_method.self_info);

            for x_init = x_inits
                [x1_p, x2_p, value_p, iter, t, success] = gradient_method.descend(obj_func, x_init(1), x_init(2));
                Util.print_descend_msg(x_init(1), x_init(2), x1_p, x2_p, value_p, iter, t, success);
            end

            fprintf('\n')
        end

        function print_descend_msg(x1_init, x2_init, x1_p, x2_p, value_p, iter, t, success)
            fprintf('| Method: Backtracking | Initial point: [%f, %f]^T | Solution: [%s, %s]^T | Value: %s | Iter: %d | Time: %.2f second(s) | Success: %d |\n', x1_init, x2_init, Util.smart_num2str(x1_p), Util.smart_num2str(x2_p), Util.smart_num2str(value_p), iter, t, success);
        end

        function str = smart_num2str(num)

            if (isinf(num))
                str = num2str(num);
            elseif (abs(num) > 1e8)
                str = num2str(num, '%.2e');
            else
                str = num2str(num, '%.2f');
            end

        end

        function obj_func_exp = get_target_func_exp
            syms x1 x2
            f1 = -1 + x1 + ((5 - x2) * x2 - 2) * x2;
            f2 = -1 + x1 + ((x2 + 1) * x2 - 10) * x2;
            obj_func_exp = f1^2 +f2^2;
        end

        function points = get_points_on_circle(center, r, num)
            seq = 1:num;
            deg = 1 / num * 2 * pi * (seq - 1);
            points = center + r * [cos(deg); sin(deg)];
        end

    end

end
