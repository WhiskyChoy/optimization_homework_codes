classdef Util

    methods (Static)

        function check_zero_to_one(param)

            if (param < 0 || param > 1)
                eid = 'Param:wrongplot_Range';
                msg = sprintf('The parameter %.f is outside of the plot_range (0, 1)', param);
                throwAsCaller(MException(eid, msg))
            end

        end

        function test_all_x_inits(obj_func, gradient_method, x_inits)
            fprintf(gradient_method.self_info);

            for x_init = x_inits
                [x_p, value_p, norm_p, iter, t, success] = gradient_method.descend(obj_func, x_init);
                Util.print_descend_msg(x_init(1), x_init(2), x_p(1), x_p(2), value_p, norm_p, iter, t, success);
            end

            fprintf('\n')
        end

        function print_descend_msg(x1_init, x2_init, x1_p, x2_p, value_p, norm_p, iter, t, success)
            fprintf('| Method: Backtracking | Initial point: [%f, %f]^T | Solution: [%s, %s]^T | Value: %s | Norm: %s | Iter: %d | Time: %.2f second(s) | Success: %d |\n', x1_init, x2_init, Util.smart_num2str(x1_p), Util.smart_num2str(x2_p), Util.smart_num2str(value_p), Util.smart_num2str(norm_p), iter, t, success);
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

        function seq = sample_seq(seq)
            save_num = 100;
            sample_rate = 0.1;
            seq_length = length(seq);

            if (seq_length > save_num)
                sep = ceil(sample_rate * (seq_length - save_num));
                seq = [seq(:, 1:save_num), seq(:, save_num + 1:sep:end)];
            end

        end

        function draw_contour_path_image(obj_func, paths, graph_index, plot_range, level_list)

            arguments
                obj_func function_handle
                paths cell
                graph_index
                plot_range
                level_list
            end

            figure(graph_index);
            % plot_range: xmin xmax ymin ymax
            % fc = fcontour(obj_func, plot_range);
            % If we limit the plot_range of countour we can't see other part of the function with contour
            fc = fcontour(obj_func);

            if ~isempty(level_list)
                fc.LevelList = level_list;
            end

            colorbar;
            hold on;

            for i = 1:length(paths)
                path = paths{i};
                % disp(path);
                plot(path(1, 1), path(2, 1), 'x');
                plot(path(1, :), path(2, :));
                plot(path(1, end), path(2, end), 'o');
                hold on;
            end

            axis equal;
            if ~isempty(plot_range)
                axis(plot_range);
            end
        end

        function draw_matrix_rows(A, graph_index)
            figure(graph_index);

            num_rows = size(A, 1);

            for i = 1:num_rows
                plot(A(i,1),A(i,2), '+');
                hold on;
            end

        end

        function draw_one_dim(x, graph_index)
            figure(graph_index);
            plot(x);
            hold on;
        end

        function draw_one_dim_log(x, graph_index, color)
            figure(graph_index);
            loglog(x, color);
            hold on;
        end
    end

end
