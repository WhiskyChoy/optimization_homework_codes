function A3_5_1
    beta1 = 1e-6;
    beta2 = 1e-6;
    p = .1;
    s = 1;
    sigma = .5;
    gamma = .1;
    tol = 1e-6;
    iter_max = 1e4;

    alpha_func = @(k)(1e-2 / (k + 2)^(1/7));

    eposilon = 1e-7;
    m = 25;

    backtracking = Backtracking(s, sigma, gamma, tol, iter_max);
    dinimishing = Dinimishing(alpha_func, tol, iter_max);
    lm_ada = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max);

    gradient_methods = {backtracking, dinimishing, lm_ada};

    globalized_newton = GlobalizedNewton(beta1, beta2, p, s, sigma, gamma, tol, iter_max);

    colors = ['r', 'g', 'b', 'm'];

    m_s = [2, 3, 5];
    n_s = [2, 5, 8];

    num_m_s = size(m_s, 2);
    num_n_s = size(n_s, 2);

    get_init_x = @(n)zeros(n, 1);

    obj_s = cell(num_m_s, num_n_s);

    for i = 1:num_m_s

        for j = 1:num_n_s
            A = randi(100, m_s(i), n_s(j));
            obj = ObjectiveFunctionAx(A);
            obj_s{i, j} = obj;
        end

    end

    for i = 1:num_m_s

        for j = 1:num_n_s
            obj = obj_s{i, j};
            fprintf('Use the matrix A(%d X %d) below for experiment:\n\n', m_s(i), n_s(j));
            disp(obj.A)
            init_x = get_init_x(n_s(j));

            [~, value_p, norm_p, iter, t, success, ~, norm_gs, ~, use_newtons] = globalized_newton.descend(obj, init_x, false, true, false);
            show_descent_msg(globalized_newton.name, value_p, norm_p, iter, t, success);

            index = (i - 1) * num_n_s + j;

            Util.draw_one_dim(use_newtons, num_m_s * num_n_s + index);

            color_index = 1;
            Util.draw_one_dim_log(norm_gs, index, colors(color_index));

            for gradient_method_cell = gradient_methods
                color_index = color_index + 1;
                gradient_method = gradient_method_cell{1};
                [~, value_p, norm_p, iter, t, success, ~, norm_gs, ~] = gradient_method.descend(obj, init_x, false, true, false);
                show_descent_msg(gradient_method.name, value_p, norm_p, iter, t, success);
                Util.draw_one_dim_log(norm_gs, index, colors(color_index));
            end

            fprintf('\n\n')
        end

    end

end

function show_descent_msg(method_name, value_p, norm_p, iter, t, success)
    fprintf('| Method: %s | Function Value: %s | Norm: %s | Iter: %d | Time: %.2f second(s) | Success: %d |\n', method_name, Util.smart_num2str(value_p), Util.smart_num2str(norm_p), iter, t, success);
end
