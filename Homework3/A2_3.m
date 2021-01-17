function A2_3
    tic

    obj_func_exp = Util.get_target_func_exp;
    obj = ObjectiveFunctionR2(obj_func_exp, false);

    tol = 1e-5;
    iter_max = 1e4;

    s = 1;
    sigma = 0.5;
    gamma = 0.1;

    alpha_func = @(k)(1e-2 / (k + 2)^(1/7));

    eposilon = 1e-7;
    m = 25;

    n = 15;

    circle_center = [-5; 0];
    circle_r = 5;

    beta1 = 1e-6;
    beta2 = 1e-6;
    p = .1;

    backtracking = Backtracking(s, sigma, gamma, tol, iter_max);
    dinimishing = Dinimishing(alpha_func, tol, iter_max);
    lm_ada_0 = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max);
    globalized_newton = GlobalizedNewton(beta1, beta2, p, s, sigma, gamma, tol, iter_max);

    solvers = {backtracking, dinimishing, lm_ada_0, globalized_newton};

    x_inits = Util.get_points_on_circle(circle_center, circle_r, n);

    grad_path_tables = cell(1, length(solvers));

    for i = 1:length(grad_path_tables)
        grad_path_tables{i} = cell(1, n);

        for j = 1:length(x_inits)
            x_init = x_inits(:, j);
            [~, ~, ~, ~, ~, ~, steps] = solvers{i}.descend(obj, x_init, true);

            grad_path_tables{i}{j} = Util.sample_seq(steps);

        end

    end

    for i = 1:length(grad_path_tables)
        paths = grad_path_tables{i};
        % use @(x,y)(obj.f([x;y])) will lead to the warning of 'please vectorized your function'
        plot_range = [-15, 5, -6, 6];
        level_list = [1e-5, 1e-3, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7];
        Util.draw_contour_path_image(@obj.f2d, paths, i, plot_range, level_list);
    end

    toc
end
