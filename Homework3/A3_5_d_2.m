function A3_5_d_2
    tic

    m = 5;
    n = 2;
    A = -.5 + rand(m, n);
    fprintf('We use the matrix A listed below for experiment:\n\n')
    disp(A)
    obj = ObjectiveFunctionAx(A);

    tol = 1e-5;
    iter_max = 1e4;

    s = 1;
    sigma = 0.5;
    gamma = 0.1;

    alpha_func = @(k)(1e-2 / (k + 2)^(1/7));

    eposilon = 1e-7;
    m = 25;

    beta1 = 1e-6;
    beta2 = 1e-6;
    p = .1;

    backtracking = Backtracking(s, sigma, gamma, tol, iter_max);
    dinimishing = Dinimishing(alpha_func, tol, iter_max);
    lm_ada_0 = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max);
    globalized_newton = GlobalizedNewton(beta1, beta2, p, s, sigma, gamma, tol, iter_max);

    solvers = {backtracking, dinimishing, lm_ada_0, globalized_newton};

    x_init = zeros(n, 1);

    grad_path_tables = cell(1, length(solvers));

    for i = 1:length(grad_path_tables)
        [~, ~, ~, ~, ~, ~, steps] = solvers{i}.descend(obj, x_init, true);
        grad_path_tables{i} = Util.sample_seq(steps);
    end

    for i = 1:length(grad_path_tables)
        % {A{i}} -> A(i)
        paths = grad_path_tables(i);
        % use @(x,y)(obj.f([x;y])) will lead to the warning of 'please vectorized your function'
        plot_range = [-.5, .5, -.5, .5];
        level_list = [];
        Util.draw_contour_path_image(@obj.f2d, paths, i, plot_range, level_list);
        Util.draw_matrix_rows(A, i);
    end

    toc
end