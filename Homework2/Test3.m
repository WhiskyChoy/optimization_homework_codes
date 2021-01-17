function Test3
    tic

    obj_func_exp = Util.get_target_func_exp;
    obj = ObjectiveFunction(obj_func_exp, false);

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

    backtracking = Backtracking(s, sigma, gamma, tol, iter_max);
    dinimishing = Dinimishing(alpha_func, tol, iter_max);
    lm_ada_0 = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max);

    gradient_methods = {backtracking, dinimishing, lm_ada_0};

    x_inits = Util.get_points_on_circle(circle_center, circle_r, n);

    grad_path_tables = cell(1, length(gradient_methods));

    for i = 1:length(grad_path_tables)
        grad_path_tables{i} = cell(1, n);

        for j = 1:length(x_inits)
            x_init = x_inits(:, j);
            [~, ~, ~, ~, ~, ~, steps] = gradient_methods{i}.descend(obj, x_init(1), x_init(2), true);

            grad_path_tables{i}{j} = sample_seq(steps);

        end

    end

    for i = 1:length(grad_path_tables)
        paths = grad_path_tables{i};
        draw_image(@obj.f, paths, i);
    end

    toc
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

function draw_image(obj_func, paths, graph_index)
    figure(graph_index);
    % range: xmin xmax ymin ymax
    range = [-15, 5, -6, 6];
    % fc = fcontour(obj_func, range);
    % If we limit the range of countour we can't see other part of the function with contour
    fc = fcontour(obj_func);
    fc.LevelList = [1e-5, 1e-3, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7];
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
    axis(range);
end
