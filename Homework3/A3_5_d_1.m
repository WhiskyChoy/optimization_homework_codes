function A3_5_d_1
    tic

    beta1 = 1e-6;
    beta2 = 1e-6;
    p = .1;
    s = 1;
    sigma = .5;
    gamma = .1;
    tol = 1e-6;
    iter_max = 1e4;

    gn_solver = GlobalizedNewton(beta1, beta2, p, s, sigma, gamma, tol, iter_max);
    disp(gn_solver.self_info());

    m_s = [2, 3, 5];
    n = 2;

    num_m_s = size(m_s, 2);

    obj_s = cell(num_m_s);

    for i = 1:num_m_s
        A = randi(10+randi(90), m_s(i), n);
        obj = ObjectiveFunctionAx(A);
        obj_s{i} = obj;
    end

    for i = 1:num_m_s
        obj = obj_s{i};
        x_init = zeros(2, 1);
        [~, ~, ~, ~, ~, ~, steps] = gn_solver.descend(obj, x_init, true);
        paths = {steps};
        plot_range = [-1.5, .5, -1.5, .5];
        disp(obj.A);
        Util.draw_contour_path_image(@obj.f2d, paths, i, plot_range, []);
    end

    toc
end
