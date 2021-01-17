function A2_1n2
    %-------------------------------Common Part------------------------------%
    syms x1 x2
    f1 = -1 + x1 + ((5 - x2) * x2 - 2) * x2;
    f2 = -1 + x1 + ((x2 + 1) * x2 - 10) * x2;
    obj_func_exp = f1^2 +f2^2;

    obj = ObjectiveFunctionR2(obj_func_exp, true);

    tol = 1e-5;
    iter_max = 1e4;

    s = 1;
    sigma = 0.5;
    gamma = 0.1;

    x_inits = [[-5; 0], [0; 0], [5; 0]];

    %-------------------------- Solution to A2.3(a) --------------------------%
    backtracking = Backtracking(s, sigma, gamma, tol, iter_max);

    alpha_func_1 = @(k)(1e-2 / (k + 2)^(1/7));
    alpha_func_2 = @(k)(1e-2 / log(k + 2));

    dinimishing_1 = Dinimishing(alpha_func_1, tol, iter_max);
    dinimishing_2 = Dinimishing(alpha_func_2, tol, iter_max);

    Util.test_all_x_inits(obj, backtracking, x_inits);
    Util.test_all_x_inits(obj, dinimishing_1, x_inits);
    Util.test_all_x_inits(obj, dinimishing_2, x_inits);

    %-------------------------- Solution to A2.3(b) --------------------------%

    eposilon = 1e-7;
    m = 25;

    lm_ada_0 = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max);
    Util.test_all_x_inits(obj, lm_ada_0, x_inits);

    ep_seq = [1e-5, 1e-6, 1e-7, 1e-8];
    m_seq = [5, 10, 15, 25, 50, 100, 200];

    [iter_count_tables, t_count_tables] = deal(cell(1, length(x_inits)));

    for i = 1:length(x_inits)
        [iter_count_tables{i}, t_count_tables{i}] = deal(zeros(length(ep_seq), length(m_seq)));
        x_init = x_inits(:,i);

        for j = 1:length(ep_seq)

            for k = 1:length(m_seq)
                lm_ada = LMAda(ep_seq(j), m_seq(k), s, sigma, gamma, tol, iter_max);
                [~, ~, ~, iter, t, success] = lm_ada.descend(obj, x_init);

                if (~success)
                    [iter, t] = deal(-1);
                end

                iter_count_tables{i}(j, k) = iter;
                t_count_tables{i}(j, k) = t;
            end

        end

    end

    fprintf('For the rows we chose eposilon from:\n\n');
    disp(ep_seq);
    fprintf('For the columns we choose memory size from:\n\n');
    disp(m_seq);

    for i = 1:length(iter_count_tables)
        x_init = x_inits(:,i);
        fprintf('With initial point [%f, %f]^T, we have the table which takes the count of [iteration]  listed below:\n\n', x_init(1), x_init(2));
        disp(iter_count_tables{i})
    end

    for i = 1:length(t_count_tables)
        x_init = x_inits(:,i);
        fprintf('With initial point [%f, %f]^T, we have the table which takes the count of [time] table listed below:\n\n', x_init(1), x_init(2));
        disp(t_count_tables{i})
    end

end
