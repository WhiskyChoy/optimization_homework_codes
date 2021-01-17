classdef LMAda < Backtracking

    properties
        % The "name" property is already defined in Backtracking
        % name = 'Limited-Memory Version of Adagrad'
        eposilon (1, 1) double = 1e-7
        m (1, 1) int32 = 25
    end

    methods

        function self = LMAda(eposilon, m, s, sigma, gamma, tol, iter_max)

            arguments
                eposilon (1, 1) double = 1e-7
                m (1, 1) int32 = 25
                s (1, 1) double = 1
                sigma (1, 1) double {Util.check_zero_to_one} = 0.5
                gamma (1, 1) double {Util.check_zero_to_one} = 0.1
                tol (1, 1) double = 1e-5
                iter_max (1, 1) int32 = 1e4
            end

            self = self@Backtracking(s, sigma, gamma, tol, iter_max);
            self.name = 'Limited-Memory Version of Adagrad';
            self.eposilon = eposilon;
            self.m = m;
        end

        function [x_p, value_p, norm_p, iter, t, success, steps, norm_gs, values] = descend(self, obj, init_x, record_steps, record_norm_gs, record_values)

            arguments
                self
                obj (1, 1) ObjectiveFunction
                init_x double
                record_steps = false
                record_norm_gs = false
                record_values = false
            end

            x_rows = size(init_x, 1);
            g_mem = zeros(x_rows, self.m);
            get_D_k = @(g_mem, eposilon)(diag(sqrt(sum(g_mem.^2, 2) + eposilon).^ - 1, 0));

            tic

            iter = 0;
            x = init_x;
            size_x = size(x, 1);
            g = obj.g(x);
            norm_g = norm(g, 2);

            if (record_steps)
                steps = zeros(size_x, self.iter_max);
            end

            if (record_norm_gs)
                norm_gs = zeros(1, self.iter_max);
            end

            if (record_values)
                values = zeros(1, self.iter_max);
            end

            while (norm_g > self.tol && iter < self.iter_max - 1)
                alpha = self.s;
                g_mem = change_mat_circularly(g_mem, g, iter);
                D_k = get_D_k(g_mem, self.eposilon);
                x_next = x - alpha * D_k * g;

                while (obj.f(x_next) - obj.f(x) > -self.gamma * alpha * g' * D_k * g)
                    alpha = alpha * self.sigma;
                    x_next = x - alpha * D_k * g;
                end

                if (record_steps)
                    steps(:,iter+1) = x;
                end

                if (record_norm_gs)
                    norm_gs(:, iter+1) = norm_g;
                end

                if (record_values)
                    values(:, iter+1) = obj.f(x);
                end

                x = x_next;
                g = obj.g(x);
                norm_g = norm(g, 2);
                iter = iter + 1;

            end

            % after one step, the iter becomes one

            t = toc;

            x_p = x;
            value_p = obj.f(x_p);
            norm_p = norm_g;
            success = norm_g <= self.tol;
            iter = iter+1;

            if (record_steps)
                steps(:, iter) = x;
                steps = steps(:, 1:iter);
            else
                steps = [];
            end

            if (record_norm_gs)
                norm_gs(:,iter) = norm_g;
                norm_gs = norm_gs(:,1:iter);
            else
                norm_gs = [];
            end

            if (record_values)
                norm_gs(:, iter) = value;
                values = values(:,1:iter);
            else
                values = [];
            end
        end

        function str = self_info(self)
            str = sprintf('---- Method Info ---> name: %s, tol: %f, iter_max: %d, s: %f, sigma: %f, gamma: %f, eposilon: %f, m: %d ----\n', self.name, self.tol, self.iter_max, self.s, self.sigma, self.gamma, self.eposilon, self.m);
        end

    end

end

function mat = change_mat_circularly(mat, new_v, iter)
    mat(:, mod(iter, length(mat)) + 1) = new_v;
end
