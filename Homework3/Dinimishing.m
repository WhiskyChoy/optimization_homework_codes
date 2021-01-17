classdef Dinimishing < Solver

    properties
        name = 'Dinimishing'
        alpha_func
    end

    methods

        function self = Dinimishing(alpha_func, tol, iter_max)

            arguments
                alpha_func (1, 1) = @(k)(1 / (1 + k));
                tol (1, 1) double = 1e-5
                iter_max (1, 1) int32 = 1e4
            end

            self = self@Solver(tol, iter_max);
            self.alpha_func = alpha_func;
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
                % In matlab the "=" between two matrix will pass the value. No need to worry about deep copy.
                alpha = self.alpha_func(iter);

                if (record_steps)
                    steps(:, iter + 1) = x;
                end

                if (record_norm_gs)
                    norm_gs(:, iter+1) = norm_g;
                end

                if (record_values)
                    values(:, iter+1) = obj.f(x);
                end

                x = x - alpha * g;
                g = obj.g(x);
                norm_g = norm(g, 2);
                iter = iter + 1;
            end

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
            str = sprintf('---- Method Info ---> name: %s, tol: %f, iter_max: %d, alpha_func: %s ----\n', self.name, self.tol, self.iter_max, func2str(self.alpha_func));
        end

    end

end
