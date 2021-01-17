classdef Backtracking < GradientMethod

    properties
        name = 'Backtracking'
        s (1, 1) double = 1
        sigma (1, 1) double {Util.check_zero_to_one} = 0.5
        gamma (1, 1) double {Util.check_zero_to_one} = 0.1
    end

    methods

        function self = Backtracking(s, sigma, gamma, tol, iter_max)

            arguments
                s (1, 1) double = 1
                sigma (1, 1) double {Util.check_zero_to_one} = 0.5
                gamma (1, 1) double {Util.check_zero_to_one} = 0.1
                tol (1, 1) double = 1e-5
                iter_max (1, 1) int32 = 1e4
            end

            self = self@GradientMethod(tol, iter_max);
            self.s = s;
            self.sigma = sigma;
            self.gamma = gamma;
        end

        function [x1_p, x2_p, value_p, iter, t, success, steps] = descend(self, obj, init_x1, init_x2, record)

            arguments
                self
                obj (1, 1) ObjectiveFunction
                init_x1 (1, 1) double
                init_x2 (1, 1) double
                record = false
            end

            tic

            iter = 0;
            x = [init_x1; init_x2];
            g = obj.g(x(1), x(2));
            norm_x = norm(g, 2);

            if (record)
                steps = zeros(2, self.iter_max);
            end

            while (norm_x > self.tol && iter < self.iter_max)
                % In matlab the "=" between two matrix will pass the value. No need to worry about deep copy.
                alpha = self.s;
                x_next = x - alpha * g;

                while (obj.f(x_next(1), x_next(2)) - obj.f(x(1), x(2)) > -self.gamma * alpha * norm_x^2)
                    alpha = alpha * self.sigma;
                    x_next = x - alpha * g;
                end

                if (record)
                    steps(:,iter+1) = x;
                end

                x = x_next;
                g = obj.g(x(1), x(2));
                norm_x = norm(g, 2);
                iter = iter + 1;

            end


            t = toc;

            x1_p = x(1);
            x2_p = x(2);
            value_p = obj.f(x1_p, x2_p);
            success = norm_x <= self.tol;
            if (record)
                steps = steps(:, 1:iter);
            end
        end

        function str = self_info(self)
            str = sprintf('---- Method Info ---> name: %s, tol: %f, iter_max: %d, s: %f, sigma: %f, gamma: %f ----\n', self.name, self.tol, self.iter_max, self.s, self.sigma, self.gamma);
        end

    end

end
