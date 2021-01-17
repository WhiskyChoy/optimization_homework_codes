classdef GlobalizedNewton < Backtracking

    properties
        beta1 (1, 1) double = 1e-6
        beta2 (1, 1) double = 1e-6
        p (1, 1) double = .1
    end

    methods

        function self = GlobalizedNewton(beta1, beta2, p, s, sigma, gamma, tol, iter_max)

            arguments
                beta1 (1, 1) double = 1e-6
                beta2 (1, 1) double = 1e-6
                p (1, 1) double = .1
                s (1, 1) double = 1
                sigma (1, 1) double {Util.check_zero_to_one} = 0.5
                gamma (1, 1) double {Util.check_zero_to_one} = 0.1
                tol (1, 1) double = 1e-5
                iter_max (1, 1) int32 = 1e4
            end

            self = self@Backtracking(s, sigma, gamma, tol, iter_max);
            self.name = 'Globalized Newton Method';
            self.beta1 = beta1;
            self.beta2 = beta2;
            self.p = p;
        end

        function [x_p, value_p, norm_p, iter, t, success, steps, norm_gs, values, use_newtons] = descend(self, obj, init_x, record_steps, record_norm_gs, record_values)

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
            [d, use_newton] = get_d(obj, self.beta1, self.beta2, self.p, x);

            if (record_steps)
                steps = zeros(size_x, self.iter_max);
            end

            if (record_norm_gs)
                norm_gs = zeros(1, self.iter_max);
            end

            if (record_values)
                values = zeros(1, self.iter_max);
            end

            use_newtons = zeros(1, self.iter_max);

            while (norm_g > self.tol && iter < self.iter_max - 1)

                alpha = self.s;
                x_next = x + alpha * d;

                while (obj.f(x_next) - obj.f(x) > self.gamma * alpha * g' * d)
                    alpha = alpha * self.sigma;
                    x_next = x + alpha * d;
                end

                if (record_steps)
                    steps(:, iter + 1) = x;
                end

                if (record_norm_gs)
                    norm_gs(:, iter + 1) = norm_g;
                end

                if (record_values)
                    values(:, iter+1) = obj.f(x);
                end

                use_newtons(:, iter+1) = use_newton;

                x = x_next;
                g = obj.g(x);
                norm_g = norm(g, 2);
                [d, use_newton] = get_d(obj, self.beta1, self.beta2, self.p, x);

                iter = iter + 1;
            end

            t = toc;

            x_p = x;
            value_p = obj.f(x_p);
            norm_p = norm_g;
            iter = iter+1;



            success = norm_g <= self.tol;

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

            use_newtons(:, iter) = use_newton;
            use_newtons = use_newtons(:, 1:iter);

        end

        function str = self_info(self)
            str = sprintf('---- Method Info ---> name: %s, tol: %f, iter_max: %d, beta1: %f, beta2: %f, p: %f, s: %f, sigma: %f, gamma: %f----\n', self.name, self.tol, self.iter_max, self.beta1, self.beta2, self.p, self.s, self.sigma, self.gamma);
        end

    end

end

function [d, use_newton]= get_d(obj, beta1, beta2, p, x)

    arguments
        obj (1, 1) ObjectiveFunction
        beta1 (1, 1) double
        beta2 (1, 1) double
        p (1, 1) double
        x double
    end

    % actually we should use hessian', but it's symmetric
    hessian = obj.h(x);
    g = obj.g(x);

    % A if full rank if det(A) not equals 0
    if (det(hessian) == 0)
        d = -g;
        use_newton = false;
        return;
    end

    % quicker than g * inv (hessian)
    % A \ b is to inverse A and 
    s = - hessian \ g;
    lhs = -g' * s;
    rhs = min(beta1, beta2 * norm(s, 2) ^ p) * norm(s, 2) ^ 2;

    if (lhs < rhs)
        d = -g;
        use_newton = false;
        return;
    end

    d = s;
    use_newton = true;
end
