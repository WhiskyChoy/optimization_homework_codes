classdef ObjectiveFunction

    methods (Abstract)
        res = f(self, x)
        res = g(self, x)
        res = h(self, x)
    end
end
