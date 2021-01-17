classdef StationaryPoint

    properties
        x1 (1, 1) 
        x2 (1, 1) 
        value (1, 1) 
        float_x1 (1, 1) double
        float_x2 (1, 1) double
        float_value (1, 1) double
        type StationaryType
    end

    methods

        function self = StationaryPoint(x1, x2, value, float_x1, float_x2, float_value, s_type)
            self.x1 = x1;
            self.x2 = x2;
            self.value = value;
            self.float_x1 = float_x1;
            self.float_x2 = float_x2;
            self.float_value = float_value;
            self.type = s_type;
        end

    end

end
