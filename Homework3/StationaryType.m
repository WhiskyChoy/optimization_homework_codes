classdef StationaryType

    properties
        value (1, 1) string
    end

    methods

        function self = StationaryType(value)
            self.value = value;
        end

        % This is the converter. I don't know why it doesn't work.
        % function str = string(StationaryTypeObj)
        %     str = string(StationaryTypeObj.value);
        % end

    end

    enumeration
        Unknown ('unknown')
        NeedAly ('need further analysis')
        GlobalMin ('global minimizer')
        LocalMin ('local minimizer')
        GlobalMax ('global maximizer')
        LocalMax ('local maximizer')
        SaddlePoint ('saddle point')
        Error ('error exists')
    end

end
