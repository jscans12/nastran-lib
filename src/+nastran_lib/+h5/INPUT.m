classdef INPUT < handle
%INPUT Access to NASTRAN input deck info
    properties (SetAccess = immutable)
        %H5_group h5 group object for results access
        h5_group
    end
    methods
        function this = INPUT(h5_group)
        %INPUT Class constructor
        
            % Set immutable data
            this.h5_group = h5_group;
            
        end
    end
end

