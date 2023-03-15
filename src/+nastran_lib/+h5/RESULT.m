classdef RESULT < handle
%RESULT Access to NASTRAN simulation results
    properties (SetAccess = immutable, GetAccess = private)
        %H5_group h5 group object for results access
        h5_group
        %SOLUTION The solution sequence
        SOLUTION
    end
    properties (SetAccess = immutable)
        %OUTPUTS structure of output data
        OUTPUTS
    end
    methods
        % Constructor
        function this = RESULT(h5_group,solver_id)
        %RESULT Class constructor
        
            % Set immutable data
            this.h5_group = h5_group;
            
            % Get the domains
            DOMAINS = struct2table(this.h5_group.get_dataset('DOMAINS').get_data);
            
            % Detect solver type
            this.SOLUTION = solver_id;
            
            % Set the OUTPUTS struct
            this.OUTPUTS = struct;
            for i = 1:numel(this.h5_group.group_names)
                group_name = this.h5_group.group_names{i};
                this.OUTPUTS.(group_name) = group_to_struct(this.h5_group.get_group(group_name),...
                                                            DOMAINS,...
                                                            this.SOLUTION);
            end
            
        end
    end
end
function struct_out = group_to_struct(group_in,DOMAINS,SOLUTION)
%GROUP_TO_STRUCT Recursive formula to convert a results group into
%a structure
    struct_out = struct;
    for i = 1:numel(group_in.group_names)
        group_name = group_in.group_names{i};
        struct_out.(group_name) = group_to_struct(group_in.get_group(group_name),DOMAINS,SOLUTION);
    end
    for i = 1:numel(group_in.dataset_names)
        NAME = group_in.dataset_names{i};
        struct_out.(NAME) = nastran_lib.h5.OUTPUT(group_in.get_dataset(NAME),...
                                                  DOMAINS,...
                                                  SOLUTION,...
                                                  NAME);
    end
end

