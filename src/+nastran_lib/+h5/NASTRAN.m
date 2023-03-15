classdef NASTRAN < handle
%NASTRAN Access to NASTRAN simulation run
    properties (SetAccess = immutable, GetAccess = protected)
        %H5_group h5 group object for simulation data access
        h5_group
    end
    properties (SetAccess = immutable)
        %ATTRIBUTES Simulation attributes
        ATTRIBUTES
    end
    properties (Access = private)
        INPUT_
        RESULT_
    end
    properties (Dependent)
        %INPUT NASTRAN input deck information
        INPUT
        %RESULT NASTRAN results data
        RESULT
    end
    methods
        % Getters and Setters
        function INPUT = get.INPUT(this)
            if isempty(this.INPUT_)
                this.INPUT_ = nastran_lib.h5.INPUT(this.h5_group.get_group('INPUT'));
            end
            INPUT = this.INPUT_;
        end
        function RESULT = get.RESULT(this)
            if isempty(this.RESULT_)
                this.RESULT_ = nastran_lib.h5.RESULT(this.h5_group.get_group('RESULT'),this.ATTRIBUTES.SOL);
            end
            RESULT = this.RESULT_;
        end
        function this = NASTRAN(h5_group)
        %NASTRAN Class constructor
        
            % Set immutable data
            this.h5_group = h5_group;
            this.ATTRIBUTES = h5_group.attributes;
            
        end
    end
end

