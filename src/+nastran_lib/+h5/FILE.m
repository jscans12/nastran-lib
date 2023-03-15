classdef FILE < handle
%FILE Access to file-level information
    properties (SetAccess = immutable, GetAccess = protected)
        %H5_PATH Path to the NASTRAN h5 file
        h5_path
        %H5_file h5 file object for results access
        h5_file
    end
    properties (Access = private)
        INDEX_
        NASTRAN_
    end
    properties (Dependent)
        %INDEX The file index, as an h5 object
        INDEX
        %NASTRAN Access to NASTRAN simulation data
        NASTRAN
    end
    methods
        % Getters and Setters
        function INDEX = get.INDEX(this)
            if isempty(this.INDEX_)
                this.INDEX_ = this.h5_file.get_group('INDEX');
            end
            INDEX = this.INDEX_;
        end
        function NASTRAN = get.NASTRAN(this)
            if isempty(this.NASTRAN_)
                this.NASTRAN_ = nastran_lib.h5.NASTRAN(this.h5_file.get_group('NASTRAN'));
            end
            NASTRAN = this.NASTRAN_;
        end
        % Constructor
        function this = FILE(h5_path)
        %FILE Class constructor
            
            % Set file data to the object
            this.h5_path = h5_path;
            this.h5_file = h5io.file(h5_path);
            
        end
    end
end

