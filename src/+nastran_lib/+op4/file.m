classdef file < handle
%FILE Read/Write access to OP4 files
    properties
        %FILENAME Name of the file being accessed
        filename char
        %MATRICES Matrix objects
        matrices nastran_lib.op4.matrix
        %BINARY Is this a binary OP4?
        binary logical
    end
    methods
        % Getters and setters
        function set.filename(this,filename)
            if ~ischar(filename)
                error('path must be a character array');
            end
            this.filename = filename;
            if exist(filename,'file')
                this.read;
            end
        end
        % Constructor
        function this = file(filename)
        %FILE Class constructor
        
            % Optionally set file name
            if nargin > 0
                this.filename = filename;
            end
            
        end
        % Read/write
        function read(this)
        %READ Reads an OP4 file
            
            % Input checks
            if isempty(this.filename)
                error('Must set filename to class before reading');
            end
            
            % Open the file
            file_id = fopen(this.filename,'r');
            if file_id == -1
                error('File not found: %s\n',this.filename);
            end

            % Initialize output
            matrices_ = [];
            
            % Try to detect whether or not the file is binary
            this.binary = false;
            bin_test = fread(file_id,1,'int32');
            if bin_test == 24
                this.binary = true;
                fprintf(1,'This file was detected as a binary OP4\n');
            else
                fprintf(1,'This file was detected as an ASCII OP4\n');
            end
            frewind(file_id);

            % Loop all lines
            while true
                
                % Loop condition for binary files
                if this.binary
                    curr_pos = ftell(file_id);
                    nrec_start = fread(file_id,1,'int32');
                    if isempty(nrec_start)
                        break
                    end
                    fseek(file_id,curr_pos,'bof');
                
                % Loop condition for ASCII
                else
                    if feof(file_id)
                        break
                    end
                end
                
                % Read a matrix
                curr_matrix = nastran_lib.op4.matrix;
                curr_matrix.read(file_id,this.binary);

                % Set to object array
                if isempty(matrices_)
                    matrices_ = curr_matrix;
                else
                    matrices_(end+1,1) = curr_matrix; %#ok<AGROW>
                end
                
            end
            
            % Set to class
            this.matrices = matrices_;

            % Close the file
            fclose(file_id);
            
        end
        function write(this,matrices,permission)
        %WRITE Writes an OP4 file
        
            % Check if binary
            if this.binary
                error('Binary file write is not yet coded')
            end
            
            % Input checks
            if isempty(this.filename)
                error('Must set filename to class before reading');
            end
            if ~isa(matrices,'nastran_lib.op4.matrix')
                error('Matrices must be of type nastran_lib.op4.matrix');
            end
            if nargin < 3
                permission = 'W'; 
            end
            
            % Open the file
            file_id = fopen(this.filename,permission);
            if file_id == -1
                error('File not found: %s\n',this.filename);
            end
            
            % Write each matrix to the file
            for i = 1:numel(matrices)
                matrices(i).write(file_id);
            end
            
            % Close the file
            fclose(file_id);
            
            % Refresh the object in memory
            this.read;
            
        end
    end
end

