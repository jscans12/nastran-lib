classdef matrix < handle
%MATRIX An OP4 matrix, see "DMAP Programmers Guide" for more info
%
%Supports the most common formats:
% - sparse binary string-header
% - sparse binary regular string (BIGMAT)
% - sparse ASCII string-header
%
%Unsupported:
% - nonsparse binary format
% - nonsparse ASCII format
% - sparse ASCII regular string (BIGMAT)
%
%WARNING | Other limitations:
% - Not thoroughly tested, use at your own risk
% - Exponent in floats must have a letter (e.g. 1.0E+3)
% - Only set up to read superelement OP4 files
% - Only real numbers
% - Stores all floats in double precision (MATLAB limitation)
    properties (SetAccess = private)
        %NAME Matrix name
        name char
        %VALUE Stored value of the matrix
        value
        %NF Form of the matrix
        NF(1,1) int32
        %NTYPE Type of the matrix
        NTYPE(1,1) int32
        %BIGMAT Matrix is stored as BIGMAT
        BIGMAT(1,1) logical
        %NUM_FORMAT Numeric format info
        num_format
    end
    properties (Dependent)
        %NCOL Number of columns in the matrix
        NCOL
        %NR Number of rows in the matrix
        NR
        %IS_SPARSE The matrix is stored sparse
        is_sparse
    end
    methods
        % Getters and setters
        function NCOL = get.NCOL(this)
            NCOL = size(this.value,2);
        end
        function NR = get.NR(this)
            NR = size(this.value,1);
        end
        function is_sparse = get.is_sparse(this)
            is_sparse = issparse(this.value);
        end
        % Populate from scratch
        function populate(this,name,value,NF,NTYPE,BIGMAT,num_format)
        %POPULATE Fill the object with existing data
            
            this.name       = name;
            this.value      = value;
            this.NF         = NF;
            this.NTYPE      = NTYPE;
            this.BIGMAT     = BIGMAT;
            this.num_format = num_format;
            
        end
        % Read/write
        function read(this,file_id,isbinary,position,origin)
        %READ Read a matrix from an OP4 file
            
            % Defaults
            if nargin < 3
                isbinary = false;
            end
            if nargin < 4
                position = ftell(file_id);
            end
            if nargin < 5
                origin = 'bof';
            end
            
            % Go to position
            fseek(file_id,position,origin);
            
            % Get the matrix header
            if isbinary
                
                % Parse
                nrec_start = fread(file_id,1,'int32');
                NCOL_      = fread(file_id,1,'int32');
                NR_        = fread(file_id,1,'int32');
                this.NF    = fread(file_id,1,'int32');
                this.NTYPE = fread(file_id,1,'int32');
                this.name  = strtrim(char(fread(file_id,8,'uchar'))');
                nrec_end   = fread(file_id,1,'int32');
                if nrec_start ~= nrec_end
                    error('Unknown file read error');
                end
                
                % Doesn't apply
                this.num_format = 'binary';
                
                % Check if complex
                if this.NTYPE == 1
                    error('Single precision binary is current unsupported');
                end
                
            else
            
                % Get record
                curr_line = int32(sscanf(fgetl(file_id),'%8d%8d%8d%8d%16c%16c'));

                % Parse
                this.name  = strtrim(char(curr_line(5:12))');
                NCOL_  = curr_line(1);
                NR_    = curr_line(2);
                this.NF    = curr_line(3);
                this.NTYPE = curr_line(4);
                FFMT  = strtrim(char(curr_line(12:end))');
                
                % Get FORTRAN format for ASCII string
                this.num_format = nastran_lib.references.fortran_ascii(FFMT);
                float_format = ['%' num2str(this.num_format.field_width) 'e'];
                replace_ex = ~strcmpi(this.num_format.exponent_char,'e');
            
            end

            % Check if complex
            if this.NTYPE > 2
                error('Complex numbers are current unsupported');
            end
            
            % Check if bigmat
            this.BIGMAT = false;
            if NR_ < 0
                this.BIGMAT = true;
                NR_ = -NR_;
            end
            
            % Allocate matrix
            this.value = zeros(double(NR_),double(NCOL_));
            
            % Loop through column entries
            while true

                % Different methods for binary and ASCII OP4s
                if isbinary
                    
                    % Read next record
                    nrec_start = fread(file_id,1,'int32');
                    ICOL       = fread(file_id,1,'int32');
                    IROW1      = fread(file_id,1,'int32');
                    NW         = fread(file_id,1,'int32');
                    
                    % Data size info
                    if ICOL > NCOL_
                        fread(file_id,1,'double');
                        fread(file_id,1,'int32');
                        break;
                    end
                    
                    % Check if sparse
                    is_sparse_ = false;
                    if IROW1 == 0
                        is_sparse_ = true;
                    end
                    
                    % Not supported yet
                    if ~is_sparse_
                        error('Nonsparse currently unsupported');
                    end
                    
                    % Loop through the row entries
                    curr_NW = 0;
                    while true
                        
                        % Support for BIGMAT and standard formats
                        if this.BIGMAT
                            
                            % Get L and IROW
                            L = fread(file_id,1,'int32') - 1;
                            IROW2 = fread(file_id,1,'int32');
                            curr_NW = curr_NW + 2;
                            
                        else
                        
                            % Sparse header thing
                            IS   = fread(file_id,1,'int32');
                            curr_NW = curr_NW + 1;

                            % Size per DMAP programmers guide
                            L = floor(double(IS)/65536) - 1;
                            IROW2 = IS - 65536 * (L + 1);
                        
                        end

                        % Read in values
                        nvalues = floor(double(L)/2);
                        values  = fread(file_id,nvalues,'double');
                        curr_NW = curr_NW + L;
                        
                        % Write to the matrix
                        this.value(IROW2:IROW2+numel(values)-1,ICOL) = values;
                        
                        % Break when end is reached
                        if curr_NW >= NW
                            if curr_NW ~= NW
                                error('Unknown file read error');
                            end
                            break
                        end
                    
                    end
                    
                    % Move on
                    nrec_end = fread(file_id,1,'int32');
                    if nrec_start ~= nrec_end
                        error('Unknown file read error');
                    end
                    
                else
                    
                    % Not supported yet
                    if this.BIGMAT
                        error('BIGMAT currently unsupported for ASCII files');
                    end
                    
                    % Get record
                    curr_line = fgetl(file_id);

                    % Switch depending field width
                    if numel(curr_line) > 16

                        % Get line
                        curr_line = int32(sscanf(curr_line,'%8d%8d%8d'));

                        % Data size info
                        ICOL = curr_line(1);
                        if ICOL > NCOL_
                            fgetl(file_id);
                            break;
                        end
                        IROW1 = curr_line(2);
                        NW   = curr_line(3); %#ok<NASGU>

                        % Check if sparse
                        is_sparse_ = false;
                        if IROW1 == 0
                            is_sparse_ = true;
                        end

                        % Not supported yet
                        if ~is_sparse_
                            error('Nonsparse currently unsupported');
                        end

                    else

                        % Get record
                        curr_line = int32(sscanf(curr_line,'%16d'));
                        IS = curr_line;

                        % Size per DMAP programmers guide
                        L = int32(double(IS)/65536) - 1;
                        IROW1 = IS - 65536 * (L + 1);

                        % Read all rows
                        rows_to_read = ceil(double(L)/double(this.num_format.repititions));
                        values = [];
                        exp_char = this.num_format.exponent_char;
                        for row_index = 1:rows_to_read
                            curr_line = fgetl(file_id);
                            if replace_ex
                                curr_line = strrep(curr_line,exp_char,'e');
                            end
                            curr_line = sscanf(curr_line,float_format);
                            values = [values; curr_line]; %#ok<*AGROW>
                        end
                        this.value(IROW1:IROW1+numel(values)-1,ICOL) = values;

                    end

                end

            end
            
            % I don't have a great feel for how much overhead a sparse
            % matrix requires in MATLAB, but for now let's say any matrix
            % under 50% denity is stored as sparse.
            if nnz(this.value) < 0.5 * NR_ * NCOL_
                this.value = sparse(this.value);
            end

        end
        function write(this,file_id)
        %WRITE Write to an OP4 file
            
            % Verify the format
            if this.BIGMAT
                error('BIGMAT currently unsupported');
            end
            if this.NTYPE > 2
                error('Complex numbers are current unsupported');
            end
            if ~this.is_sparse
                error('Nonsparse currently unsupported');
            end
            
            % Value write format
            value_format = ['%' num2str(this.num_format.field_width) '.' num2str(this.num_format.dec_places) 'E'];
            
            % Matrix header
            fprintf(file_id,'%8d%8d%8d%8d%-8s%-16s\n',this.NCOL,this.NR,this.NF,this.NTYPE,this.name,this.num_format.str_fmt);
            
            % Loop through columns
            for ICOL = 1:this.NCOL
                
                % Get current column value
                curr_col = this.value(:,ICOL);
                
                % Find start and end of nonsparse sections
                zero_rows = full(curr_col == 0);
                row_starts = find([~zero_rows(1); diff(zero_rows) < 0]);
                row_ends   = find([diff(zero_rows) > 0; ~zero_rows(end)]);
                
                % Check for code issues since I don't have 100% confidence
                % in myself here
                if numel(row_starts) ~= numel(row_ends)
                    error('Check code');
                end
                
                % Only write if nonzero
                if ~all(zero_rows)

                    % Print column header
                    fprintf(file_id,'%8d%8d%8d\n',ICOL,0,sum(~zero_rows)+numel(row_starts));

                    % Write each group of non-sparse rows
                    for row_find_ind = 1:numel(row_starts)

                        % Start and end of row (use DMAP nomenclature)
                        IROW    = row_starts(row_find_ind);
                        row_end = row_ends(row_find_ind);

                        % Calculate L and IS from DMAP Programmer's Guide
                        L = row_end - IROW + 1;
                        IS = IROW + 65536 * (L + 1);

                        % Row header
                        fprintf(file_id,'%11d\n',IS);

                        % Print each set of nonsparse sequential rows
                        value_start_indx = IROW;
                        while true
                            value_end_indx = min(row_end,value_start_indx + this.num_format.repititions - 1);
                            fprintf(file_id,value_format,full(curr_col(value_start_indx:value_end_indx)));
                            fprintf(file_id,'\n');
                            value_start_indx = value_end_indx + 1;
                            if value_start_indx > row_end
                                break
                            end
                        end

                    end
                    
                end
                
            end
            
            % Print trailer, not really sure this is right
            fprintf(file_id,'%8d%8d%8d\n',this.NCOL+1,1,1);
            fprintf(file_id,value_format,1);
            fprintf(file_id,'\n');
        
        end
        % Other
        function varargout = spy(this,varargin)
        %SPY Visualize sparsity pattern of matrix
        %
        %See Also: spy
            
            % Plot
            h = figure;
            spy(this.value,varargin{:});
            title(this.name);
            grid('on');
            
            % Output handle
            if nargout == 1
                varargout{1} = h;
            end
            if nargout > 1
                error('Too many output arguments');
            end
            
        end
    end
end

