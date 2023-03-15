classdef OUTPUT < handle
%OUTPUT Container to hold an output dataset
    properties (SetAccess = immutable, GetAccess = private)
        %DOMAINS Domain info
        DOMAINS
        %DATASET The h5 dataset object
        dataset
    end
    properties (SetAccess = immutable)
        %SOLUTION The solution sequence used to generate this data
        SOLUTION
        %NAME The name of the dataset
        NAME
    end
    properties (SetAccess = protected)
        %DATA data as buffered
        data
    end
    properties (Dependent)
        %BUFFERED Flag set for buffered data
        buffered
    end
    properties
        %CONDENSE_RANGE Range of IDs to condense. Condensing assumes the
        %final digit of an ID is the DOF
        condense_range int64
    end
    methods
        % Getters and setters
        function buffered = get.buffered(this)
            buffered = ~isempty(this.data);
        end
        function set.condense_range(this,condense_range)
            if ~isequal(condense_range,this.condense_range)
                if ~isempty(condense_range)
                    if ~isequal(size(condense_range),[2,1])
                        error('condense_range must be 2 element column vector');
                    end
                end
                this.flush;
                this.condense_range = condense_range;
            end
        end
        % Constructor
        function this = OUTPUT(dataset,DOMAINS,SOLUTION,NAME)
        %OUTPUT Class constructor
            
            % Set properties
            this.dataset  = dataset;
            this.DOMAINS  = DOMAINS;
            this.SOLUTION = SOLUTION;
            this.NAME     = NAME;
            
        end
        % Data access
        function table_out = summarize(this,TYPE,SUBCASE,SE,ID,DOF,varargin)
        %SUMMARIZE Summarize a downselection of this output
        %
        %INPUTS:
        %------
        %
        %TYPE     - The type of output, support for solvers listed below
        %               SOL 101: 'max','min','absmax','rss','vonmises'
        %               SOL 111: 'max','min','absmax','rss','time','psd','vrs','vonmises'
        %               SOL 112: 'max','min','absmax','rss','time','psd','srs','vonmises'
        %SUBCASE  - subcase ID(s), empty for all
        %SE       - superelement ID(s), empty for all
        %ID       - entity ID(s), empty for all
        %DOF      - degree(s) of freedom, empty for all. RSS and von mises
        %           require an input. RSS can be entered as a cellstr
        %           of 1-2 RSS DOF, or a cell array of multiple cellstrs
        %           for multiple RSS requests
        %varargin - extra arguments for some output types, * if required
        %           srs: freq*, Q*
        %           vrs: freq*, Q*
        %           psd  (with SOL 112): window, overlap, trace
        %           time (with SOL 111): duration*, n_segments*
        %           max/min/absmax (with SOL 111): sigma (default = 3)
            
            % Make sure the object is buffered
            this.buffer;
        
            % Sometimes ID is called ID, sometimes EID (for elements)
            ID_name = 'ID';
            if any(strcmpi('EID',this.data.Properties.VariableNames))
                ID_name = 'EID';
            end
            
            % Input checks
            if ~ischar(TYPE)
                error('TYPE must be a character vector');
            end
            isRSS = strcmpi(TYPE,'rss');
            
            % Detect DOF names
            DOF_names = setdiff(this.data.Properties.VariableNames,{ID_name,'SUBCASE','TIME_FREQ_EIGR','SE','RANDOM'});
            
            % Default inputs
            if isempty(SUBCASE)
                SUBCASE = unique(this.data.SUBCASE);
            end
            if isempty(SE)
                SE = unique(this.data.SE);
            end
            if isempty(ID)
                ID = unique(this.data.(ID_name));
            end
            if isempty(DOF)
                switch lower(TYPE)
                    case {'rss','vonmises'}
                        error('User must specify DOF for this output');
                    otherwise
                        DOF = DOF_names;
                end
            end
            
            % Input checks
            if ~isnumeric(SUBCASE) || ~isvector(SUBCASE)
                error('SUBCASE must be a numeric vector');
            end
            if ~isnumeric(SE) || ~isvector(SE)
                error('SE must be a numeric vector');
            end
            if ~isnumeric(ID) || ~isvector(ID)
                error('ID must be a numeric vector');
            end
            if ischar(DOF)
                DOF = {DOF};
            elseif ~isRSS && ~iscellstr(DOF) %#ok<ISCLSTR>
                error('DOF must be a char or a cellstr');
            end
            if isRSS
                if iscellstr(DOF) %#ok<ISCLSTR>
                    DOF = {DOF};
                end
                if ~all(cellfun(@iscellstr,DOF))
                    error('DOF must be cell of cellstrs');
                end
            end
            
            % Parse optional arguments
            optional_args = struct;
            for i = 1:2:numel(varargin)
                arg_name = lower(varargin{i});
                arg_val  = varargin{i+1};
                if ~ischar(arg_name)
                    error('Optional arguments are out of order');
                end
                optional_args.(arg_name) = arg_val;
            end
            
            % Ensure uniqueness of data
            SUBCASE = unique(SUBCASE);
            SE      = unique(SE);
            ID      = unique(ID);
            if ~isRSS
                DOF = unique(DOF);
            end
            
            % Get indices of rows used
            rows_used = ismember(this.data.SUBCASE,SUBCASE) ...
                      & ismember(this.data.SE,SE) ...
                      & ismember(this.data.(ID_name),ID);
            if this.SOLUTION == 111
                rows_used = rows_used & this.data.RANDOM <= 1;
            end
            
            % Get groups for pivoted table data
            [row_groups,table_base] = findgroups(this.data(rows_used,{'SUBCASE','SE',ID_name}));
            table_tmp = table(repmat({this.NAME},height(table_base),1),repmat({TYPE},height(table_base),1),'VariableNames',{'OUTPUT','TYPE'});
            table_base = [table_base table_tmp];
            
            % Instantiate output table
            table_out = cell2table(cell(0,7),'VariableNames',{'SUBCASE','SE',ID_name,'OUTPUT','TYPE','DOF','VALUE'});
            
            % DOF loop
            for i = 1:numel(DOF)
                
                % Current DOF
                iDOF = DOF{i};
                
                % Methods based on solution type
                switch this.SOLUTION
                    
                    % Static
                    case 101
                        
                        % Output based on type requested
                        switch lower(TYPE)
                            case 'max'
                                rowFn = @(x) x;
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'min'
                                rowFn = @(x) x;
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'absmax'
                                rowFn = @(x) abs(x);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'rss'
                                switch numel(iDOF)
                                    case 2
                                        rowFn = @(x,y) sqrt(x.^2+y.^2);
                                        iData = splitapply(rowFn,this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),row_groups);
                                    case 3
                                        rowFn = @(x,y,z) sqrt(x.^2+y.^2+z.^2);
                                        iData = splitapply(rowFn,this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),...
                                                                 this.data.(iDOF{3})(rows_used),row_groups);
                                    otherwise
                                        error('Only 2 and 3 DOF RSS''s are supported');
                                end
                            case 'vonmises'
                                switch lower(this.NAME)
                                    case {'tria3','quad4'}
                                        rowFn = @(Sx,Sy,Txy) nastran_lib.references.vonMises_time(Sx,Sy,0,Txy,0,0);
                                        switch lower(iDOF)
                                            case 'top'
                                                XYname = 'XY1';
                                                if any(strcmpi('TXY1',this.data.Properties.VariableNames))
                                                    XYname = 'TXY1';
                                                end
                                                iData = splitapply(rowFn,this.data.X1(rows_used),...
                                                                         this.data.Y1(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            case 'bottom'
                                                XYname = 'XY2';
                                                if any(strcmpi('TXY2',this.data.Properties.VariableNames))
                                                    XYname = 'TXY2';
                                                end
                                                iData = splitapply(rowFn,this.data.X2(rows_used),...
                                                                         this.data.Y2(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            otherwise
                                                error('%s is not an acceptable Von Mises DOF, must be TOP or BOTTOM',iDOF);
                                        end
                                    otherwise
                                        error('Element type %s not coded for Von Mises stress',this.NAME);
                                end
                            otherwise
                                error('Output type %s is unknown or unsupported for SOL %d',TYPE,this.SOLUTION);
                        end

                    % Modal frequency response
                    case 111
                        
                        % Set sigma
                        sigma = 3;
                        if isfield(optional_args,'window')
                            sigma = optional_args.sigma;
                        end
                        
                        % PSDF used to calculate responses from transfer
                        % functions
                        if isfield(optional_args,'psdf')
                            if ~isa(optional_args.psdf,'nastran_lib.references.types.psd')
                                error('Optional PSDF input must be of type nastran_lib.references.types.psd');
                            end
                            freq = unique(this.data.TIME_FREQ_EIGR(rows_used));
                                psd_interp = nastran_lib.references.interp1_log(optional_args.psdf.Frequency,...
                                                                                optional_args.psdf.Data,...
                                                                                freq);
                                psd_interp(isnan(psd_interp)) = 0;
                                PSDF = nastran_lib.references.types.psd(freq,psd_interp);
                        else
                            PSDF = [];
                        end
                        
                        % Output based on type requested
                        switch lower(TYPE)
                            case 'max'
                                if isempty(PSDF)
                                    rowFn = @(x,y)  sigma*nastran_lib.references.types.psd(x,y).RMS;
                                else
                                    rowFn = @(x,y)  sigma*nastran_lib.references.types.psd(x,(abs(y).^2).*PSDF.Data).RMS;
                                end
                                iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                         this.data.(iDOF)(rows_used),row_groups);
                            case 'min'
                                if isempty(PSDF)
                                    rowFn = @(x,y) -sigma*nastran_lib.references.types.psd(x,y).RMS;
                                else
                                    rowFn = @(x,y) -sigma*nastran_lib.references.types.psd(x,(abs(y).^2).*PSDF.Data).RMS;
                                end
                                iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                         this.data.(iDOF)(rows_used),row_groups);
                            case 'absmax'
                                if isempty(PSDF)
                                    rowFn = @(x,y)  sigma*nastran_lib.references.types.psd(x,y).RMS;
                                else
                                    rowFn = @(x,y)  sigma*nastran_lib.references.types.psd(x,(abs(y).^2).*PSDF.Data).RMS;
                                end
                                iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                         this.data.(iDOF)(rows_used),row_groups);
                            case 'rss'
                                switch numel(iDOF)
                                    case 2
                                        if isempty(PSDF)
                                            rowFn = @(freq,x,y) sqrt((sigma*nastran_lib.references.types.psd(freq,x).RMS).^2 ...
                                                                   + (sigma*nastran_lib.references.types.psd(freq,y).RMS).^2);
                                        else
                                            rowFn = @(freq,x,y) sqrt((sigma*nastran_lib.references.types.psd(freq,(abs(x).^2).*PSDF.Data).RMS).^2 ...
                                                                   + (sigma*nastran_lib.references.types.psd(freq,(abs(y).^2).*PSDF.Data).RMS).^2);
                                        end
                                        iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                                 this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),row_groups);
                                    case 3
                                        if isempty(PSDF)
                                            rowFn = @(freq,x,y,z) sqrt((sigma*nastran_lib.references.types.psd(freq,x).RMS).^2 ...
                                                                     + (sigma*nastran_lib.references.types.psd(freq,y).RMS).^2 ...
                                                                     + (sigma*nastran_lib.references.types.psd(freq,z).RMS).^2);
                                        else
                                            rowFn = @(freq,x,y,z) sqrt((sigma*nastran_lib.references.types.psd(freq,(abs(x).^2).*PSDF.Data).RMS).^2 ...
                                                                     + (sigma*nastran_lib.references.types.psd(freq,(abs(y).^2).*PSDF.Data).RMS).^2 ...
                                                                     + (sigma*nastran_lib.references.types.psd(freq,(abs(z).^2).*PSDF.Data).RMS).^2);
                                        end
                                        iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                                 this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),...
                                                                 this.data.(iDOF{3})(rows_used),row_groups);
                                    otherwise
                                        error('Only 2 and 3 DOF RSS''s are supported');
                                end
                            case 'vrs'
                                error('Not coded');
                                % https://vibrationdata.wordpress.com/2012/10/15/sdof-response-to-an-acceleration-psd-base-input/
                            case 'psd'
                                if isempty(PSDF)
                                    rowFn = @(x,y) nastran_lib.references.types.psd(x,abs(y));
                                else
                                    rowFn = @(x,y) nastran_lib.references.types.psd(x,(abs(y).^2).*PSDF.Data);
                                end
                                iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                         this.data.(iDOF)(rows_used),row_groups);
                            case 'time'
                                if isempty(PSDF)
                                    rowFn = @(x,y) nastran_lib.references.types.psd(x,y).timehistory(optional_args.duration,optional_args.n_segments);
                                else
                                    rowFn = @(x,y) nastran_lib.references.types.psd(x,(abs(y).^2).*PSDF.Data).timehistory(optional_args.duration,optional_args.n_segments);
                                end
                                iData = splitapply(rowFn,this.data.TIME_FREQ_EIGR(rows_used),...
                                                         this.data.(iDOF)(rows_used),row_groups);
                            case 'vonmises'
                                switch lower(this.NAME)
                                    case {'tria3_cplx','quad4_cplx'}
                                        rowFn = @(Sx,Sy,Txy) sigma*nastran_lib.references.types.psd(freq,...
                                                                                                    nastran_lib.references.vonMises_freq(Sx,...
                                                                                                                                         Sy,...
                                                                                                                                         [],...
                                                                                                                                         Txy,...
                                                                                                                                         [],...
                                                                                                                                         [],...
                                                                                                                                         PSDF.Data)).RMS;
                                        switch lower(iDOF)
                                            case 'top'
                                                XYname = 'XY1';
                                                if any(strcmpi('TXY1',this.data.Properties.VariableNames))
                                                    XYname = 'TXY1';
                                                end
                                                iData = splitapply(rowFn,this.data.X1(rows_used),...
                                                                         this.data.Y1(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            case 'bottom'
                                                XYname = 'XY2';
                                                if any(strcmpi('TXY2',this.data.Properties.VariableNames))
                                                    XYname = 'TXY2';
                                                end
                                                iData = splitapply(rowFn,this.data.X2(rows_used),...
                                                                         this.data.Y2(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            otherwise
                                                error('%s is not an acceptable Von Mises DOF, must be TOP or BOTTOM',iDOF);
                                        end
                                    otherwise
                                        error('Element type %s not coded for Von Mises stress',this.NAME);
                                end
                            otherwise
                                error('Output type %s is unknown or unsupported for SOL %d',TYPE,this.SOLUTION);
                        end

                    % Modal transient response
                    case 112

                        % Output based on type requested
                        switch lower(TYPE)
                            case 'max'
                                rowFn = @(x) max(x);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'min'
                                rowFn = @(x) min(x);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'absmax'
                                rowFn = @(x) max(abs(x));
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),row_groups);
                            case 'rss'
                                switch numel(iDOF)
                                    case 2
                                        rowFn = @(x,y) max(sqrt(x.^2+y.^2));
                                        iData = splitapply(rowFn,this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),row_groups);
                                    case 3
                                        rowFn = @(x,y,z) max(sqrt(x.^2+y.^2+z.^2));
                                        iData = splitapply(rowFn,this.data.(iDOF{1})(rows_used),...
                                                                 this.data.(iDOF{2})(rows_used),...
                                                                 this.data.(iDOF{3})(rows_used),row_groups);
                                    otherwise
                                        error('Only 2 and 3 DOF RSS''s are supported');
                                end
                            case 'srs'
                                rowFn = @(x,y) this.calc_srs(nastran_lib.references.types.timehistory(x,y),optional_args.freq,...
                                                                                                           optional_args.q);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),...
                                                         this.data.TIME_FREQ_EIGR(rows_used),row_groups);
                            case 'psd'
                                if ~isfield(optional_args,'window')
                                    optional_args.window = [];
                                end
                                if ~isfield(optional_args,'overlap')
                                    optional_args.overlap = 0.5;
                                end
                                if ~isfield(optional_args,'trace')
                                    optional_args.trace = 'mean';
                                end
                                rowFn = @(x,y) nastran_lib.references.types.timehistory(x,y).psd(optional_args.window,...
                                                                                                 optional_args.overlap,...
                                                                                                 optional_args.trace);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),...
                                                         this.data.TIME_FREQ_EIGR(rows_used),row_groups);
                            case 'time'
                                rowFn = @(x,y) nastran_lib.references.types.timehistory(x,y);
                                iData = splitapply(rowFn,this.data.(iDOF)(rows_used),...
                                                         this.data.TIME_FREQ_EIGR(rows_used),row_groups);
                            case 'vonmises'
                                switch lower(this.NAME)
                                    case {'tria3','quad4'}
                                        rowFn = @(Sx,Sy,Txy) max(nastran_lib.references.vonMises_time(Sx,Sy,0,Txy,0,0));
                                        switch lower(iDOF)
                                            case 'top'
                                                XYname = 'XY1';
                                                if any(strcmpi('TXY1',this.data.Properties.VariableNames))
                                                    XYname = 'TXY1';
                                                end
                                                iData = splitapply(rowFn,this.data.X1(rows_used),...
                                                                         this.data.Y1(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            case 'bottom'
                                                XYname = 'XY2';
                                                if any(strcmpi('TXY2',this.data.Properties.VariableNames))
                                                    XYname = 'TXY2';
                                                end
                                                iData = splitapply(rowFn,this.data.X2(rows_used),...
                                                                         this.data.Y2(rows_used),...
                                                                         this.data.(XYname)(rows_used),row_groups);
                                            otherwise
                                                error('%s is not an acceptable Von Mises DOF, must be TOP or BOTTOM',iDOF);
                                        end
                                    otherwise
                                        error('Element type %s not coded for Von Mises stress',this.NAME);
                                end
                            otherwise
                                error('Output type %s is unknown or unsupported for SOL %d',TYPE,this.SOLUTION);
                        end

                    otherwise
                        error('Outputs are not yet supported for SOL %d',this.SOLUTION);

                end
                
                % And append to the table
                if isRSS
                    dofName = sprintf('%s,',iDOF{:});
                    dofName = dofName(1:end-1);
                else
                    dofName = iDOF;
                end
                table_tmp = table(repmat({dofName},numel(iData),1),iData,'VariableNames',{'DOF','VALUE'});
                table_tmp = [table_base table_tmp]; %#ok<AGROW>
                
                % Append to output
                table_out = [table_out;table_tmp]; %#ok<AGROW>

            end
            
        end
        % Buffering
        function buffer(this)
        %BUFFER Buffer the data
        
            % Don't buffer if it's already buffered
            if this.buffered
                return
            end
            
            % Get the solution data
            solution_data = this.dataset.get_data;
            
            % Special processing for complex numbers
            if contains(this.NAME,'_cplx','IgnoreCase',true)
                fnames = fieldnames(solution_data);
                for i = 1:numel(fnames)
                    fname = fnames{i};
                    if numel(fname) > 1 && strcmpi(fname(end),'r')
                        new_name = fname(1:end-1);
                        imag_name = sprintf('%sI',new_name);
                        imag_indx = strcmpi(imag_name,fnames);
                        if any(imag_indx)
                            solution_data.(fname) = complex(solution_data.(fname),solution_data.(imag_name));
                            solution_data = rmfield(solution_data,imag_name);
                            solution_data = renameStructField(solution_data,fname,new_name);
                        end
                    end
                end
            end
            
            % Deal with condensed ranges
            if ~isempty(this.condense_range)
                if ~isequal(fieldnames(solution_data),{'ID','X','Y','Z','RX','RY','RZ','DOMAIN_ID'}')
                    error('Condensing currently only defined for nodal motion');
                end
                
                % Index where IDs are condensed
                condense_index = solution_data.ID >= this.condense_range(1) & solution_data.ID <= this.condense_range(2);
                
                % New condensed IDs
                condense_ids = int64(zeros(size(solution_data.ID)));
                condense_ids(condense_index) = int64(floor(double(solution_data.ID(condense_index))/10));
                
                % New condensed DOF
                condense_dof = int64(zeros(size(solution_data.ID)));
                condense_dof(condense_index) = solution_data.ID(condense_index) - condense_ids(condense_index)*10;
                
                % Get updated ID and domain for new data
                new_data.ID        = condense_ids(condense_dof == 1);
                new_data.X         = solution_data.X(condense_dof == 1);
                new_data.Y         = solution_data.X(condense_dof == 2);
                new_data.Z         = solution_data.X(condense_dof == 3);
                new_data.RX        = solution_data.X(condense_dof == 4);
                new_data.RY        = solution_data.X(condense_dof == 5);
                new_data.RZ        = solution_data.X(condense_dof == 6);
                new_data.DOMAIN_ID = solution_data.DOMAIN_ID(condense_dof == 1);
                
                % Update struct
                solution_data.ID        = [solution_data.ID(~condense_index);new_data.ID];
                solution_data.X         = [solution_data.X(~condense_index);new_data.X];
                solution_data.Y         = [solution_data.Y(~condense_index);new_data.Y];
                solution_data.Z         = [solution_data.Z(~condense_index);new_data.Z];
                solution_data.RX        = [solution_data.RX(~condense_index);new_data.RX];
                solution_data.RY        = [solution_data.RY(~condense_index);new_data.RY];
                solution_data.RZ        = [solution_data.RZ(~condense_index);new_data.RZ];
                solution_data.DOMAIN_ID = [solution_data.DOMAIN_ID(~condense_index);new_data.DOMAIN_ID];
                
            end
            
            % Convert to a table
            solution_data = struct2table(solution_data);
            
            % Output is based on the solution sequence
            switch this.SOLUTION
                
                % Linear statics
                case 101
                    
                    DOMAINS_used = this.DOMAINS(:,{'ID','SUBCASE','SE'});
                    
                % Modal frequency and transient response
                case {111,112}
                    
                    rows_keep = this.DOMAINS.MODE == 0;
                    DOMAINS_used = this.DOMAINS(rows_keep,{'ID','SUBCASE','TIME_FREQ_EIGR','SE','RANDOM'});
                    
                otherwise
                    
                    DOMAINS_used = this.DOMAINS;
                    
            end
            
            % Table join
            data_dom_col = find(strcmpi('DOMAIN_ID',solution_data.Properties.VariableNames));
            dom_dom_col  = find(strcmpi('ID',DOMAINS_used.Properties.VariableNames));
            out = innerjoin(solution_data,DOMAINS_used,...
                            'LeftKeys',data_dom_col,...
                            'RightKeys',dom_dom_col);
            out = removevars(out,'DOMAIN_ID');
            
            % Buffer the data
            this.data = out;
            
        end
        function flush(this)
        %FLUSH Flush the data buffer
        
            % Just clear the data out
            this.data = [];
            
        end
    end
    methods (Static,Access = private)
        function srs = calc_srs(th,freq,q)
            
            srs      = struct;
            if isrow(freq)
                srs.freq = freq';
            else
                srs.freq = freq;
            end
            srs.SRS  = th.srs(q,freq);
            
        end
    end
end

