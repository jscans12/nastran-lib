classdef deck < handle
%DECK The NASTRAN deck
    properties
        %FILE_MGMT File management section
        file_mgmt string
        %CASE_CTRL Case control section
        case_ctrl string
        %BULK_ENTRIES All bulk data entries
        bulk_entries nastran_lib.bulk.entry
    end
    methods
        % Reader
        function read_file(this,filename)
        %READ_FILE Read in raw text from a bulk file
            
            % Ensure the filename is a string
            if ~ischar(filename)
                error('filename must be a char');
            end
            
            % Check that the file exists
            if ~exist(filename,'file')
                error('File not found: %s\n',filename);
            end
            
            % Read lines from the file
            file_lines = this.get_lines(filename);
            file_lines = splitlines(string(file_lines));
            
            % Remove empty lines and comments
            temp = strip(file_lines);
            blank_lines = (temp == "");
            comment_lines = strfind(temp,'$');
            comment_lines = cellfun(@(x) any(x == 1), comment_lines);
            rem_lines = blank_lines | comment_lines;
            file_lines(rem_lines) = [];
            
            % Find important keys in the file
            temp = strip(lower(file_lines));
            CEND_line    = strfind(temp,'cend');
            CEND_line    = find(cellfun(@(x) any(x == 1), CEND_line));
            BBULK_line   = strfind(temp,'begin bulk');
            BBULK_line   = find(cellfun(@(x) any(x == 1), BBULK_line));
            ENDDATA_line = strfind(temp,'enddata');
            ENDDATA_line = find(cellfun(@(x) any(x == 1), ENDDATA_line));
            
            % Flag to turn off warning
            ed_warn_flag = true;
            
            % Add file management section, only if there's a CEND card
            if ~isempty(CEND_line)
                this.file_mgmt = file_lines(1:CEND_line-1);
                this.file_mgmt = strtrim(this.file_mgmt);
                if isempty(ENDDATA_line) && ed_warn_flag
                    warning('This deck has no ENDDATA card');
                    ed_warn_flag = false;
                end
            end
            
            % Add case control, there needs to at least be a CEND card and
            % a BEGIN BULK card
            if ~isempty(CEND_line)
                if ~isempty(BBULK_line)
                    this.case_ctrl = file_lines(CEND_line+1:BBULK_line-1);
                else
                    this.case_ctrl = file_lines(CEND_line+1:end);
                end
                % Remove comments
                this.case_ctrl = cellfun(@nastran_lib.bulk.deck.remove_comment,...
                                         this.case_ctrl,...
                                         'UniformOutput',false);
                this.case_ctrl = strtrim(this.case_ctrl);
                if isempty(ENDDATA_line) && ed_warn_flag
                    warning('This deck has no ENDDATA card');
                    ed_warn_flag = false;
                end
            end
            
            % Add bulk data
            bulk_str = [];
            if ~isempty(BBULK_line)...
            || (isempty(BBULK_line) && isempty(CEND_line))
                if ~isempty(ENDDATA_line)
                    bulk_str = file_lines(BBULK_line+1:ENDDATA_line-1);
                elseif ~isempty(BBULK_line)
                    bulk_str = file_lines(BBULK_line+1:end);
                else
                    bulk_str = file_lines(1:end);
                end
                % Remove comments
                bulk_str = cellfun(@nastran_lib.bulk.deck.remove_comment,...
                                   bulk_str,...
                                   'UniformOutput',false);
                bulk_str = deblank(bulk_str);
                if isempty(ENDDATA_line) && ~isempty(BBULK_line) && ed_warn_flag
                    warning('This deck has no ENDDATA card');
                end
                if ~isempty(ENDDATA_line) && isempty(BBULK_line) && isempty(CEND_line)
                    warning('This deck appears to contain only bulk data but also has an ENDDATA card');
                end
            end
            if ~isempty(bulk_str)
                this.bulk_entries = nastran_lib.bulk.deck.parse_bulk(bulk_str);
            end
            
        end
        % Writer
        function write_file(this,filename)
        %WRITE_FILE Write bulk entries to a file
            
            % Ensure the filename is a string
            if ~ischar(filename)
                error('filename must be a char');
            end
            
            % Open the file for writing
            file_id = fopen(filename,'W');
            if file_id == -1
                error('Could not open for writing: %s\n',filename);
            end
            
            % Loop to write each entry
            for i = 1:numel(this.bulk_entries)
                this.bulk_entries(i).write(file_id);
            end
            
            % Close the file
            fclose(file_id);
            
        end
        % Helpers
        function index = get_entry_index(this,entry_name)
            index = strcmpi(entry_name,{this.bulk_entries.name});
        end
    end
    methods (Static, Access = private)
        function file_lines = get_lines(filename)
        %GET_LINES Get all lines of text from a NASTRAN deck, while
        %observing NASTRAN's INCLUDE logic
            
            % Default output
            file_lines = '';
            
            % Open the file
            file_id = fopen(filename);
            
            % Loop until EOF reached
            while true
                
                % Get the current line
                curr_line = fgetl(file_id);
                curr_line = sprintf('%s\n',curr_line);
                
                % If it's a continuation entry, append the contents of the
                % referenced file recursively
                if numel(curr_line) > 8 && contains(lower(curr_line(1:8)),'include')
                    include_filename = strrep(curr_line(9:end-1),'''','');
                    include_filename = nastran_lib.bulk.deck.resolve_fullpath(filename,include_filename);
                    file_lines = [file_lines, nastran_lib.bulk.deck.get_lines(include_filename)]; %#ok<AGROW>
                    
                % Otherwise append to this
                else
                    file_lines = [file_lines, curr_line]; %#ok<AGROW>
                end
                
                % Stop when EOF is reached
                if feof(file_id)
                    break
                end
                
            end
            
            % Close before exiting the function
            fclose(file_id);
            
        end
        function full_path = resolve_fullpath(parent_filename,child_filename)
        %RESOLVE_FULLPATH Get the full path for an include statement
            
            % If the filename has a : or starts with \\ then it's
            % already a full path
            if contains(child_filename,':') ...
            || (numel(child_filename) >= 2 && strcmp(child_filename(1:2),'\\'))
                full_path = child_filename;
            else
                master_folder = fileparts(parent_filename);
                master_folder = sprintf('%s%s',master_folder,filesep);
                full_path = nastran_lib.references.GetFullPath(sprintf('%s%s',master_folder,child_filename));
            end
            
        end
        function text_out = remove_comment(text_in)
        %REMOVE_COMMENT Remove all comments from a line
            
            text_out = text_in;
            comment_ind = strfind(text_in,'$');
            if isempty(comment_ind)
                return
            end
            comment_ind = comment_ind(1);
            if comment_ind == 1
                text_out = '';
                return
            end
            text_out = text_in(1:comment_ind-1);
            
        end
        function bulk_entries = parse_bulk(bulk_str)
        %PARSE_BULK Parse bulk entries into "entry" objects
            
            % Loop through all lines in the BULK section and add to a cell
            % array
            bulk_cell = cell(0);
            for i = 1:length(bulk_str)
                
                % Get the current line
                curr_line = bulk_str{i};
                if isempty(curr_line)
                    continue
                end
                
                % Interpret the current line
                [curr_line,cont_key,is_cont,format] = nastran_lib.bulk.deck.interpret_bulk_entry(curr_line);
                
                % Append continuations
                if ~is_cont
                    curr_entry = nastran_lib.bulk.entry(deblank(curr_line{1}),curr_line(2:end),cont_key,format);
                    bulk_cell{end+1} = curr_entry; %#ok<AGROW>
                else
                    if isempty(cont_key)
                        entry_ind = numel(bulk_cell);
                    else
                        entry_ind = find(cellfun(@(x) strcmpi(x.continuation_key,cont_key)),1,'last');
                    end
                    bulk_cell{entry_ind}.fields = [bulk_cell{entry_ind}.fields curr_line(2:end)];
                end
                
            end
            
            % Convert to object array
            nbulk = numel(bulk_cell);
            bulk_entries = bulk_cell{1};
            if nbulk > 1
                bulk_entries(nbulk,1) = bulk_cell{end};
                if nbulk > 2
                    for i = 2:nbulk-1
                        bulk_entries(i) = bulk_cell{i};
                    end
                end
            end
            
        end
        function [bulk_cell,cont_key,is_cont,format] = interpret_bulk_entry(bulk_entry)
            
            % Default outputs
            cont_key = '';
            is_cont  = false;
            
            % Free field, easiest to parse. Note that NASTRAN enforces
            % spaces as delimeters so we need to handle that here.
            if contains(bulk_entry,',')
                bulk_entry = strrep(bulk_entry,' ',',');
                bulk_cell = strsplit(bulk_entry,',');
                format = 'free';
                return
            end
            
            % Length of the current line
            ncol = numel(bulk_entry);
            
            % Large field
            if contains(bulk_entry(1:min(8,ncol)),'*')
                field_width = 16;
                bulk_cell = cell(1,6);
                format = 'large';
                
            % Small field
            else
                field_width = 8;
                bulk_cell = cell(1,10);
                format = 'small';
            end
            
            % Get first and last entries
            bulk_cell{1} = bulk_entry(1:min(8,ncol));
            if ncol > 72
                bulk_cell{end} = bulk_entry(73:min(80,ncol));
            end
            
            % Get fields
            if ncol > 8
                field_text = bulk_entry(9:min(72,ncol));
                ncol_field = numel(field_text);
                nfield = ceil(ncol_field / field_width);
                field_text = [field_text, repmat(' ',1,(nfield * field_width) - ncol_field)];
                field_cell = cellstr(reshape(field_text,field_width,[])');
                bulk_cell(2:numel(field_cell)+1) = cellfun(@strtrim,field_cell,'UniformOutput',false);
            end
            
            % Is this a continuation?
            is_cont = any(bulk_cell{1}(1) == [' ','*','+']);

            % Find the continuation key
            if is_cont
                ncol_f1 = numel(bulk_cell{1});
                if ncol_f1 > 1
                    cont_key = strrep(bulk_cell{1}(2:min(8,ncol_f1)),' ','');
                end
            else
                ncol_fend = numel(bulk_cell{end});
                if ncol_fend > 1
                    cont_key = strrep(bulk_cell{end}(2:min(8,ncol_fend)),' ','');
                end
            end
            
            % Remove final field
            bulk_cell = bulk_cell(1:end-1);
            
        end
    end
end

