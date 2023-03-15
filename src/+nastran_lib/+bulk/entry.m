classdef entry < handle
%ENTRY A bulk data entry
    properties
        %NAME The name of the entry
        name char
        %FIELDS All fields for the entry, stored as a cell array
        fields cell
        %CONTINUATION_KEY ID for following continuation entries
        continuation_key char
        %FORMAT The format of this entry
        format char
    end
    methods
        function this = entry(name,fields,continuation_key,format)
        %ENTRY Class constructor
            
            if nargin > 0
                this.name = name;
            end
            if nargin > 1
                this.fields = fields;
            end
            if nargin > 2
                this.continuation_key = continuation_key;
            end
            if nargin > 3
                this.format = format;
            end
            
        end
        function write(this,file_id)
        %WRITE Write this field to a file
            
            % Print name
            name_ = this.format_entry_small(this.name);
            fprintf(file_id,'%s',name_);
            
            % Stop here if no fields
            n_fields = numel(this.fields);
            if n_fields == 0
                return
            end
            
            % Print each field
            field_ind = 0;
            while field_ind < n_fields
                
                % Move in groups of 8
                for i = 1:8
                    
                    % Get index and break if end is reached
                    field_ind = field_ind + 1;
                    if field_ind > n_fields
                        break
                    end
                    
                    % Print current field
                    field_ = this.format_entry_small(this.fields{field_ind});
                    fprintf(file_id,'%s',field_);
                    
                end
                
                % Print continuation entries if there are more fields
                if n_fields > field_ind
                    cont_entry = '+       ';
                    fprintf(file_id,'%s\n%s',cont_entry,cont_entry);
                    
                % Otherwise print carriage return
                else
                    fprintf(file_id,'\n');
                end
                
            end
            
        end
    end
    methods (Static,Access = private)
        function entry_fmt = format_entry_small(entry)
        %FORMAT_ENTRY_SMALL Format an entry for NASTRAN bulk small field
        %width
            
            % Get the data type
            curr_class = class(entry);
            
            % Empty is easy
            if isempty(entry)
                entry_fmt = '        ';
                return
            end
            
            % Print method depends on type
            switch curr_class
                
                % Method for integers and logicals
                case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','logical'}
                    if entry > 99999999 || entry < -9999999
                        error('Integer field is too large to print');
                    end
                    entry_fmt = sprintf('%8d',entry);
                    
                % Method for floats... this one is pretty complicated and
                % can probably be optimized. Idea is to print as compact as
                % possible
                case {'double','single'}
                    
                    exp_rep = 0;
                    if entry == 0
                        entry_fmt = '0.000000';
                    else
                        exponent_1 = floor(log10(abs(entry)));
                        if abs(exponent_1) > 999
                            error('Floating point number is too large to print');
                        end
                        if entry > 0
                            if exponent_1 > 99
                                entry_fmt = sprintf('%0.2E',entry);
                                exp_rep = 1;
                            elseif exponent_1 > 9
                                entry_fmt = sprintf('%0.3E',entry);
                                exp_rep = 1;
                            elseif exponent_1 > 6
                                entry_fmt = sprintf('%0.4E',entry);
                                exp_rep = 2;
                            elseif exponent_1 == 6
                                entry_fmt = sprintf('%7.0f.',entry);
                            elseif exponent_1 >= 0
                                n_bd = exponent_1 + 1;
                                n_ad = 6 - exponent_1;
                                entry_fmt = sprintf(sprintf('%%%d.%df',n_bd,n_ad),entry);
                            elseif exponent_1 >= -3
                                entry_fmt = strrep(sprintf('%0.7f',entry),'0.','.');
                            elseif exponent_1 >= -9
                                entry_fmt = sprintf('%0.4E',entry);
                                exp_rep = 3;
                            elseif exponent_1 >= -99
                                entry_fmt = sprintf('%0.3E',entry);
                                exp_rep = 4;
                            else
                                entry_fmt = sprintf('%0.2E',entry);
                                exp_rep = 4;
                            end
                        else
                            if exponent_1 > 99
                                entry_fmt = sprintf('%0.1E',entry);
                                exp_rep = 1;
                            elseif exponent_1 > 9
                                entry_fmt = sprintf('%0.2E',entry);
                                exp_rep = 1;
                            elseif exponent_1 > 5
                                entry_fmt = sprintf('%0.3E',entry);
                                exp_rep = 2;
                            elseif exponent_1 == 5
                                entry_fmt = sprintf('%6.0f.',entry);
                            elseif exponent_1 >= 0
                                n_bd = exponent_1 + 1;
                                n_ad = 5 - exponent_1;
                                entry_fmt = sprintf(sprintf('%%%d.%df',n_bd,n_ad),entry);
                            elseif exponent_1 >= -3
                                entry_fmt = strrep(sprintf('%0.6f',entry),'0.','.');
                            elseif exponent_1 >= -9
                                entry_fmt = sprintf('%0.3E',entry);
                                exp_rep = 3;
                            elseif exponent_1 >= -99
                                entry_fmt = sprintf('%0.2E',entry);
                                exp_rep = 4;
                            else
                                entry_fmt = sprintf('%0.1E',entry);
                                exp_rep = 4;
                            end
                        end
                        
                        % Make sure the value ended up as expected, since
                        % rounding can cause weird issues sometimes
                        if numel(entry_fmt) ~= 8
                            entry_test = str2double(entry_fmt);
                            exponent_2 = floor(log10(abs(entry_test)));
                            if ~isequal(exponent_1,exponent_2)
                                entry_fmt = nastran_lib.bulk.entry.format_entry_small(entry_test);
                            end
                        end

                        % Remove exponent text
                        switch exp_rep
                            case 1
                                entry_fmt = strrep(entry_fmt,'E+' ,'+');
                            case 2
                                entry_fmt = strrep(entry_fmt,'E+0','+');
                            case 3
                                entry_fmt = strrep(entry_fmt,'E-0','-');
                            case 4
                                entry_fmt = strrep(entry_fmt,'E-' ,'-');
                        end
                        
                    end
                    
                % Method for character arrays
                case {'char'}
                    field_width = numel(entry);
                    if field_width > 8
                        error('Character array is too large to print');
                    end
                    entry_fmt = [entry,repmat(' ',1,8-field_width)];
                
                % Nothing else is defined so far
                otherwise
                    error('Type %s is not coded\n',curr_class);
            end
            
            % Double-check code output
            if numel(entry_fmt) ~= 8
                error('Code is broken, contact developer');
            end
            
        end
    end
end

