classdef psd
%PSD Container for a Power Spectral Density
    properties
        %NAME The name of the dataset
        Name(1,:) char
    end
    properties (SetAccess = private)
        %DATA Power spectral density vector
        Data(:,1) double
        %FREQUENCY Frequency vector
        Frequency(:,1) double
        %FREQUENCYINFO Info about the frequency vector
        FrequencyInfo struct
        %POWERINFO Info about the Power spectral density vector
        DataInfo struct
    end
    methods
        % Getters and setters
        function this = set.Frequency(this,Frequency)
            if isrow(Frequency)
                Frequency = Frequency';
            end
            if any(isnan(Frequency))
                error('Frequency cannot be NaN');
            end
            this.Frequency = Frequency;
        end
        function this = set.Data(this,Data)
            if isrow(Data)
                Data = Data';
            end
            if any(isnan(Data)) || any(Data<0)
                error('Power must be positive and numeric');
            end
            this.Data = Data;
        end
        % Constructor
        function this = psd(Frequency,Data,varargin)
        %PSD Class constructor
        %
        %INPUTS:
        %------
        %
        %Frequency       - Frequency vector
        %Data            - Power spectral density vector
        %
        %OPTIONS:
        %-------
        %Enter these as 'option Name',value
        %
        %Name            - Name of psd
        %Frequency Units - Units for Frequency
        %Data Units      - Units for Power
        
            % Defaults
            this.Name = 'unnamed';
            this.FrequencyInfo = struct;
            this.FrequencyInfo.Units = 'freq';
            this.DataInfo = struct;
            this.DataInfo.Units = 'unnamed';
            
            % Parse varargs
            for i = 1:2:numel(varargin)
                arg_title = varargin{i};
                arg_value = varargin{i+1};
                switch lower(arg_title)
                    case 'name'
                        this.Name = arg_value;
                    case 'frequency units'
                        this.FrequencyInfo.Units = arg_value;
                    case 'power units'
                        this.DataInfo.Units = arg_value;
                    otherwise
                        error('%s is an unknown argument',arg_title);
                end
            end
            
            % Set class data
            if nargin > 0
                this = this.set_Frequency_Data(Frequency,Data);
            else
                this = this.set_Frequency_Data([],[]);
            end
            
        end
        % Metrics
        function rms_out = RMS(this)
        %RMS Calculate the RMS of the PSD via numerical integration
        
            rms_out = sqrt(trapz(this.Frequency,abs(this.Data)));
            
        end
        % Other methods
        function varargout = plot(this,varargin)
        %PLOT Plot the PSD
            
            % Plot the PSD figure
            if isreal(this.Data)
                loglog(this.Frequency,abs(this.Data,varargin{:}));
                f = gcf;
                ax = gca;
                title(ax,sprintf('PSD Plot: %s',this.Name));
                xlabel(ax,sprintf('Frequency (%s)',this.FrequencyInfo.Units));
                ylabel(ax,sprintf('%s (%s)',this.Name,this.DataInfo.Units));
                grid(ax,'on');
            else
                subplot(2,1,1); loglog(this.Frequency,abs(this.Data,varargin{:}));
                f = gcf;
                ax1 = gca;
                title(ax1,sprintf('PSD Plot: %s',this.Name));
                xlabel(ax1,sprintf('Frequency (%s)',this.FrequencyInfo.Units));
                ylabel(ax1,sprintf('%s (%s)',this.Name,this.DataInfo.Units));
                grid(ax1,'on');
                subplot(2,1,2); semilogx(this.Frequency,angle(this.Data,varargin{:}));
                ax2 = gca;
                xlabel(ax2,sprintf('Frequency (%s)',this.FrequencyInfo.Units));
                ylabel(ax2,'Phase (rad)');
                grid(ax2,'on');
            end
            
            % Handle output arguments
            if nargout == 1
                varargout{1} = f;
            end
            if nargout > 1
                error('Max of 1 output argument');
            end
            
        end
        function th_out = timehistory(this,duration,nseg)
        %TIMEHISTORY Synthesize a time history to satisfy a PSD
        %
        %INPUTS:
        %------
        %
        %duration - (optional) Duration of PSD, default 1 unit of time
        %nseg     - (optional) number of segments, user will be prompted to choose 
        %           a value if this argument in not provided
        %
        %OUTPUTS:
        %-------
        %
        %th_out   - timehistory object
        %
        %See Also nastran_lib.references.psd2time
            
            % Defaults
            if nargin < 2
                duration = 1;
            end
            if nargin < 3
                nseg = [];
            end
            
            % Perform a psd2time, adapted from Tom Irvine's codebase
            [amp_out,time_out] = nastran_lib.references.psd2time(this.Frequency,...
                                                                 abs(this.Data),...
                                                                 duration,...
                                                                 nseg);
            
            % Create output time history
            th_out = nastran_lib.references.types.timehistory(amp_out,time_out,'Name',this.Name);
            
            % Parse time units
            time_units = '';
            if strcmpi(this.FrequencyInfo.Units,'Hz')
                time_units = 'seconds';
            end
            if numel(this.FrequencyInfo.Units) > 2
                if strcmpi(this.FrequencyInfo.Units(1:2),'1/')
                    time_units = this.FrequencyInfo.Units(3:end);
                end
            end
            th_out.TimeInfo.Units = time_units;
            
            % Parse data units
            data_units = '';
            carrot_ind = strfind(this.DataInfo.Units,'^');
            if ~isempty(carrot_ind)
                carrot_ind = carrot_ind(1);
                if carrot_ind > 1
                    data_units = this.DataInfo.Units(1:carrot_ind-1);
                    data_units = strrep(data_units,'(','');
                    data_units = strrep(data_units,')','');
                end
            end
            th_out.DataInfo.Units = data_units;
                                                            
        end
        function this = set_Frequency_Data(this,Frequency,Data)
        %SET_FREQUENCY_POWER Set class properties with error checking
            
            % Input checks
            if ~(numel(Frequency) == numel(Data))
                error('Frequency and Data vectors must contain the same number of elements');
            end
            
            % Set class data
            this.Frequency = Frequency;
            this.Data      = Data;
            
        end
        function [psd1,psd2] = set_common_domain(this,arg2)
        %SET_COMMON_DOMAIN Put both PSDs on a common frequency vector
            
            % Error check
            if ~isa(arg2,'nastran_lib.references.types.psd')
                error('Both inputs must be of type nastran_lib.references.types.psd');
            end
            
            % Get frequency vector for each input
            freq1 = this.Frequency;
            freq2 = arg2.Frequency;
            
            % Default
            psd1 = this;
            psd2 = arg2;
            
            % Easy case
            if isequal(freq1,freq2)
                return
            end
            
            % Find common frequency vector
            freq_common = unique([freq1;freq2]);
            
            % Use the first input as the basis for interpolation
            power1 = nastran_lib.references.interp1_log(freq1,this.Data,freq_common);
            power2 = nastran_lib.references.interp1_log(freq2,arg2.Data,freq_common);
            
            % Find valid data points
            valid_inds = ~(isnan(power1) | isnan(power2));
            freq_common = freq_common(valid_inds);
            power1 = power1(valid_inds);
            power2 = power2(valid_inds);
            
            % Update output arguments
            psd1.set_Frequency_Data(freq_common,power1);
            psd2.set_Frequency_Data(freq_common,power2);
            
        end
        % Arithmetic
        function psd_out = uplus(this)
            
            psd_out = this;
            
        end
        function psd_out = uminus(this)
            
            psd_out = this;
            for i = 1:numel(psd_out)
                psd_out(i) = psd_out(i).set_Frequency_Data(psd_out(i).Frequency,-psd_out(i).Data);
                psd_out(i).Name = sprintf('-%s',psd_out(i).Name);
            end
            
        end
        function psd_out = plus(arg1,arg2)
            
            % Handle two PSDs
            if isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    [psd1,psd2] = set_common_domain(arg1(i),arg2(i));
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd1.Frequency,psd1.Data+psd2.Data);
                    psd_out(i).Name = sprintf('%s+%s',arg1(i).Name,arg2(i).Name);
                end
                
            % First input is a PSD
            elseif isa(arg1,'nastran_lib.references.types.psd') && ~isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,psd_out.Data+arg2);
                    psd_out(i).Name = sprintf('%s+%0.1e',psd_out(i).Name,arg2);
                end
                
            % Second input is a PSD
            elseif ~isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg2;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,arg1+psd_out.Data);
                    psd_out(i).Name = sprintf('%0.1e+%s',arg1,psd_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function psd_out = minus(arg1,arg2)
            
            % Handle two PSDs
            if isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    [psd1,psd2] = set_common_domain(arg1(i),arg2(i));
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd1.Frequency,psd1.Data-psd2.Data);
                    psd_out(i).Name = sprintf('%s-%s',arg1(i).Name,arg2(i).Name);
                end
                
            % First input is a PSD
            elseif isa(arg1,'nastran_lib.references.types.psd') && ~isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,psd_out.Data-arg2);
                    psd_out(i).Name = sprintf('%s-%0.1e',psd_out(i).Name,arg2);
                end
                
            % Second input is a PSD
            elseif ~isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg2;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,arg1-psd_out.Data);
                    psd_out(i).Name = sprintf('%0.1e-%s',arg1,psd_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function psd_out = times(arg1,arg2)
            
            % Handle two PSDs
            if isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    [psd1,psd2] = set_common_domain(arg1(i),arg2(i));
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd1.Frequency,psd1.Data.*psd2.Data);
                    psd_out(i).Name = sprintf('%s*%s',arg1(i).Name,arg2(i).Name);
                end
                
            % First input is a PSD
            elseif isa(arg1,'nastran_lib.references.types.psd') && ~isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,psd_out.Data.*arg2);
                    psd_out(i).Name = sprintf('%s*%0.1e',psd_out(i).Name,arg2);
                end
                
            % Second input is a PSD
            elseif ~isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg2;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,arg1.*psd_out.Data);
                    psd_out(i).Name = sprintf('%0.1e*%s',arg1,psd_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function psd_out = rdivide(arg1,arg2)
            
            % Handle two PSDs
            if isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    [psd1,psd2] = set_common_domain(arg1(i),arg2(i));
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd1.Frequency,psd1.Data./psd2.Data);
                    psd_out(i).Name = sprintf('%s/%s',arg1(i).Name,arg2(i).Name);
                end
                
            % First input is a PSD
            elseif isa(arg1,'nastran_lib.references.types.psd') && ~isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg1;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,psd_out.Data./arg2);
                    psd_out(i).Name = sprintf('%s/%0.1e',psd_out(i).Name,arg2);
                end
                
            % Second input is a PSD
            elseif ~isa(arg1,'nastran_lib.references.types.psd') && isa(arg2,'nastran_lib.references.types.psd')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.psd');
                end
                psd_out = arg2;
                for i = 1:numel(psd_out)
                    psd_out(i) = psd_out(i).set_Frequency_Data(psd_out.Frequency,arg1./psd_out.Data);
                    psd_out(i).Name = sprintf('%0.1e/%s',arg1,psd_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function psd_out = power(this,arg2)
            
            if ~isscalar(arg2)
                error('second argument must be a scalar');
            end
            
            psd_out = this;
            for i = 1:numel(psd_out)
                psd_out(i) = psd_out(i).set_Frequency_Data(psd_out(i).Frequency,psd_out(i).Data .^ arg2);
                psd_out(i).Name = sprintf('%s.^%0.1e',psd_out(i).Name,arg2);
            end
            
        end
        function psd_out = sqrt(this)
            
            psd_out = this;
            for i = 1:numel(psd_out)
                psd_out(i) = psd_out(i).set_Frequency_Data(psd_out(i).Frequency,sqrt(psd_out(i).Data));
                psd_out(i).Name = sprintf('sqrt(%s)',psd_out(i).Name);
            end
            
        end
    end
end

