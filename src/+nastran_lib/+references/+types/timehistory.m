classdef timehistory < timeseries
%TIMEHISTORY Stores timehistory data and provides a framework for
%mathematical operations
    properties (Dependent)
        %SAMPLE_RATE Sample rate for this channel
        sample_rate;
    end
    methods
        % Getters / setters
        function sample_rate = get.sample_rate(this)
            dt = this.TimeInfo.Increment;
            if isnan(dt)
                dt = median(diff(this.Time));
            end
            sample_rate = 1/dt;
        end
        % Constructor
        function this = timehistory(varargin)
        %TIMEHISTORY Construct an instance of this class
        
            % Use base class constructor
            this = this@timeseries(varargin{:});
            
        end
        % Data processing
        function th_out = filt_butter(this,type,cutoff,order)
        %FILT_BUTTER Apply a butterworth filter
           
            % Get the filter coefficients
            nyquist_freq = this.sample_rate / 2;
            [b,a] = butter(order,cutoff/nyquist_freq,type);
            
            % Change data to match
            th_out = this;
            th_out.Data = filtfilt(b,a,this.Data);
            th_out.Name = sprintf('filter_%s',this.Name);
            
        end
        function th_out = filt_fir1(this,type,cutoff)
        %FILT_BUTTER Apply a butterworth filter
           
            % Get the filter coefficients
            nyquist_freq = this.sample_rate / 2;
            order = floor(nyquist_freq);
            if mod(order,2) == 1 
                order = order - 1;
            end
            [b,a] = fir1(order,cutoff/nyquist_freq,type);
            
            % Change data to match
            th_out = this;
            th_out.Data = filtfilt(b,a,this.Data);
            th_out.Name = sprintf('filter_%s',this.Name);
            
        end
        function srs_aa = srs(this,Q,freq)
        %SRS Acceleration SRS using Smallwood ramp invariant digital recursive
        %filtering relationship
        %
        %The highest natural frequency should be less than the Nyquist frequency.
        %Ideally, the highest natural frequency should less that one-tenth of
        %the sample rate.
        %
        %srs = [ natural frequency(Hz), peak pos(G), absolute value of peak neg(G)]
        %
        %INPUTS:
        %------
        %
        %Q    - damping quality factor
        %freq - array of natural frequencies Hz (pick your own)
        %
        %--------------------------------------------------------------------------
        %
        % ver 1.2  May 15, 2014
        %
        %By Tom Irvine
        %Adapted for this toolbox J.Scanlon 
            
            % Call SRS function
            dt = 1 / this.sample_rate;
            srs_aa = nastran_lib.references.srs(this.Data,dt,Q,freq);
            
        end
        function psd_out = psd(this,window,overlap,trace)
        %PSD Compute a PSD using Welch’s power spectral density estimate
        %
        %INPUTS:
        %------
        %
        %window  - length of each PSD segment (in seconds)
        %          default = 1s
        %overlap - percent overlap (as a ratio - e.g. enter 0.5 for 50%)
        %          default = 50%
        %trace   - either mean, maxhold, or minhold
        %          default = 'mean'
            
            % Default inputs
            if nargin < 4
                trace = 'mean';
            end
            if nargin < 3
                overlap = 0.5;
            end
            if nargin < 2
                window = [];
            end
            
            % Call PSD function
            dt = 1 / this.sample_rate;
            [pxx,freq] = nastran_lib.references.psd(this.Data,dt,window,overlap,trace);
            
            % Convert to PSD object
            time_units = this.TimeInfo.Units;
            if isempty(time_units)
                time_units = 'time';
            end
            freq_units = sprintf('(1/%s)',time_units);
            if any(strcmpi(time_units,{'s','sec','seconds','second'}))
                freq_units = 'Hz';
            end
            data_units = this.DataInfo.Units;
            if isempty(data_units)
                data_units = 'units';
            end
            power_units = sprintf('%s^2/%s',data_units,freq_units);
            name = this.Name;
            psd_out = nastran_lib.references.types.psd(freq,pxx,...
                                                       'Name',name,...
                                                       'Frequency Units',freq_units,...
                                                       'Power Units',power_units);
            
        end
        function th_out = detrend(this,varargin)
            
            th_out = this;
            th_out.Data = detrend(this.Data,varargin{:});
            th_out.Name = sprintf('detrend(%s)',this.Name);
            
        end
        function th_out = set_uniform_dt(this,dt)
            
            new_time = [this.Time(1):dt:this.Time(end)]'; %#ok<NBRAK>
            new_data = interp1(this.Time,this.Data,new_time);
            
            th_out = nastran_lib.references.types.timehistory(new_data,'Name',this.Name);
            th_out.DataInfo.Units = this.DataInfo.Units;
            th_out = setuniformtime(th_out,'Interval',dt);
            th_out.TimeInfo.StartDate = this.TimeInfo.StartDate;
            
        end
        function th_out = mov_rms(this,window,overlap)
        %MOV_RMS Calculate a moving RMS
        %
        %INPUTS:
        %------
        %
        %window  - Duration of the moving window
        %overlap - Percent overlap of windows, as a ratio 
        %
        %Adapted from FEX code written by A. Bolu Ajiboye
        %https://www.mathworks.com/matlabcentral/fileexchange/11871-signal-rms
        
            % Base signal is the data for this timeseries
            signal = this.Data;
            
            % Window length calculation
            n_window  = round(window*this.sample_rate);
            n_overlap = round(n_window*overlap);
            
            % Zeropad signal
            delta = n_window - n_overlap;
            indices = 1:delta:length(signal);
            if length(signal) - indices(end) + 1 < n_window
                indices = indices(1:find(indices+n_window-1 <= length(signal), 1, 'last'));
            end
            time = zeros(length(indices),1);
            data = zeros(length(indices),1);
            
            % Square the samples
            signal = signal.^2;
            index = 0;
            for i = indices
                index = index+1;
                % Average and take the square root of each window
                data(index) = sqrt(mean(signal(i:i+n_window-1)));
                time(index) = mean(this.Time(i:i+n_window-1));
            end
            
            % Output TH
            name_out = sprintf('mov_rms(%s)',this.Name);
            th_out = nastran_lib.references.types.timehistory(data,time,'Name',name_out);
            th_out.DataInfo.Units = this.DataInfo.Units;
            th_out.TimeInfo.StartDate = this.TimeInfo.StartDate;
            th_out.TimeInfo.Units = this.TimeInfo.Units;
            
        end
        % Acoustics
        function sound(this,scale)
        %SOUND Play a sound from this data
        
            if nargin < 2
                scale = 1;
            end
            sound(this.Data.*scale,this.sample_rate);
            
        end
        function th_out = filt_a(this)
        %FILT_A Apply a-weighting filter to the data
        %
        %Note: The A-weighting filter's coefficients 
        %are acccording to IEC 61672-1:2002 standard 
        %
        %https://www.mathworks.com/matlabcentral/fileexchange/46819-a-weighting-filter-with-matlab
            
            % translate variables
            x  = this.Data;
            fs = this.sample_rate;
            
            % determine the signal length
            xlen = length(x);
            
            % number of unique points
            NumUniquePts = ceil((xlen+1)/2);
            
            % FFT
            X = fft(x);
            
            % fft is symmetric, throw away the second half
            X = X(1:NumUniquePts);
            
            % frequency vector with NumUniquePts points
            f = (0:NumUniquePts-1)*fs/xlen;
            
            % A-weighting filter coefficients
            c1 = 12194.217^2;
            c2 = 20.598997^2;
            c3 = 107.65265^2;
            c4 = 737.86223^2;
            
            % evaluate the A-weighting filter in the frequency domain
            f   = f.^2;
            num = c1*(f.^2);
            den = (f+c2) .* sqrt((f+c3).*(f+c4)) .* (f+c1);
            A   = 1.2589*num./den;
            
            % filtering in the frequency domain
            XA = X(:).*A(:);
            
            % reconstruct the whole spectrum
            if rem(xlen, 2)                     % odd xlen excludes the Nyquist point
                XA = [XA; conj(XA(end:-1:2))];
            else                                % even xlen includes the Nyquist point
                XA = [XA; conj(XA(end-1:-1:2))];
            end
            
            % IFFT
            xA = real(ifft(XA));
            
            % represent the filtered signal in the form of the original one
            xA = reshape(xA, size(x));
            
            % Change data to match
            th_out = this;
            th_out.Data = xA;
            th_out.Name = sprintf('filterA_%s',this.Name);

        end
        function oaspl = oaspl(this,ref)
        %OASPL Calculate overall sound pressure level (in dB) of given 
        %pressure signal 'p_Pa'.
        %
        % INPUTS: 
        % ref        =  reference pressure in units matching p_Pa or simply 'air' or
        %               'water' if p_Pa is in pascals. 
        
            % Translate variables
            p_Pa = this.Data;
            
            % Check inputs
            if nargin < 2
                % Default to air
                ref = 'air';
            end
            
            % Define reference pressure
            switch ref
                case {'air','Air','AIR','gas','Gas','GAS'}
                    p_ref = 20*1e-6; % reference pressure in air is typically 20 uPa

                case {'water','Water','WATER','liquid','Liquid','LIQUID','SALTWATER','saltwater','Saltwater'}
                    p_ref = 1*1e-6; % reference pressure in water is typically 1 uPa

                otherwise
                    if ~isnumeric(ref)
                        error('ref must be numeric');
                    end
                    p_ref = ref; % reference pressure can be any user-defined value 'ref'
            end
            
            % Calculate the SPL
            oaspl = nastran_lib.references.SPL(p_Pa,p_ref);
            
        end
        function spectrogram(this)
        %SPECTROGRAM Get a spectrogram of the signal
            
            % Translate data
            x = this.Data;
            fs = this.sample_rate;
            
            % Set spectrogram parameters
            winlen = 1024;
            win = blackman(winlen, 'periodic');
            hop = round(winlen/4);
            nfft = round(2*winlen);
            
            % Get the spectrogram
            figure;
            spectrogram(x, win, winlen-hop, nfft, fs, 'yaxis');

        end
        % Arithmetic
        function th_out = uplus(this)
            
            th_out = this;
            
        end
        function th_out = uminus(this)
            
            th_out = this;
            for i = 1:numel(th_out)
                th_out(i).Data = -th_out(i).Data;
                th_out(i).Name = sprintf('-%s',th_out(i).Name);
            end
            
        end
        function th_out = plus(arg1,arg2)
            
            % Handle two timehistories
            if isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                th_out = arg1;
                for i = 1:numel(arg1)
                    [ts1,ts2] = synchronize(arg1(i),arg2(i),'intersection');
                    th_out(i) = ts1;
                    th_out(i).Data = ts1.Data + ts2.Data;
                    th_out(i).Name = sprintf('%s+%s',ts1.Name,ts2.Name);
                end
            
            % First input is a timehistory
            elseif isa(arg1,'nastran_lib.references.types.timehistory') && ~isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg1;
                for i = 1:numel(arg1)
                    th_out(i).Data = th_out(i).Data + arg2;
                    th_out(i).Name = sprintf('%s+%0.1e',th_out(i).Name,arg2);
                end
                
            % Second input is a timehistory
            elseif ~isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg2;
                for i = 1:numel(arg2)
                    th_out(i).Data = arg1 + th_out(i).Data;
                    th_out(i).Name = sprintf('%0.1e+%s',arg1,th_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function th_out = minus(arg1,arg2)
            
            % Handle two timehistories
            if isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                th_out = arg1;
                for i = 1:numel(arg1)
                    [ts1,ts2] = synchronize(arg1(i),arg2(i),'intersection');
                    th_out(i) = ts1;
                    th_out(i).Data = ts1.Data - ts2.Data;
                    th_out(i).Name = sprintf('%s-%s',ts1.Name,ts2.Name);
                end
            
            % First input is a timehistory
            elseif isa(arg1,'nastran_lib.references.types.timehistory') && ~isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg1;
                for i = 1:numel(arg1)
                    th_out(i).Data = th_out(i).Data - arg2;
                    th_out(i).Name = sprintf('%s-%0.1e',th_out(i).Name,arg2);
                end
                
            % Second input is a timehistory
            elseif ~isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg2;
                for i = 1:numel(arg2)
                    th_out(i).Data = arg1 - th_out(i).Data;
                    th_out(i).Name = sprintf('%0.1e-%s',arg1,th_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function th_out = times(arg1,arg2)
            
            % Handle two timehistories
            if isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                th_out = arg1;
                for i = 1:numel(arg1)
                    [ts1,ts2] = synchronize(arg1(i),arg2(i),'intersection');
                    th_out(i) = ts1;
                    th_out(i).Data = ts1.Data .* ts2.Data;
                    th_out(i).Name = sprintf('%s*%s',ts1.Name,ts2.Name);
                end
            
            % First input is a timehistory
            elseif isa(arg1,'nastran_lib.references.types.timehistory') && ~isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg1;
                for i = 1:numel(arg1)
                    th_out(i).Data = th_out(i).Data .* arg2;
                    th_out(i).Name = sprintf('%s*%0.1e',th_out(i).Name,arg2);
                end
                
            % Second input is a timehistory
            elseif ~isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg2;
                for i = 1:numel(arg2)
                    th_out(i).Data = arg1 .* th_out(i).Data;
                    th_out(i).Name = sprintf('%0.1e*%s',arg1,th_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function th_out = rdivide(arg1,arg2)
            
            % Handle two timehistories
            if isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                th_out = arg1;
                for i = 1:numel(arg1)
                    [ts1,ts2] = synchronize(arg1(i),arg2(i),'intersection');
                    th_out(i) = ts1;
                    th_out(i).Data = ts1.Data ./ ts2.Data;
                    th_out(i).Name = sprintf('%s/%s',ts1.Name,ts2.Name);
                end
            
            % First input is a timehistory
            elseif isa(arg1,'nastran_lib.references.types.timehistory') && ~isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg2)
                    error('arg2 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg1;
                for i = 1:numel(arg1)
                    th_out(i).Data = th_out(i).Data ./ arg2;
                    th_out(i).Name = sprintf('%s/%0.1e',th_out(i).Name,arg2);
                end
                
            % Second input is a timehistory
            elseif ~isa(arg1,'nastran_lib.references.types.timehistory') && isa(arg2,'nastran_lib.references.types.timehistory')
                
                if ~isscalar(arg1)
                    error('arg1 must be a scalar or another nastran_lib.references.types.timehistory');
                end
                th_out = arg2;
                for i = 1:numel(arg2)
                    th_out(i).Data = arg1 ./ th_out(i).Data;
                    th_out(i).Name = sprintf('%0.1e/%s',arg1,th_out(i).Name);
                end
                
            % Other situations result in error
            else
                error('Unknown input argument type');
            end
            
        end
        function th_out = power(this,arg2)
            
            if ~isscalar(arg2)
                error('second argument must be a scalar');
            end
            
            th_out = this;
            for i = 1:numel(th_out)
                th_out(i).Data = th_out(i).Data .^ arg2;
                th_out(i).Name = sprintf('%s.^%0.1e',th_out(i).Name,arg2);
            end
            
        end
        function th_out = sqrt(this)
            
            th_out = this;
            for i = 1:numel(th_out)
                th_out(i).Data = th_out(i).Data .^ 0.5;
                th_out(i).Name = sprintf('sqrt(%s)',th_out(i).Name);
            end
            
        end
        % Calculus
        function th_out = integral(this)
            
            th_out = this;
            for i = 1:numel(th_out)
                th_out(i).Data = cumtrapz(th_out(i).Time,th_out(i).Data);
                th_out(i).Name = sprintf('integral(%s)',th_out(i).Name);
            end
            
        end
    end
end

