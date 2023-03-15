function varargout = psd2time(freq_in,psd_in,duration,nseg)
%PSD2TIME Synthesize a time history to satisfy a PSD
%
%This program synthesizes a time history to satisfy a PSD.
%
%INPUTS:
%------
%
%freq_in  - frequency vector input
%psd_in   - PSD vector input
%duration - (optional) Duration of PSD, default 1 unit of time
%nseg     - (optional) number of segments, user will be prompted to choose 
%           a value if this argument in not provided
%
%OUTPUTS:
%-------
%
%amp_out  - amplitude vector output
%time_out - time vector output
%
%
%Adapted from:
%
%http://www.vibrationdata.com/matlab2/generic_psd_syn.zip
%By Tom Irvine   Email: tom@vibrationdata.com 
%version 1.0     August 15, 2013
    
    % Defaults
    if nargin < 3
        duration = 1;
    end
    if nargin < 4
        nseg = [];
    end

    % Check for zero frequency
    if freq_in(1) < 1.0e-09
        freq_in(1) = 10^(floor(-0.1+log10(freq_in(2))));
    end
    
    % Check for zero amplitude
    if psd_in(1) < 1.0e-30
        noct = log(freq_in(2)/freq_in(1)) / log(2);
        psd_in(1) = (noct/4) * psd_in(2); % 6 db/octave 
    end
    
    % New freq and amplitude without zero freqs
    nfreq_in = numel(freq_in);
    freq     = zeros(nfreq_in+1,1);
    amp      = zeros(nfreq_in+1,1);
    k = 1;
    for i = 1:nfreq_in
        if (freq_in(i) > 0)
            amp(k)  = psd_in(i);
            freq(k) = freq_in(i);
            k = k+1;
        end
    end
    
    % Issue a warning if duration too low
    min_dur = 10/freq(1);
    if duration < min_dur
        warning(['Recommended duration > %9.6g sec \n',...
                 'If recommended duration is too long, then increase initial ',...
                 'specification frequency.\n'],min_dur);
    end
    
    % Get slope and GRMS
    freq(nfreq_in+1) = freq(nfreq_in)*2^(1/48);
    amp(nfreq_in+1)  = amp(nfreq_in);
    grms = calculate_grms(freq,amp,nfreq_in);
    
    % Sample rate should be 2*10x Nyquist frequency
    sr = 20*max(freq);
    dt = 1/sr;
    spec_grms = grms;
    
    % Fix vector orientation
    freq = fix_size(freq);
    amp  = fix_size(amp);
    
    % Generate White Noise 
    [np,dt,df,white_noise] = PSD_syn_white_noise(duration,dt);
    
    % Interpolated PSD spec
    spec = interpolate_PSD_spec(np,freq,amp,df);
    
    % Core synthesis
    mmm       = round(np/2);
    fmax      = max(freq);
    nsegments = 1;
    sq_spec   = sqrt(spec);
    [psd_th,nL] = PSD_syn_FFT_core(nsegments,mmm,np,fmax,df,sq_spec,white_noise);
    
    % Scale the result
    [TT,psd_th] = PSD_syn_scale_time_history(psd_th,grms,nL,dt,duration);
    
    % Verification
    [freq,amp,NW,df,mr_choice,h_choice,den,mH] = PSD_syn_verify(TT,psd_th,spec_grms,dt,nseg);
    
    % Correction
    nnt = 3;
    [amp_out,full_,time_out] = generic_psd_syn_correction(nnt,amp,spec_grms,mH,NW,freq_in,psd_in,df,mr_choice,h_choice,den,freq,TT);
    
    % Fix sizes
    time_out = fix_size(time_out);
    amp_out = fix_size(amp_out);
    full_ = fix_size(full_);
    
    % Output options
    switch nargout
        case 0
            % If no outputs, plot the result
            
            % Time history
            figure;
            plot(time_out,amp_out);
            title('Time History Generated');
            xlabel('Time');
            ylabel('Amplitude');
            grid on;
            
            % PSD comparison
            figure;
            loglog(freq,full_,'b');
            hold on;
            loglog(freq_in,psd_in,'r');
            loglog(freq_in,psd_in*sqrt(2),'k');
            loglog(freq_in,psd_in/sqrt(2),'k');
            legend ('Synthesis','Specification','+ 1.5 dB tol','- 1.5 dB tol');
            title('PSD Comparison');
            xlabel('Freq');
            ylabel('PSD');
            grid on;
            xlim([min(freq_in)/10,max(freq_in*10)]);
            ylim([min(psd_in)/10,max(psd_in*10)]);
       
        case 1
            varargout{1} = amp_out;
        case 2
            varargout{1} = amp_out;
            varargout{2} = time_out;
        otherwise
            error('User can only specify 0-2 output arguments');
    end
    
end
function grms = calculate_grms(freq,amp,nsz)
    
    % Basic info
    n  = length(amp);
    
    % Get slope
    slope = zeros(n-1,1);
    for i = 1:(n-1)
        if freq(i) > 1.0e-12
            slope(i) = log(amp(i+1)/amp(i)) / log(freq(i+1)/freq(i));
        else
            slope(i) = 0;
        end
    end
    
    % Get ra
    ra = 0;
    for i = 1:nsz-1
        if slope(i) < -1.0001 ||  slope(i) > -0.9999
            ra = ra + (amp(i+1)*freq(i+1)-amp(i)*freq(i))/(slope(i)+1);
        else  
            ra = ra + amp(i)*freq(i)*log(freq(i+1)/freq(i));
        end
    end
    
    % Get GRMS
    grms = sqrt(ra);
    
end
function [np,dt,df,white_noise] = PSD_syn_white_noise(tmax,dt)
    
    np   = ceil(tmax/dt);
    noct = ceil(log(np/2)/log(2));
    np   = 2^noct;
    dt   = tmax/(np-1);
    df   = 1/(np*dt);

    white_noise = randn(np,1);
    white_noise = fix_size(white_noise);
    white_noise = white_noise-mean(white_noise);

end
function spec = interpolate_PSD_spec(np,freq,amp,df)

    fft_freq = linspace(0,(np-1)*df,np);
    fft_freq = fft_freq';
    
    ls = length(freq);
    x = zeros((ls+4),1);
    y = (1.0e-30)*ones((ls+4),1);
    x(3:(ls+2)) = freq;
    y(3:(ls+2)) = amp;
    
    x(1) = 1.0e-90;
    x(2) = x(3)*0.999;
    x(ls+3) = x(ls+2)*1.00001;
    x(ls+4) = x(ls+3)*10000;
    
    x = log10(x);
    Y = log10(y);
    
    fft_freq(1)=1.0e-80;
    
    xi = log10(fft_freq);
    
    yi = interp1(x,Y,xi);
    
    spec = 10.^yi;
    
    for i=1:length(yi)
        if ~(spec(i) >= 0 && spec(i) <= 1.0e+20)
            spec(i) = 0;
        end
    end
    
end
function [amp,full_,tim] = generic_psd_syn_correction(nnt,amp,spec_grms,mH,NW,freq_spec,amp_spec,df,mr_choice,h_choice,den,freq,TT)

    tpi   = 2*pi;
    delta = freq_spec(1)/30;
    n     = length(amp);
    tim   = double(TT(1:n));
    tim   = fix_size(tim);
    mmm   = 2^fix(log(n/NW)/log(2));
    for kvn = 1:nnt

        ratio = (spec_grms/std(amp));
        amp   = amp*ratio;    
        full_ = zeros(mH,1);
        nov   = 0;
        for ijk = 1:(2*NW-1)
            
            amp_seg = zeros(mmm,1);
            amp_seg(1:mmm) = amp((1+ nov):(mmm+ nov));  
            nov = nov+fix(mmm/2);
            
            mag_seg = FFT_core(amp_seg,mmm,mH,mr_choice,h_choice);
            mag_seg = fix_size(mag_seg);
            
            full_ = full_+mag_seg(1:mH);
            
        end
        full_ = full_/den;
        
        nlen = length(full_);
        for i = 1:nlen-1
            if freq(i) <= freq_spec(1) && freq_spec(1) <= freq(i+1)
                x    = freq_spec(1) - freq(i);
                c2   = x / (freq(i+1) - freq(i));
                c1   = 1 - c2;
                psd1 = c1*full_(i) + c2*full_(i+1);
                break
            end
        end  
        
        if psd1 < amp_spec(1) && kvn < nnt
            
            ca  = sqrt(2)*sqrt(amp_spec(1)*df-psd1*df);
            pha = tpi*rand;
            fff = freq_spec(1)+(-0.5+rand)*delta;
            
            switch kvn
                case 1
                    fff = freq_spec(1);
                case 2
                    fff = freq_spec(1)+delta/2;
                case 3
                    fff=freq_spec(1)-delta/2;
            end
            
            amp = amp + ca*sin(tpi*fff*tim+pha);
            
        end
    end
    
    amp = amp*(spec_grms/std(amp));
    amp = fix_size(amp);
    
end
function [Y_real,nL] = PSD_syn_FFT_core(nsegments,mmm,np,fmax,df,sq_spec,white_noise)

    for i=1:nsegments
        
        YF = fft(white_noise,np);
        Y  = sq_spec.*YF;
        
        % Make symmetric
        for j = 1:mmm
            fff = (j-1)*df;
            if fff >= fmax
                Y(j) = 0;
            end
        end
        
        for j = 2:np             
            Y(j) = complex(real(Y(j)),-imag(Y(j)));
        end

        Y_real = real(ifft(Y));
           
    end
    
    nL = length(Y_real);

end
function [TT,psd_th] = PSD_syn_scale_time_history(psd_th,grms,nL,dt,tmax)

    TT = linspace(0,nL*dt,nL);
    
    if max(TT) > tmax
        
        nnn        = round(tmax/dt);
        new_TT     = TT(1:nnn);
        new_psd_th = psd_th(1:nnn);
        
        TT     = new_TT;
        psd_th = new_psd_th;
        
    end
    
    psd_th = detrend(psd_th);
    
    % scale th
    mu      = mean(psd_th);
    stddev  = std(psd_th);
    grmsout = sqrt(mu^2+stddev^2);
    scale   = grms/grmsout;
    psd_th  = psd_th*scale;
    
end
function [freq,amp,NW,df,mr_choice,h_choice,den,mH] = PSD_syn_verify(TT,psd_th,spec_grms,dt,nseg)

    % mean removal yes
    mr_choice = 1;
    
    % Hanning Window
    h_choice  = 2;
    
    amp = double(psd_th);
    n   = length(amp);
    
    tim = double(TT(1:n));
    tim = tim';
    
    dtmin = min(diff(tim));
    dtmax = max(diff(tim));
    
    if (dtmax-dtmin)/dt > 0.01
        warning('Time step is not constant.');
    end
    
    [df,mmm,NW] = PSD_advise(dt,n,nseg);
    
    % begin overlap
    mH = (mmm/2)-1;
    
    amp = amp*(spec_grms/std(amp));
    
    fmax = (mH-1)*df;
    freq = linspace(0,fmax,mH);
    
    den = df*(2*NW-1);
    
end
function mag_seg = FFT_core(amp_seg,mmm,mH,mr_choice,h_choice)

    mu = mean(amp_seg);
    if mr_choice == 1
       amp_seg = amp_seg - mu;
    end
    
    if h_choice == 2
        alpha = linspace(0,2*pi,mmm);
        H     = 0.5*(1-cos(alpha));
        ae    = sqrt(8./3.);
        
        sz = size(H);
        if sz(2) > sz(1)
            H = H';
        end
        
        amp_seg = ae*amp_seg.*H;
        
    end
	
    Y = fft(amp_seg,mmm);
	
    Ymag       = abs(Y);
    mag_seg    = [];
    mag_seg(1) = sqrt(Ymag(1))/mmm;
    mag_seg(1) = mag_seg(1)^2; 
	
    mag_seg(2:mH) = 2*Ymag(2:mH)/mmm;
    mag_seg(2:mH) = (mag_seg(2:mH).^2)/2.;
    
end
function[df,mmm,nseg]=PSD_advise(dt,n,nseg)
    
    NC    = 0;
    ss    = zeros(1,1000);
    seg   = zeros(1,1000);
    i_seg = zeros(1,1000);
    for i = 1:1000
        nmp = 2^(i-1);
        if nmp <= n
            ss(i)    = 2^(i-1);
            seg(i)   = n/ss(i);
            i_seg(i) = fix(seg(i));
            NC = NC+1;
        else
            break;
        end
    end
    
    % Display info about segments if none entered
    if isempty(nseg)
        fprintf(1,'\n Number of    Samples per    Time per                        \n');
        fprintf(1,  ' Segments     Segment        Segment         df         dof  \n');
    end
    
    for i = 1:NC
        j = NC+1-i;
        if j > 0
            if i_seg(j) > 0
                tseg = dt*ss(j);
                ddf  = 1./tseg;
                % Display info about segments if none entered
                if isempty(nseg)
                    fprintf(1,' %8d  %8d      %11.5f    %9.4f   %8d\n',i_seg(j),ss(j),tseg,ddf,2*i_seg(j));
                end
            end
        end
        if i == 12
            break;
        end
    end
    
    looping = true;
    while looping
        
        % Prompt for segments if none entered
        if isempty(nseg)
            fprintf(1,'\n');
            nseg = input(' Choose the Number of Segments from the first column:  ');
            fprintf(1,'\n');
        end
        
        for j = 1:length(i_seg)
            if nseg == i_seg(j)
                looping = false;
            end
        end
        
    end
	
    mmm = 2^fix(log(n/nseg)/log(2));
    df  = 1./(mmm*dt);
    
end
function a = fix_size(a)
    
    sz = size(a);
    if sz(2) > sz(1)
        a = a';
    end
    
end
