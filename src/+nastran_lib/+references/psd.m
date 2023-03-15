function [pxx,freq] = psd(data,dt,window,overlap,trace)
%PSD Compute a PSD using Welch’s power spectral density estimate
%
%INPUTS:
%------
%
%data    - signal to be processed
%dt      - time step
%window  - length of each PSD segment (in seconds)
%          default = longest possible segments to obtain as close to but
%          not exceed 8 segments with 50% overlap
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
    
    % Trace verification
    if ~any(strcmp(trace,{'mean','maxhold','minhold'}))
        error('trace for PSD is invalid');
    end
    
    % Convert to info that's usable by the pwelch function
    window_nsamples  = floor(window/dt);
    overlap_nsamples = floor(window_nsamples * overlap);

    % Run through the pwelch function
    [pxx,freq] = pwelch(data,window_nsamples,overlap_nsamples,[],1/dt,trace);

end

